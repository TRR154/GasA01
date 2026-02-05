/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of the code gasa01, which contains code              */
/*    for Project A01                                                        */
/*    "Global Methods for Stationary and Instationary Gastransport"          */
/*    of the SFB/TRR 154 "Mathematical Modelling, Simulation and             */
/*    Optimization on the Example of Gas Networks".                          */
/*                                                                           */
/*    Copyright (C) 2015-     Discrete Optimization Group, TU Darmstadt      */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_pipe_ode.c
 * @brief  constraint handler for pipe_ode constraints
 * @author Oliver Habeck
 * @author Alexandra Stern
 */

/* #define SCIP_DEBUG */
/* #define UNDERESTIMATOR_DEBUG */
/* #define UNDERESTIMATOR_DEBUG_MIN_OUTPUT */
/* #define OVERESTIMATOR_DEBUG */
/* #define OVERESTIMATOR_DEBUG_MIN_OUTPUT */
/* #define STEPSIZE_DEBUG */
/* #define FLOW_TIGHTENING_DEBUG */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <scip/misc.h>
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
#include <scip/scip_expr.h>
#endif
#include "scip/scipdefplugins.h"       /* needed for the secondary SCIP instance */
#include "scip/interrupt.h"
#include "elements_gas.h"
#include "ode_functions_gas.h"
#include "cons_pipe_ode.h"
#include "comparisons_gas.h"
#include "data_structs_cons_pipe_ode.h"

#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
#include <scip/type_expr.h>
#include <scip/pub_expr.h>
#endif

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "pipe_ode"
#define CONSHDLR_DESC          "constraint for odes arising from pipes in a gas network"
#define CONSHDLR_ENFOPRIORITY        -100 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY       -100 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ            100 /**< frequency for using all instead of only the useful constraints in separation,
                                           *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS           TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
/* TODO: remove properties which are never used because the corresponding routines are not supported */
#define CONSHDLR_SEPAPRIORITY        -100 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ               1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA          FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ               1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP          FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING         SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */

#define CONSHDLR_MAXPREROUNDS          -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_PRESOLTIMING     (SCIP_PRESOLTIMING_EXHAUSTIVE | SCIP_PRESOLTIMING_FINAL)
/* #define CONSHDLR_PRESOLTIMING     SCIP_PRESOLTIMING_MAX */

/* default parameter values */
#define DEFAULT_WITH_FLOWTIGHTENING  TRUE /**< if propagation based flowTightening should be used */
#define DEFAULT_WITH_OBBT            TRUE /**< if OBBT should be used */
#define DEFAULT_OBBT_MINGAP          10.0 /**< only perform obbt for flow variable if bounds differ by at least obbt_mingap */
#define DEFAULT_NLP_REPRESENTATION   TRUE /**< whether the pipes should be represented in the NLP */
#define DEFAULT_NLP_CONVEX           TRUE /**< whether the pipes are represented by a convex relaxation */

/* event handler properties */
#define EVENTHDLR_NAME         "ODE"
#define EVENTHDLR_DESC         "bound change event handler for pipe_ode constraints"


#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
/* expression handler properties */
#define EXPRHDLR_NAME         "pipe"
#define EXPRHDLR_DESC         "pipe ODE solution operator"
#define EXPRHDLR_PRECEDENCE   70000
#endif


/*
 *  Data structures
 */

/* typedef conshdlrdata and consdata see data_structs_cons_pipe_ode.h */

/*
 * Local methods
 */

/** enumerate the different types the solution can be infeasible */
enum
{
   FEASIBLE = -1,                                 /**< solution is feasible */
   MACHBOUND = 1,                                 /**< machbound is not fulfilled at the output */
   PIN_TOO_BIG = 2,                               /**< the input pressure is too big */
   PIN_TOO_LOW = 3                                /**< the input pressure is too low */
};

/** enumerate reasons for bound change of pressure and flow variables */
enum
{
   FPROP_UPPER = 1,               /**< Forward propagation fo upper pressure bound */
   FPROP_LOWER = 2,               /**< Forward propagation of lower pressure bound */
   BPROP_UPPER = 3,               /**< Backward propagation of upper pressure bound */
   BPROP_LOWER = 4,               /**< Backward propagation of lower pressure bound */
   PIN_GT_POUT = 5,               /**< Lower bound on p_in is greater than upper bound on p_out, ie, no negative flow */
   PIN_LT_POUT = 6,               /**< Upper Bound on p_in is less than lower bound on p_out, ie, no positive flow */
   MACH_BOUND = 7,                /**< One vertex of the pressure-flow-box does not fulfill the "mach"bound */
   POS_LB = 8,                    /**< change of lb of positive flow binvar */
   POS_UB = 9,                    /**< change of ub of positive flow binvar */
   NEG_LB = 10,                   /**< change of lb of negative flow binvar */
   NEG_UB = 11                    /**< change of ub of negative flow binvar */
};

/** function which enforces the "mach-bound" on the pressure-flow box. */
static
SCIP_RETCODE machBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of the constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_RESULT*          result,             /**< result code */
   int*                  nchgbds             /**< total number of variable bounds tightened by this presolver */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool infeasible;
   SCIP_Bool tightened;
   SCIP_Real p_min;
   SCIP_Real q_max;
   SCIP_Real q_var_lb;
   SCIP_Real q_var_ub;
   SCIP_Real numerator;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( cons != NULL );
   assert( result != NULL );
   assert( nchgbds != NULL );

   consdata = SCIPconsGetData(cons);

   q_var_lb = SCIPcomputeVarLbLocal(scip, consdata->q_var);
   q_var_ub = SCIPcomputeVarUbLocal(scip, consdata->q_var);

   numerator = (SCIP_Real) conshdlrdata->machbound;

   /* Lower Bounds on pressure values */
   if ( SCIPisFeasPositive(scip, q_var_lb) || SCIPisFeasNegative(scip, q_var_ub) )
   {
      p_min = (5.0/numerator) * consdata->c * MAX(q_var_lb, -q_var_ub) / (consdata->A * 1e5);

      if ( SCIPisFeasGT(scip, p_min, SCIPcomputeVarLbLocal(scip, consdata->p_out_var)) )
      {
         SCIP_CALL( SCIPinferVarLbCons(scip, consdata->p_out_var, p_min, cons, MACH_BOUND, TRUE, &infeasible, &tightened) );
         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         if ( tightened )
         {
            ++(*nchgbds);
            consdata->propvarind |= PIPEODE_POUT_LB;
            SCIPdebugMessage("New lower bound on p_out: <%f> \n",p_min);
            *result = SCIP_REDUCEDDOM;
         }
      }

      if ( SCIPisFeasGT(scip, p_min, SCIPcomputeVarLbLocal(scip, consdata->p_in_var)) )
      {
         SCIP_CALL( SCIPinferVarLbCons(scip, consdata->p_in_var, p_min, cons, MACH_BOUND, TRUE, &infeasible, &tightened) );
         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         if ( tightened )
         {
            ++(*nchgbds);
            consdata->propvarind |= PIPEODE_PIN_LB;
            SCIPdebugMessage("New lower bound on p_in: <%f> \n",p_min);
            *result = SCIP_REDUCEDDOM;
         }
      }
   }

   /* bounds on the flow always possible */
   q_max = MIN(SCIPcomputeVarUbLocal(scip, consdata->p_out_var), SCIPcomputeVarUbLocal(scip, consdata->p_in_var)) * 1e5; /*lint !e666 */
   q_max *= (numerator * consdata->A) / (5.0 * consdata->c);
   if ( SCIPisFeasLT(scip, q_max, q_var_ub) )
   {
      SCIP_CALL( SCIPinferVarUbCons(scip, consdata->q_var, q_max, cons, MACH_BOUND, TRUE, &infeasible, &tightened) );
      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if ( tightened )
      {
         ++(*nchgbds);
         consdata->propvarind |= PIPEODE_Q_UB;
         SCIPdebugMessage("New upper bound on q: <%f> \n", q_max);
         *result = SCIP_REDUCEDDOM;
      }
   }

   if ( SCIPisFeasGT(scip, -q_max, q_var_lb) )
   {
      SCIP_CALL( SCIPinferVarLbCons(scip, consdata->q_var, -q_max, cons, MACH_BOUND, TRUE, &infeasible, &tightened) );
      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if ( tightened )
      {
         ++(*nchgbds);
         consdata->propvarind |= PIPEODE_Q_LB;
         SCIPdebugMessage("New lower bound on q: <%f> \n", -q_max);
         *result = SCIP_REDUCEDDOM;
      }
   }

   return SCIP_OKAY;
}

/** function to propagate pressure bounds from both sides of the pipe to the other side */
static
SCIP_RETCODE boundPropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of the constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_RESULT*          result,             /**< result code */
   int*                  nchgbds             /**< total number of variable bounds tightened by this presolver */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* qvar;
   SCIP_VAR* pin_var;
   SCIP_VAR* pout_var;

   SCIP_Real qvar_ub;
   SCIP_Real qvar_lb;

   SCIP_Real new_p_bound;
   SCIP_Real old_p_bound;

   unsigned int propindex = PIPEODE_NONE;
   int dir = 0;

   SCIP_Bool infeasible;
   SCIP_Bool tightened;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( cons != NULL );
   assert( result != NULL );
   assert( nchgbds != NULL );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   consdata = SCIPconsGetData(cons);
   qvar = consdata->q_var;
   qvar_lb = SCIPcomputeVarLbLocal(scip, qvar);
   qvar_ub = SCIPcomputeVarUbLocal(scip, qvar);

   if ( probdata->noFlowBinvars )
   {
      if ( (consdata->propvarind & PIPEODE_NEG_UB) || (consdata->propvarind & PIPEODE_NEG_LB) || (consdata->propvarind & PIPEODE_POS_UB) || (consdata->propvarind & PIPEODE_POS_LB))
      {
         SCIPerrorMessage("Bound change indicator set to a change of a flow direction variable, although none exists!\n");
         return SCIP_ERROR;
      }
   }
   else
   {
      if ( consdata->propvarind & PIPEODE_NEG_UB )
      {
         consdata->propvarind = consdata->propvarind ^ PIPEODE_NEG_UB;

         if ( SCIPisFeasNegative(scip, qvar_ub) )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisNegative(scip, qvar_lb) )
         {
            SCIP_CALL( SCIPinferVarLbCons(scip, qvar, 0.0, cons, NEG_UB, TRUE, &infeasible, &tightened) );
            assert( ! infeasible );

            if ( tightened )
            {
               ++(*nchgbds);
               SCIPdebugMessage("Only nonnegative flow on pipe <%s> possible.\n", SCIPconsGetName(cons));

               qvar_lb = 0.0;
               consdata->propvarind |= PIPEODE_Q_LB;
            }
         }
      }

      if ( consdata->propvarind & PIPEODE_NEG_LB )
      {
         consdata->propvarind = consdata->propvarind ^ PIPEODE_NEG_LB;

         if ( SCIPisFeasPositive(scip, qvar_lb) )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         if ( SCIPisPositive(scip, qvar_ub) )
         {
            SCIP_CALL( SCIPinferVarUbCons(scip, qvar, 0.0, cons, NEG_LB, TRUE, &infeasible, &tightened) );
            assert( ! infeasible );

            if ( tightened )
            {
               ++(*nchgbds);
               SCIPdebugMessage("Only nonpositive flow on pipe <%s> possible.\n", SCIPconsGetName(cons));
               qvar_ub = 0.0;
               consdata->propvarind |= PIPEODE_Q_UB;
            }
         }
      }

      if ( consdata->propvarind & PIPEODE_POS_UB )
      {
         consdata->propvarind = consdata->propvarind ^ PIPEODE_POS_UB;

         if ( SCIPisFeasPositive(scip, qvar_lb) )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisPositive(scip, qvar_ub) )
         {
            SCIP_CALL( SCIPinferVarUbCons(scip, qvar, 0.0, cons, POS_UB, TRUE, &infeasible, &tightened) );
            assert( ! infeasible );

            if ( tightened )
            {
               ++(*nchgbds);
               SCIPdebugMessage("Only nonpositive flow on pipe <%s> possible.\n", SCIPconsGetName(cons));
               qvar_ub = 0.0;
               consdata->propvarind |= PIPEODE_Q_UB;
            }
         }
      }

      if ( consdata->propvarind & PIPEODE_POS_LB )
      {
         consdata->propvarind = consdata->propvarind ^ PIPEODE_POS_LB;

         if ( SCIPisFeasNegative(scip, qvar_ub) )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisNegative(scip, qvar_lb) )
         {
            SCIP_CALL( SCIPinferVarLbCons(scip, qvar, 0.0, cons, POS_LB, TRUE, &infeasible, &tightened) );
            assert( ! infeasible );

            if ( tightened )
            {
               ++(*nchgbds);
               SCIPdebugMessage("Only nonnegative flow on pipe <%s> possible.\n", SCIPconsGetName(cons));
               qvar_lb = 0.0;
               consdata->propvarind |= PIPEODE_Q_LB;
            }
         }
      }
   }

   /* maybe we can fix the direction of the flow */
   if ( SCIPisFeasNegative(scip, qvar_lb) && SCIPisFeasPositive(scip, qvar_ub)  )
   {
      if ( pressureIsGT(SCIPcomputeVarLbLocal(scip, consdata->p_in_var), SCIPcomputeVarUbLocal(scip, consdata->p_out_var), conshdlrdata->pressure_feastol) )
      {
         /* negative flow is not possible! */
         SCIP_CALL( SCIPinferVarLbCons(scip, qvar, 0.0, cons, PIN_GT_POUT, TRUE, &infeasible, &tightened) );
         assert( ! infeasible );

         if ( tightened )
         {
            ++(*nchgbds);
            SCIPdebugMessage("Only nonnegative flow on pipe <%s> possible.\n", SCIPconsGetName(cons));

            qvar_lb = 0.0;
            /* qvar_lb == 0, i.e. there is no use in propagating bounds with the lower flow bound */
            if ( consdata->propvarind & PIPEODE_PIN_UB )
               consdata->propvarind = consdata->propvarind ^ PIPEODE_PIN_UB;
            if ( consdata->propvarind & PIPEODE_POUT_LB )
               consdata->propvarind = consdata->propvarind ^ PIPEODE_POUT_LB;
            if ( consdata->propvarind & PIPEODE_Q_LB )
               consdata->propvarind = consdata->propvarind ^ PIPEODE_Q_LB;
         }

         if ( ! probdata->noFlowBinvars )
         {
            /* fix binary flow direction variable */
            SCIP_CALL( SCIPinferVarUbCons(scip, consdata->negFlowBinvar, 0.0, cons, PIN_GT_POUT, FALSE, &infeasible, &tightened) );

            if ( infeasible )
            {
               SCIPdebugMessage("NegFlowBinvar of pipe <%s> fixed to 1.0 although only positive flow is possible.\n", SCIPconsGetName(cons));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            if ( tightened )
            {
               ++(*nchgbds);
            }
         }
      }

      if ( pressureIsLT(SCIPcomputeVarUbLocal(scip, consdata->p_in_var), SCIPcomputeVarLbLocal(scip, consdata->p_out_var), conshdlrdata->pressure_feastol) )
      {
         /* positive flow is not possible! */
         SCIP_CALL( SCIPinferVarUbCons(scip, consdata->q_var, 0.0, cons, PIN_LT_POUT, FALSE, &infeasible, &tightened) );
         assert( ! infeasible );

         if ( tightened )
         {
            ++(*nchgbds);
            SCIPdebugMessage("Only nonpositive flow on pipe <%s> possible.\n", SCIPconsGetName(cons));

            qvar_ub = 0.0;

            /* qvar_ub == 0, i.e. there is no use in propagating bounds with the upper flow bound */
            if ( consdata->propvarind & PIPEODE_Q_UB )
               consdata->propvarind = consdata->propvarind ^ PIPEODE_Q_UB;
            if ( consdata->propvarind & PIPEODE_PIN_LB )
               consdata->propvarind = consdata->propvarind ^ PIPEODE_PIN_LB;
            if ( consdata->propvarind & PIPEODE_POUT_UB )
               consdata->propvarind = consdata->propvarind ^ PIPEODE_POUT_UB;
         }

         if ( ! probdata->noFlowBinvars )
         {
            /* fix binary flow direction variable */
            SCIP_CALL( SCIPinferVarUbCons(scip, consdata->posFlowBinvar, 0.0, cons, PIN_LT_POUT, FALSE, &infeasible, &tightened) );

            if ( infeasible )
            {
               SCIPdebugMessage("PosFlowBinvar of pipe <%s> fixed to 1.0 although only negative flow is possible.\n", SCIPconsGetName(cons));
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }

            if ( tightened )
               ++(*nchgbds);
         }
      }
   }

   /* reset indicator for binary variables */
   if ( consdata->propvarind & PIPEODE_NEG_LB )
      consdata->propvarind ^= PIPEODE_NEG_LB;
   if ( consdata->propvarind & PIPEODE_NEG_UB )
      consdata->propvarind ^= PIPEODE_NEG_UB;
   if ( consdata->propvarind & PIPEODE_POS_LB )
      consdata->propvarind ^= PIPEODE_POS_LB;
   if ( consdata->propvarind & PIPEODE_POS_UB )
      consdata->propvarind ^= PIPEODE_POS_UB;

   if ( qvar_ub <= 0.0 )
   {
      dir = -1;
      pin_var = consdata->p_out_var;
      pout_var = consdata->p_in_var;
      qvar_ub = - SCIPcomputeVarLbLocal(scip, qvar);
      qvar_lb = - SCIPcomputeVarUbLocal(scip, qvar);
   }
   else if ( qvar_lb >= 0.0 )
   {
      dir = 1;
      pin_var = consdata->p_in_var;
      pout_var = consdata->p_out_var;
   }
   else
   {
      qvar_lb *= (-1);
      pin_var = consdata->p_in_var;
      pout_var = consdata->p_out_var;
   }

   /* choose propindex such that pin/pout coincides with the local definition of in and out */
   if ( dir == -1 )
   {
      if ( consdata->propvarind & PIPEODE_Q_UB )
      {
         propindex |= PIPEODE_Q_LB;
         consdata->propvarind ^= PIPEODE_Q_UB;
      }
      if ( consdata->propvarind & PIPEODE_Q_LB )
      {
         propindex |= PIPEODE_Q_UB;
         consdata->propvarind ^= PIPEODE_Q_LB;
      }
      if ( consdata->propvarind & PIPEODE_POUT_UB )
      {
         propindex |= PIPEODE_PIN_UB;
         consdata->propvarind ^= PIPEODE_POUT_UB;
      }
      if ( consdata->propvarind & PIPEODE_POUT_LB )
      {
         propindex |= PIPEODE_PIN_LB;
         consdata->propvarind ^= PIPEODE_POUT_LB;
      }
      if ( consdata->propvarind & PIPEODE_PIN_UB )
      {
         propindex |= PIPEODE_POUT_UB;
         consdata->propvarind ^= PIPEODE_PIN_UB;
      }
      if ( consdata->propvarind & PIPEODE_PIN_LB )
      {
         propindex |= PIPEODE_POUT_LB;
         consdata->propvarind ^= PIPEODE_PIN_LB;
      }
      assert( consdata->propvarind == PIPEODE_NONE );
   }
   else
   {
      propindex = consdata->propvarind;
      consdata->propvarind = PIPEODE_NONE;
   }

   /*
    *  Upper Bound of Flow has changed
    */
   if ( propindex & PIPEODE_Q_UB )
   {
      /* check for new upper bound on pin */
      new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata,consdata, SCIPcomputeVarUbLocal(scip, pout_var), qvar_ub, -1);
      if ( new_p_bound > -0.5 )
      {
         old_p_bound = SCIPcomputeVarUbLocal(scip, pin_var);
         if ( SCIPisFeasLT(scip, new_p_bound, SCIPcomputeVarLbLocal(scip, pin_var)) )
         {
            /* there is no feasible solution! */
            SCIPdebugMessage("boundPropagation cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisFeasLT(scip, new_p_bound, old_p_bound) )
         {
            if ( dir == -1 )
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, pin_var, new_p_bound, cons, FPROP_UPPER, FALSE, &infeasible, &tightened) );
            }
            else
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, pin_var, new_p_bound, cons, BPROP_UPPER, FALSE, &infeasible, &tightened) );
            }

            assert( ! infeasible );
            if ( tightened )
            {
               ++(*nchgbds);
               propindex |= PIPEODE_PIN_UB;
               SCIPdebugMessage("New upper bound on %s: <%f> \n", SCIPvarGetName(pin_var), new_p_bound);
            }
         }
      }

      /* check for new lower bound on pout */
      new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata, consdata, SCIPcomputeVarLbLocal(scip, pin_var), qvar_ub, 1);
      if ( new_p_bound > -0.5 )
      {
         old_p_bound = SCIPcomputeVarLbLocal(scip, pout_var);
         if ( SCIPisFeasGT(scip, new_p_bound, SCIPcomputeVarUbLocal(scip, pout_var)) )
         {
            /* there is no feasible solution! */
            SCIPdebugMessage("boundPropagation cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisFeasGT(scip, new_p_bound, old_p_bound) )
         {
            if ( dir == -1 )
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, pout_var, new_p_bound, cons, BPROP_LOWER, FALSE, &infeasible, &tightened) );
            }
            else
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, pout_var, new_p_bound, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );
            }

            assert( ! infeasible );
            if ( tightened )
            {
               ++(*nchgbds);
               propindex |= PIPEODE_POUT_LB;
               SCIPdebugMessage("New lower bound on %s: <%f> \n",SCIPvarGetName(pout_var), new_p_bound);
            }
         }
      }

      /* we do not need to propagate the upper bound on pout again*/
      if ( propindex & PIPEODE_POUT_UB )
         propindex ^= PIPEODE_POUT_UB;

      /* and we do not need to propagate the lower bound on pin again */
      if ( propindex & PIPEODE_PIN_LB )
         propindex ^= PIPEODE_PIN_LB;

      propindex ^= PIPEODE_Q_UB;
   }

   /*
    *  Lower Bound of Flow has changed
    */
   if ( propindex & PIPEODE_Q_LB )
   {
      /* check for new upper bound on pout */
      if ( dir == 0 )
         new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata, consdata, SCIPcomputeVarUbLocal(scip, pin_var), qvar_lb, -1);
      else
         new_p_bound = computeModifiedEulerScheme(scip, conshdlrdata, consdata, SCIPcomputeVarUbLocal(scip, pin_var), qvar_lb, 1);

      if ( new_p_bound > -0.5 )
      {
         old_p_bound = SCIPcomputeVarUbLocal(scip, pout_var);
         if ( SCIPisFeasLT(scip, new_p_bound, SCIPcomputeVarLbLocal(scip, pout_var)) )
         {
            /* there is no feasible solution! */
            SCIPdebugMessage("boundPropagation cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisFeasLT(scip, new_p_bound, old_p_bound) )
         {
            if ( dir == -1 )
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, pout_var, new_p_bound, cons, BPROP_UPPER, FALSE, &infeasible, &tightened) );
            }
            else
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, pout_var, new_p_bound, cons, FPROP_UPPER, FALSE, &infeasible, &tightened) );
            }

            assert( ! infeasible );
            if ( tightened )
            {
               ++(*nchgbds);
               propindex |= PIPEODE_POUT_UB;
               SCIPdebugMessage("New upper bound on %s: <%f> \n",SCIPvarGetName(pout_var), new_p_bound);
            }
         }
      }

      /* check for new lower bound on pin */
      if ( dir == 0 )
         new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata, consdata, SCIPcomputeVarLbLocal(scip, pout_var), qvar_lb, 1);
      else
         new_p_bound = computeModifiedEulerScheme(scip, conshdlrdata, consdata, SCIPcomputeVarLbLocal(scip, pout_var), qvar_lb, -1);

      if ( new_p_bound > -0.5 )
      {
         old_p_bound = SCIPcomputeVarLbLocal(scip, pin_var);
         if ( SCIPisFeasGT(scip, new_p_bound, SCIPcomputeVarUbLocal(scip, pin_var)) )
         {
            /* there is no feasible solution! */
            SCIPdebugMessage("boundPropagation cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisFeasGT(scip, new_p_bound, old_p_bound) )
         {
            if ( dir == -1 )
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, pin_var, new_p_bound, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );
            }
            else
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, pin_var, new_p_bound, cons, BPROP_LOWER, FALSE, &infeasible, &tightened) );
            }

            assert( ! infeasible );
            if ( tightened )
            {
               ++(*nchgbds);
               propindex |= PIPEODE_PIN_LB;
               SCIPdebugMessage("New lower bound on %s: <%f> \n",SCIPvarGetName(pin_var), new_p_bound);
            }
         }
      }

      /* we already propagated the upper bound on pin */
      if ( propindex & PIPEODE_PIN_UB )
         propindex ^= PIPEODE_PIN_UB;

      /* and also the lower bound on pout */
      if ( propindex & PIPEODE_POUT_LB )
         propindex ^= PIPEODE_POUT_LB;

      propindex ^= PIPEODE_Q_LB;
   }

   /*
    *  Upper Bound of POUT has changed
    */
   if ( propindex & PIPEODE_POUT_UB )
   {
      /* check for new upper bound on pin */
      new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata, consdata, SCIPcomputeVarUbLocal(scip, pout_var), qvar_ub, -1);

      if ( new_p_bound > -0.5 )
      {
         old_p_bound = SCIPcomputeVarUbLocal(scip, pin_var);
         if ( SCIPisFeasLT(scip, new_p_bound, SCIPcomputeVarLbLocal(scip, pin_var)) )
         {
            /* there is no feasible solution! */
            SCIPdebugMessage("boundPropagation cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisFeasLT(scip, new_p_bound, old_p_bound) )
         {
            if ( dir == -1 )
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, pin_var, new_p_bound, cons, FPROP_UPPER, FALSE, &infeasible, &tightened) );
            }
            else
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, pin_var, new_p_bound, cons, BPROP_UPPER, FALSE, &infeasible, &tightened) );
            }

            assert( ! infeasible );
            if ( tightened )
            {
               ++(*nchgbds);
               SCIPdebugMessage("New upper bound on %s: <%f> \n",SCIPvarGetName(pin_var), new_p_bound);
            }
         }
      }

      propindex ^= PIPEODE_POUT_UB;
   }

   /*
    *  Lower Bound of POUT has changed
    */
   if ( propindex & PIPEODE_POUT_LB )
   {
      /* check for new lower bound on pin */
      if ( dir == 0 )
         new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata, consdata, SCIPcomputeVarLbLocal(scip, pout_var), qvar_lb, 1);
      else
         new_p_bound = computeModifiedEulerScheme(scip, conshdlrdata, consdata, SCIPcomputeVarLbLocal(scip, pout_var), qvar_lb, -1);

      if ( new_p_bound > -0.5 )
      {
         old_p_bound = SCIPcomputeVarLbLocal(scip, pin_var);
         if ( SCIPisFeasGT(scip, new_p_bound, SCIPcomputeVarUbLocal(scip, pin_var)) )
         {
            /* there is no feasible solution! */
            SCIPdebugMessage("boundPropagation cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisFeasGT(scip, new_p_bound, old_p_bound) )
         {
            if ( dir == -1 )
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, pin_var, new_p_bound, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );
            }
            else
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, pin_var, new_p_bound, cons, BPROP_LOWER, FALSE, &infeasible, &tightened) );
            }

            assert( ! infeasible );
            if ( tightened )
            {
               ++(*nchgbds);
               SCIPdebugMessage("New lower bound on %s: <%f> \n",SCIPvarGetName(pin_var), new_p_bound);
            }
         }
      }

      propindex ^= PIPEODE_POUT_LB;
   }

   /*
    *  Upper Bound of PIN has changed
    */
   if ( propindex & PIPEODE_PIN_UB )
   {
      /* check for new upper bound on pout */
      if ( dir == 0 )
         new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata, consdata, SCIPcomputeVarUbLocal(scip, pin_var), qvar_lb, -1);
      else
         new_p_bound = computeModifiedEulerScheme(scip, conshdlrdata, consdata, SCIPcomputeVarUbLocal(scip, pin_var), qvar_lb, 1);

      if ( new_p_bound > -0.5 )
      {
         old_p_bound = SCIPcomputeVarUbLocal(scip, pout_var);

         if ( SCIPisFeasLT(scip, new_p_bound, SCIPcomputeVarLbLocal(scip, pout_var)) )
         {
            /* there is no feasible solution! */
            SCIPdebugMessage("boundPropagation cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisFeasLT(scip, new_p_bound, old_p_bound) )
         {
            if ( dir == -1 )
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, pout_var, new_p_bound, cons, BPROP_UPPER, FALSE, &infeasible, &tightened) );
            }
            else
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, pout_var, new_p_bound, cons, FPROP_UPPER, FALSE, &infeasible, &tightened) );
            }
            assert( ! infeasible );

            if ( tightened )
            {
               ++(*nchgbds);
               SCIPdebugMessage("New upper bound on %s: <%f> \n",SCIPvarGetName(pout_var), new_p_bound);
            }
         }
      }

      propindex ^= PIPEODE_PIN_UB;
   }

   /*
    *  Lower Bound of PIN has changed
    */
   if ( propindex & PIPEODE_PIN_LB )
   {
      /* check for new lower bound on pout */
      new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata, consdata, SCIPcomputeVarLbLocal(scip, pin_var), qvar_ub, 1);

      if ( new_p_bound > -0.5 )
      {
         old_p_bound = SCIPcomputeVarLbLocal(scip, pout_var);
         if ( SCIPisFeasGT(scip, new_p_bound, SCIPcomputeVarUbLocal(scip, pout_var)) )
         {
            /* there is no feasible solution! */
            SCIPdebugMessage("boundPropagation cutoff.\n");
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }
         else if ( SCIPisFeasGT(scip, new_p_bound, old_p_bound) )
         {
            if ( dir == -1 )
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, pout_var, new_p_bound, cons, BPROP_LOWER, FALSE, &infeasible, &tightened) );
            }
            else
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, pout_var, new_p_bound, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );
            }

            assert( ! infeasible );
            if ( tightened )
            {
               ++(*nchgbds);
               SCIPdebugMessage("New lower bound on %s: <%f> \n",SCIPvarGetName(pout_var), new_p_bound);
            }
         }
      }

      propindex ^= PIPEODE_PIN_LB;
   }
   assert( propindex == 0 );

   return SCIP_OKAY;
}

/** function that determines the shape of the feasible region (p, q) with respect to the variable
 *  bounds and the "mach constraint" 4 * A * p > 5 * c * |q| */
static
int determineShape(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             q_var_ub,           /**< Variable of the upper flow bound */
   SCIP_Real             q_var_lb,           /**< Variable of the lower flow bound */
   SCIP_Real             p_var_ub,           /**< Variable of the upper pressure bound */
   SCIP_Real             p_var_lb,           /**< Variable of the lower pressure bound */
   SCIP_Real             A,                  /**< Cross-sectional area of the pipe */
   SCIP_Real             D,                  /**< Diameter of the pipe */
   SCIP_Real             cc,                 /**< Speed of Sound */
   int                   mach_num            /**< numerator used for bound on mach number: v/c = machbound/5 */
   )
{
   SCIP_Real mach1;
   SCIP_Real mach2;
   SCIP_Real mach3;
   SCIP_Real machbound;
   int shape = 0;
   int shape_plus = 0;
   int shape_minus = 0;

   assert( scip != NULL );
   assert( SCIPisGT(scip, q_var_ub, q_var_lb) );
   assert( SCIPisGT(scip, p_var_ub, p_var_lb) );

   machbound = ((SCIP_Real)mach_num)/5.0;

   if ( SCIPisFeasLE(scip, q_var_ub, 0.0) )
   {
      /* if flow is negative, change sign and bounds */
      q_var_ub += q_var_lb;
      q_var_lb -= q_var_ub;
      q_var_ub += q_var_lb;
      q_var_ub *= -1;
   }

   if ( SCIPisFeasGE(scip, q_var_lb, 0.0) )
   {
      /* assert that the problem is feasible */
      assert( ( cc * q_var_lb ) / ( A * 1e5 * p_var_ub ) <= machbound );
      /* p lower, q lower */
      mach1 = ( cc * q_var_lb ) / ( A * p_var_lb * 1e5);
      /* p upper, q upper */
      mach2 = ( cc * q_var_ub ) / ( A * p_var_ub * 1e5);
      /* p lower, q upper */
      mach3 = ( cc * q_var_ub ) / ( A * p_var_lb * 1e5);

      if ( SCIPisFeasLE(scip, mach3, machbound) )
      {
         shape = 1;
      }
      else
      {
         if ( SCIPisFeasGE(scip, mach1, machbound) )
         {
            /* 3 oder 4 */
            if ( SCIPisFeasGE(scip, mach2, machbound) )
            {
               shape = 3;
            }
            else
            {
               shape = 4;
            }
         }
         else
         {
            /* 2 oder 5 */
            if ( SCIPisFeasGE(scip, mach2, machbound) )
            {
               shape = 2;
            }
            else
            {
               shape = 5;
            }
         }
      }
   }

   if ( q_var_lb < 0.0 && q_var_ub > 0.0 )
   {
      shape_plus = determineShape(scip, q_var_ub, 0.0, p_var_ub, p_var_lb, A, D, cc, mach_num);
      shape_minus = determineShape(scip, -q_var_lb, 0.0, p_var_ub, p_var_lb, A, D, cc, mach_num);

      assert( shape_plus != 3 && shape_plus != 4 );
      assert( shape_minus != 3 && shape_minus != 4 );

      switch ( shape_plus )
      {
      case 1:
         switch ( shape_minus )
         {
         case 1: shape = 11; break;
         case 2: shape = 16; break;
         case 5: shape = 13; break;
         default: abort();
         }
         break;
      case 2:
         switch ( shape_minus )
         {
         case 1: shape = 15; break;
         case 2: shape = 17; break;
         case 5: shape = 18; break;
         default: abort();
         }
         break;
      case 5:
         switch ( shape_minus )
         {
         case 1: shape = 12; break;
         case 2: shape = 19; break;
         case 5: shape = 14; break;
         default: abort();
         }
         break;
      default: abort();
      }
   }

   assert( shape > 0 );

   return shape;
}


/** function that returns a new bound for the flow by the flow pointer and
 *  a boolean which indicates if we a valid new bound was found.
 *  Finds a new bound by bisection.
 */
static
SCIP_Bool bisectedFlowBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data of cons_pipe_ode */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            flow,               /**< pointer to new flow bound */
   int                   dir                 /**< -1/+1, which bound to improve */
   )
{  /*lint --e{438}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR* q;
   SCIP_VAR* pout;
   SCIP_Real pout_lb;
   SCIP_Real pin_ub;
   SCIP_Real out_pressure;
   SCIP_Real lower_flow = 0.0;
   SCIP_Real upper_flow = 0.0;
   SCIP_Real lower_in_pressure;
   SCIP_Real upper_in_pressure;
   SCIP_Real mid_flow = 0.0;
   SCIP_Real mid_in_pressure = 0.0;
   SCIP_Real scalar = 0.0;
   SCIP_Real numerator;
   SCIP_Real delta = 1.0;                    /* tolerance for computing new upper flow bound */
   SCIP_Bool foundNewBound = FALSE;
   int i = 0;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( cons != NULL );
   assert( flow != NULL );
   assert( dir == 1 || dir == -1);

   consdata = SCIPconsGetData(cons);

   q = consdata->q_var;

   numerator = (SCIP_Real) conshdlrdata->machbound;

   /* do so as if q was only positive: */
   if ( dir == 1 )
   {
      lower_flow = SCIPcomputeVarLbLocal(scip, q);
      lower_flow = MAX(0.0, lower_flow);
      upper_flow = SCIPcomputeVarUbLocal(scip, q);

      pout = consdata->p_out_var;
      pin_ub = SCIPcomputeVarUbLocal(scip, consdata->p_in_var);
   }
   else
   {
      lower_flow = - SCIPcomputeVarUbLocal(scip, q);
      lower_flow = MAX(0.0, lower_flow);
      upper_flow = - SCIPcomputeVarLbLocal(scip, q);

      pout = consdata->p_in_var;
      pin_ub = SCIPcomputeVarUbLocal(scip, consdata->p_out_var);
   }

   pout_lb = SCIPcomputeVarLbLocal(scip, pout);

   /* compute propagation of lower output pressure */
   /* pout_lb should already satisfy the mach bound, since the function machBound should be called before flowTighteneing */
   out_pressure = pout_lb;
   lower_in_pressure = computeModifiedEulerScheme(scip, conshdlrdata, consdata, out_pressure, lower_flow, -1);

   /* compute propagation with upper flow bound */
   /* due to the function machBound q_ub should be so small such that a "mach"-feasible pressure exists! */
   out_pressure = MAX( pout_lb , (5.0 * consdata->c * upper_flow ) / (numerator * consdata->A * 1e5) );
   assert( SCIPisFeasLE(scip, out_pressure, SCIPcomputeVarUbLocal(scip, pout)) );
   upper_in_pressure = computeModifiedEulerScheme(scip, conshdlrdata, consdata, out_pressure, upper_flow, -1);

   if ( SCIPisFeasLE(scip, upper_in_pressure, pin_ub + delta) )
   {
      /* flow bound cannot be reduced */
      return FALSE;
   }

   if ( SCIPisFeasLT(scip, pin_ub, lower_in_pressure) )
   {
      /* there is no feasible soltion */
      consdata->propvarind |= PIPEODE_POUT_LB; /* make sure boundPropagation detects infeasibility */
      return FALSE;
   }

   while ( i < 5 )
   {
      scalar = (pin_ub + delta/2 -  lower_in_pressure) / ( upper_in_pressure - lower_in_pressure);
      mid_flow = scalar * upper_flow + (1 - scalar) * lower_flow;
      out_pressure = MAX( pout_lb , (5.0 * consdata->c * mid_flow ) / (numerator * consdata->A * 1e5) );
      mid_in_pressure = computeModifiedEulerScheme(scip, conshdlrdata, consdata, out_pressure, mid_flow, -1);

      if ( SCIPisFeasGT(scip, mid_in_pressure, pin_ub + delta) )
      {
         upper_flow = mid_flow;
         upper_in_pressure = mid_in_pressure;
         foundNewBound = TRUE;
      }
      else if ( SCIPisFeasLE(scip, mid_in_pressure, pin_ub) )
      {
         lower_flow = mid_flow;
         lower_in_pressure = mid_in_pressure;
      }
      else
      {
         (*flow) = mid_flow;
         foundNewBound = TRUE;
         upper_in_pressure = mid_in_pressure;
         break;
      }

      ++i;
   }

   if ( i == 5 )
   {
      (*flow) = upper_flow;
      out_pressure = MAX( pout_lb , (5.0 * consdata->c * upper_flow ) / (numerator * consdata->A * 1e5) );
   }

   /* correct sign of the new flow bound: */
   (*flow) *= dir;

   return foundNewBound;
}

/** function to decrease flow intervalls during presolving!
 *  Function machBound should always be called before flowTightening!
 */
static
SCIP_RETCODE flowTightening(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS*            cons,               /**< constraint */
   int*                  nchgbds             /**< total number of variable bounds tightened by this presolver */
   )
{
   SCIP_CONSDATA* consdata;

   SCIP_VAR* qvar;
   SCIP_Real qvar_ub;
   SCIP_Real qvar_lb;
   int dir = 0;
   unsigned int propindex;

   SCIP_Real flowBound = 0.0;

   SCIP_Bool infeasible = FALSE;
   SCIP_Bool tightened = FALSE;
   SCIP_Bool foundNewBound = FALSE;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( cons != NULL );
   assert( nchgbds != NULL );

   consdata = SCIPconsGetData(cons);
   qvar = consdata->q_var;
   qvar_lb = SCIPcomputeVarLbLocal(scip, qvar);
   qvar_ub = SCIPcomputeVarUbLocal(scip, qvar);
   propindex = consdata->propvarind;

   if( SCIPisLE(scip, qvar_ub - qvar_lb, 10.0) )
   {
      return SCIP_OKAY;
   }

   if ( qvar_ub <= 0.0 )
   {
      dir = -1;
   }
   else if ( qvar_lb >= 0.0 )
   {
      dir = 1;
   }

   if ( dir == 0 || dir == 1 )
   {
      if ( (propindex & PIPEODE_PIN_UB) || (propindex & PIPEODE_POUT_LB) )
      {
         /* apply optimization of flow bound to the positive flow intervall */
         foundNewBound = bisectedFlowBound(scip, conshdlrdata, cons, &flowBound, 1);

         if ( foundNewBound )
         {
            SCIP_CALL( SCIPinferVarUbCons(scip, qvar, flowBound, cons, BPROP_LOWER, FALSE, &infeasible, &tightened) );
            assert( ! infeasible );

            if ( tightened )
            {
               ++(*nchgbds);
               consdata->propvarind |= PIPEODE_Q_UB;
            }
         }
      }
   }

   infeasible = FALSE;
   tightened = FALSE;

   if ( dir == 0 || dir == -1 )
   {
      if ( (propindex & PIPEODE_PIN_LB) || (propindex & PIPEODE_POUT_UB) )
      {
         /* apply optimization of flow bound to the positive flow intervall */
         foundNewBound = bisectedFlowBound(scip, conshdlrdata, cons, &flowBound, -1);

         if ( foundNewBound )
         {
            SCIP_CALL( SCIPinferVarLbCons(scip, qvar, flowBound, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );
            assert( ! infeasible );

            if ( tightened )
            {
               ++(*nchgbds);
               consdata->propvarind |= PIPEODE_Q_LB;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** function to generate the linear constraints for the concave overestimator of the pressure-input function.
 *  adds the overestimator to the lp:
 *                 - if the solution is separated
 *                 - if no solution is given, i.e. pointer are NULL */
static
SCIP_RETCODE generateCaveNothingFixed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler data */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            p_in_sol,           /**< pointer to the (pseudo)solution of input pressure */
   SCIP_Real*            p_out_sol,          /**< pointer to the (pseudo)solution of output pressure */
   SCIP_Real*            q_sol               /**< pointer to the (pseudo)solution of the flow */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_ROW*      cave[3];
   SCIP_Real      q_var_ub;
   SCIP_Real      q_var_lb;
   SCIP_Real      p_var_ub;
   SCIP_Real      p_var_lb;
   SCIP_Real      pout[5];
   SCIP_Real      pin[5];
   SCIP_Real      flow[5];
   SCIP_Real      rhs[3];
   SCIP_Real      n_pout[3];
   SCIP_Real      n_q[3];
   int            shape;
   int            i;
   int            j;
   int            dir;
   SCIP_Bool      addcut     = TRUE;
   SCIP_Bool      infeasible = TRUE;
   char           s[SCIP_MAXSTRLEN];
   SCIP_Real      numerator;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );

   if ( p_in_sol != NULL || p_out_sol != NULL || q_sol != NULL )
   {
      assert( p_in_sol != NULL );
      assert( p_out_sol != NULL );
      assert( q_sol != NULL );
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   numerator = (SCIP_Real) conshdlrdata->machbound;

   if ( ! flowIsNegative(SCIPcomputeVarLbLocal(scip, consdata->q_var), conshdlrdata->flow_feastol) )
   {
      q_var_lb = MAX(0.0, SCIPcomputeVarLbLocal(scip, consdata->q_var)); /*lint !e666 */
      q_var_ub = SCIPcomputeVarUbLocal(scip, consdata->q_var);
      p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_out_var);
      p_var_lb = SCIPcomputeVarLbLocal(scip, consdata->p_out_var);
      dir = 1;
      /* assert( *p_in_sol == SCIPgetVarSol(scip, consdata->p_in_var) ); */
      /* assert( *p_out_sol == SCIPgetVarSol(scip, consdata->p_out_var) ); */
   }
   else if ( ! flowIsPositive( SCIPcomputeVarUbLocal(scip, consdata->q_var), conshdlrdata->flow_feastol) )
   {
      q_var_ub = - SCIPcomputeVarLbLocal(scip, consdata->q_var);
      q_var_lb = MAX(0.0, - SCIPcomputeVarUbLocal(scip, consdata->q_var)); /*lint !e666 */
      p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_in_var);
      p_var_lb = SCIPcomputeVarLbLocal(scip, consdata->p_in_var);
      dir = -1;
      /* assert( *p_out_sol == SCIPgetVarSol(scip, consdata->p_in_var) ); */
      /* assert( *p_in_sol == SCIPgetVarSol(scip, consdata->p_out_var) ); */
   }
   else
   {
      SCIPdebugMessage("generateCaveNothingFixed works only, when the direction of the flow is fixed.\n");
      return SCIP_INVALIDDATA;
   }

   assert( SCIPisFeasGE(scip, numerator * consdata->A * p_var_ub * 1e5, 5.0 * consdata->c * q_var_lb) );

   shape = determineShape(scip, q_var_ub, q_var_lb, p_var_ub, p_var_lb, consdata->A, consdata->D, consdata->c, conshdlrdata->machbound);

   assert( shape > 0 );
   assert( shape < 6 );

   pout[1] = p_var_ub;
   flow[1] = q_var_lb;

   /* determine the vertices */
   if ( shape == 3 )
   {
      i = 3;
      pout[0] = MAX(5.0 * consdata->c * q_var_lb / (numerator * consdata->A * 1e5), p_var_lb) ;
      flow[0] = q_var_lb;
      pout[2] = p_var_ub;
      flow[2] = MIN(numerator * (consdata->A) * 1e5 * p_var_ub / (5.0 * consdata->c), q_var_ub);
   }
   else if ( shape == 5 )
   {
      i = 5;
      pout[0] = p_var_lb;
      flow[0] = q_var_lb;
      pout[2] = p_var_ub;
      flow[2] = q_var_ub;
      pout[3] = 5.0 * consdata->c * q_var_ub / (numerator * consdata->A * 1e5);
      flow[3] = q_var_ub;
      pout[4] = p_var_lb;
      flow[4] = numerator * consdata->A * 1e5 * p_var_lb / (5.0 * consdata->c);
   }
   else
   {
      i = 4;
      if ( shape == 1)
      {
         pout[0] = p_var_lb;
         flow[0] = q_var_lb;
         pout[2] = p_var_ub;
         flow[2] = q_var_ub;
         pout[3] = p_var_lb;
         flow[3] = q_var_ub;
      }
      else if ( shape == 2 )
      {
         pout[0] = p_var_lb;
         flow[0] = q_var_lb;
         pout[2] = p_var_ub;
         flow[2] = MIN(numerator * (consdata->A) * 1e5 * p_var_ub / ( 5.0 * (consdata->c)), q_var_ub);
         pout[3] = p_var_lb;
         flow[3] = numerator * (consdata->A) * 1e5 * p_var_lb / ( 5.0 * (consdata->c));
      }
      else if ( shape == 4 )
      {
         pout[0] = MAX(5.0 * consdata->c * q_var_lb / (numerator * consdata->A * 1e5), p_var_lb);
         flow[0] = q_var_lb;
         pout[2] = p_var_ub;
         flow[2] = q_var_ub;
         pout[3] = 5.0 * consdata->c * q_var_ub / (numerator * consdata->A * 1e5);
         flow[3] = q_var_ub;
      }
   }

   /* compute the input pressure with trapezoidal scheme */
   assert( i >= 3 );
   while ( i > 0 )
   {
      --i;
      pin[i] = computeTrapezoidalScheme(scip, conshdlrdata, consdata, pout[i], flow[i], -1);
   }

   /* determine the coefficients and rhs for the triangle 0-1-2 */
   n_pout[0] = - (pin[1] - pin[0]) / (pout[1] - pout[0]); /*lint !e771*/
   n_q[0]    = - (pin[2] - pin[1]) / (flow[2] - flow[1]);
   rhs[0]    = pin[1] + pout[1] * n_pout[0] + flow[1] * n_q[0];

   /*
    *  Case 1
    */
   if ( shape == 1 )
   {
      /* overestimator consists of i = 2 constraints */
      i = 2;
      if ( SCIPisFeasLE(scip, pin[3] + n_pout[0] * pout[3] + n_q[0] * flow[3], rhs[0]) )
      {
         /* determine the coefficients and rhs for the triangle 0-2-3 */
         n_pout[1] = - (pin[2] - pin[3]) / (pout[2] - pout[3]);
         n_q[1]    = - (pin[3] - pin[0]) / (flow[3] - flow[0]);
         rhs[1]    = pin[3] + pout[3] * n_pout[1] + flow[3] * n_q[1];
         if ( SCIPisFeasLE(scip, pin[1] + n_pout[1] * pout[1] + n_q[1] * flow[1], rhs[1]) )
         {
            infeasible = FALSE;
         }
      }

      if ( infeasible )
      {
         /* determine the coefficients and rhs for the triangle 0-1-3 */
         n_pout[0] = - (pin[1] - pin[0]) / (pout[1] - pout[0]);
         n_q[0]    = - (pin[3] - pin[0]) / (flow[3] - flow[0]);
         rhs[0]    = pin[0] + pout[0] * n_pout[0] + flow[0] * n_q[0];

         if ( SCIPisFeasLE(scip, pin[2] + n_pout[0] * pout[2] + n_q[0] * flow[2], rhs[0]) )
         {
            /* determine the coefficients and rhs for the triangle 1-2-3 */
            n_pout[1] = - (pin[2] - pin[3]) / (pout[2] - pout[3]);
            n_q[1]    = - (pin[2] - pin[1]) / (flow[2] - flow[1]);
            rhs[1]    = pin[2] + pout[2] * n_pout[1] + flow[2] * n_q[1];
            if ( SCIPisFeasLE(scip, pin[0] + n_pout[1] * pout[0] + n_q[1] * flow[0], rhs[1]) )
            {
               infeasible = FALSE;
            }
         }
      }

      if ( infeasible )
      {
         i = 0;
         SCIPerrorMessage("In generateCaveNothingFixed: no overerstimator found for Case 1.\n");
      }
   }
   /*
    *  Case 2
    */
   else if ( shape == 2 )
   {
      /* overestimator consists of i = 2 constraints */
      i = 2;
      if ( SCIPisFeasLE(scip, pin[3] + n_pout[0] * pout[3] + n_q[0] * flow[3], rhs[0]) )
      {
         /* determine the coefficients and rhs for the triangle 0-2-3 */
         n_pout[1] = (flow[2] - flow[0]) * (pin[3] - pin[0]) / ((pout[2] -  pout[0]) * (flow[3] - flow[0])) - (pin[2] - pin[0]) / (pout[2] - pout[0]);
         n_q[1]    = - (pin[3] - pin[0]) / (flow[3] - flow[0]);
         rhs[1]    = pin[0] + pout[0] * n_pout[1] + flow[0] * n_q[1];
         if ( SCIPisFeasLE(scip, pin[1] + n_pout[1] * pout[1] + n_q[1] * flow[1], rhs[1]) )
         {
            infeasible = FALSE;
         }
      }

      if ( infeasible )
      {
         /* determine the coefficients and rhs for the triangle 0-1-3 */
         n_pout[0] = - (pin[1] - pin[0]) / (pout[1] - pout[0]);
         n_q[0]    = - (pin[3] - pin[0]) / (flow[3] - flow[0]);
         rhs[0]    = pin[0] + pout[0] * n_pout[0] + flow[0] * n_q[0];

         if ( SCIPisFeasLE(scip, pin[2] + n_pout[0] * pout[2] + n_q[0] * flow[2], rhs[0]) )
         {
            /* determine the coefficients and rhs for the triangle 1-2-3 */
            n_pout[1] = - (flow[3] - flow[1]) * (pin[2] - pin[1]) / ((pout[1] -  pout[3]) * (flow[2] - flow[1])) + (pin[3] - pin[1]) / (pout[1] - pout[3]);
            n_q[1]    = (pin[2] - pin[1]) / (flow[2] - flow[1]);
            rhs[1]    = pin[1] + pout[1] * n_pout[1] + flow[1] * n_q[1];
            if ( SCIPisFeasLE(scip, pin[0] + n_pout[1] * pout[0] + n_q[1] * flow[0], rhs[1]) )
            {
               infeasible = FALSE;
            }
         }
      }

      if ( infeasible )
      {
         i = 0;
         SCIPerrorMessage("In generateCaveNothingFixed: no overerstimator found for Case 2.\n");
      }
   }
   /*
    *  Case 3
    */
   else if ( shape == 3 )
   {
      i = 1;
      infeasible = FALSE;
   }
   /*
    *  Case 4
    */
   else if ( shape == 4 )
   {
      /* overestimator consists of i = 2 constraints */
      i = 2;
      if ( SCIPisFeasLE(scip, pin[3] + n_pout[0] * pout[3] + n_q[0] * flow[3], rhs[0]) )
      {
         /* determine the coefficients and rhs for the triangle 0-2-3 */
         n_pout[1] = - (pin[2] - pin[3]) / (pout[2] - pout[3]);
         n_q[1] = (pout[2] - pout[0]) * (pin[2] - pin[3]) / ((flow[2] - flow[0]) * (pout[2] - pout[3])) - (pin[2] - pin[0]) / (flow[2] - flow[0]);
         rhs[1] = pin[2] + n_pout[1] * pout[2] + n_q[1] * flow[2];
         if ( SCIPisFeasLE(scip, pin[1] + n_pout[1] * pout[1] + n_q[1] * flow[1], rhs[1]) )
         {
            infeasible = FALSE;
         }
      }

      if ( infeasible )
      {
         /* determine the coefficients and rhs for the triangle 0-1-3 */
         n_pout[0] = - (pin[1] - pin[0]) / (pout[1] - pout[0]);
         n_q[0] = (pout[3] - pout[0]) * (pin[1] - pin[0]) / ((flow[3] - flow[0]) * (pout[1] - pout[0])) - (pin[3] - pin[0]) / (flow[3] - flow[0]);
         rhs[0] = pin[0] + n_pout[0] * pout[0] + n_q[0] * flow[0];
         if ( SCIPisFeasLE(scip, pin[2] + n_pout[0] * pout[2] + n_q[0] * flow[2], rhs[0]) )
         {
            /* determine the coefficients and rhs for the triangle 1-2-3 */
            n_pout[1] = - (pin[2] - pin[3]) / (pout[2] - pout[3]);
            n_q[1]    = - (pin[2] - pin[1]) / (flow[2] - flow[1]);
            rhs[1]    = pin[2] + pout[2] * n_pout[1] + flow[2] * n_q[1];
            if ( SCIPisFeasLE(scip, pin[0] + n_pout[1] * pout[0] + n_q[1] * flow[0], rhs[1]) )
            {
               infeasible = FALSE;
            }
         }
      }

      if ( infeasible )
      {
         i = 0;
         SCIPerrorMessage("In generateCaveNothingFixed: no overerstimator found for Case 4.\n");
      }
   }
   /*
    *  Case 5
    */
   else if ( shape == 5 )
   {
      /* overestimator consists of i = 3 constraints */
      i = 3;
      if ( SCIPisFeasLE(scip, pin[3] + n_pout[0] * pout[3] + n_q[0] * flow[3], rhs[0]) && SCIPisFeasLE(scip, pin[4] + n_pout[0] * pout[4] + n_q[0] * flow[4], rhs[0]) )
      {
         /* determine the coefficients and rhs for the triangle 0-2-3 */
         n_pout[1] = - (pin[2] - pin[3]) / (pout[2] - pout[3]);
         n_q[1] = (pout[2] - pout[0]) * (pin[2] - pin[3]) / ((flow[2] - flow[0]) * (pout[2] - pout[3])) - (pin[2] - pin[0]) / (flow[2] - flow[0]);
         rhs[1] = pin[2] + n_pout[1] * pout[2] + n_q[1] * flow[2];
         if ( SCIPisFeasLE(scip, pin[1] + n_pout[1] * pout[1] + n_q[1] * flow[1], rhs[1]) && SCIPisFeasLE(scip, pin[4] + n_pout[1] * pout[4] + n_q[1] * flow[4], rhs[1]) )
         {
            /* determine the coefficients and rhs for the triangle 0-3-4 */
            n_pout[2] = (flow[3] - flow[0]) * (pin[4] - pin[0]) / ((pout[3] -  pout[0]) * (flow[4] - flow[0])) - (pin[3] - pin[0]) / (pout[3] - pout[0]);
            n_q[2]    = - (pin[4] - pin[0]) / (flow[4] - flow[0]);
            rhs[2]    = pin[0] + pout[0] * n_pout[2] + flow[0] * n_q[2];
            if ( SCIPisFeasLE(scip, pin[1] + n_pout[2] * pout[1] + n_q[2] * flow[1], rhs[2]) && SCIPisFeasLE(scip, pin[2] + n_pout[2] * pout[2] + n_q[2] * flow[2], rhs[2]) )
            {
               infeasible = FALSE;
            }
         }

         if ( infeasible )
         {
            /* determine the coefficients and rhs for the triangle 2-3-4 */
            n_pout[1] = - (pin[2] - pin[3]) / (pout[2] - pout[3]);
            n_q[1] = (pout[2] - pout[4]) * (pin[2] - pin[3]) / ((flow[2] - flow[4]) * (pout[2] - pout[3])) - (pin[2] - pin[4]) / (flow[2] - flow[4]);
            rhs[1] = pin[2] + n_pout[1] * pout[2] + n_q[1] * flow[2];
            if ( SCIPisFeasLE(scip, pin[0] + n_pout[1] * pout[0] + n_q[1] * flow[0], rhs[1]) && SCIPisFeasLE(scip, pin[1] + n_pout[1] * pout[1] + n_q[1] * flow[1], rhs[1]) )
            {
               /* determine the coefficients and rhs for the triangle 0-2-4 */
               n_pout[2] = (flow[2] - flow[0]) * (pin[4] - pin[0]) / ((pout[2] -  pout[0]) * (flow[4] - flow[0])) - (pin[2] - pin[0]) / (pout[2] - pout[0]);
               n_q[2] = - (pin[4] - pin[0]) / (flow[4] - flow[0]);
               rhs[2] = pin[0] + pout[0] * n_pout[2] + flow[0] * n_q[2];
               if ( SCIPisFeasLE(scip, pin[1] + n_pout[2] * pout[1] + n_q[2] * flow[1], rhs[2]) && SCIPisFeasLE(scip, pin[3] + n_pout[2] * pout[3] + n_q[2] * flow[3], rhs[2]) )
               {
                  infeasible = FALSE;
               }
            }
         }
      }

      if ( infeasible )
      {
         /* determine the coefficients and rhs for the triangle 1-2-3 */
         n_pout[0] = - (pin[2] - pin[3]) / (pout[2] - pout[3]);
         n_q[0]    = - (pin[2] - pin[1]) / (flow[2] - flow[1]);
         rhs[0]    = pin[2] + pout[2] * n_pout[0] + flow[2] * n_q[0];
         if ( SCIPisFeasLE(scip, pin[0] + n_pout[0] * pout[0] + n_q[0] * flow[0], rhs[0]) && SCIPisFeasLE(scip, pin[4] + n_pout[0] * pout[4] + n_q[0] * flow[4], rhs[0]) )
         {
            /* determine the coefficients and rhs for the triangle 0-1-3 */
            n_pout[1] = - (pin[1] - pin[0]) / (pout[1] - pout[0]);
            n_q[1] = (pout[3] - pout[0]) * (pin[1] - pin[0]) / ((flow[3] - flow[0]) * (pout[1] - pout[0])) - (pin[3] - pin[0]) / (flow[3] - flow[0]);
            rhs[1] = pin[0] + n_pout[1] * pout[0] + n_q[1] * flow[0];

            if ( SCIPisFeasLE(scip, pin[2] + n_pout[1] * pout[2] + n_q[1] * flow[2], rhs[1]) && SCIPisFeasLE(scip, pin[4] + n_pout[1] * pout[4] + n_q[1] * flow[4], rhs[1]) )
            {
               /* determine the coefficients and rhs for the triangle 0-3-4 */
               n_pout[2] = (flow[3] - flow[0]) * (pin[4] - pin[0]) / ((pout[3] -  pout[0]) * (flow[4] - flow[0])) - (pin[3] - pin[0]) / (pout[3] - pout[0]);
               n_q[2]    = - (pin[4] - pin[0]) / (flow[4] - flow[0]);
               rhs[2]    = pin[0] + pout[0] * n_pout[2] + flow[0] * n_q[2];
               if ( SCIPisFeasLE(scip, pin[1] + n_pout[2] * pout[1] + n_q[2] * flow[1], rhs[2]) && SCIPisFeasLE(scip, pin[2] + n_pout[2] * pout[2] + n_q[2] * flow[2], rhs[2]) )
               {
                  infeasible = FALSE;
               }
            }

            if ( infeasible )
            {
               /* determine the coefficients and rhs for the triangle 0-1-4 */
               n_pout[1] = - (pin[1] - pin[0]) / (pout[1] - pout[0]);
               n_q[1]    = - (pin[4] - pin[0]) / (flow[4] - flow[0]);
               rhs[1]    = pin[0] + pout[0] * n_pout[1] + flow[0] * n_q[1];
               if ( SCIPisFeasLE(scip, pin[2] + n_pout[1] * pout[2] + n_q[1] * flow[2], rhs[1]) && SCIPisFeasLE(scip, pin[3] + n_pout[1] * pout[3] + n_q[1] * flow[3], rhs[1]) )
               {
                  SCIP_Real tmp;
                  /* determine the coefficients and rhs for the triangle 1-3-4 */
                  tmp = (pout[3] - pout[1]) * (flow[4] - flow[1]) - (pout[4] - pout[1]) * (flow[3] - flow[1]);
                  n_pout[2] = ((flow[3] - flow[1]) * (pin[4] - pin[1]) - (flow[4] - flow[1]) * (pin[3] - pin[1])) / tmp;
                  n_q[2] = ((pout[4] - pout[1]) * (pin[3] - pin[1]) - (pout[3] - pout[1]) * (pin[4] - pin[1])) / tmp;
                  rhs[2] = pin[1] + pout[1] * n_pout[2] + flow[1] * n_q[2];
                  if ( SCIPisFeasLE(scip, pin[0] + n_pout[2] * pout[0] + n_q[2] * flow[0], rhs[2]) && SCIPisFeasLE(scip, pin[2] + n_pout[2] * pout[2] + n_q[2] * flow[2], rhs[2]) )
                  {
                     infeasible = FALSE;
                  }
               }
            }
         }
      }

      if ( infeasible )
      {
         /* determine the coefficients and rhs for the triangle 0-1-4 */
         n_pout[0] = - (pin[1] - pin[0]) / (pout[1] - pout[0]);
         n_q[0]    = - (pin[4] - pin[0]) / (flow[4] - flow[0]);
         rhs[0]    = pin[0] + pout[0] * n_pout[0] + flow[0] * n_q[0];
         if ( SCIPisFeasLE(scip, pin[2] + n_pout[0] * pout[2] + n_q[0] * flow[2], rhs[0]) && SCIPisFeasLE(scip, pin[3] + n_pout[0] * pout[3] + n_q[0] * flow[3], rhs[0]) )
         {
            /* determine the coefficients and rhs for the triangle 2-3-4 */
            n_pout[1] = - (pin[2] - pin[3]) / (pout[2] - pout[3]);
            n_q[1] = (pout[2] - pout[4]) * (pin[2] - pin[3]) / ((flow[2] - flow[4]) * (pout[2] - pout[3])) - (pin[2] - pin[4]) / (flow[2] - flow[4]);
            rhs[1] = pin[2] + n_pout[1] * pout[2] + n_q[1] * flow[2];
            if ( SCIPisFeasLE(scip, pin[0] + n_pout[1] * pout[0] + n_q[1] * flow[0], rhs[1]) && SCIPisFeasLE(scip, pin[1] + n_pout[1] * pout[1] + n_q[1] * flow[1], rhs[1]) )
            {
               /* determine the coefficients and rhs for the triangle 1-2-4 */
               n_pout[2] = - (flow[4] - flow[1]) * (pin[2] - pin[1]) / ((pout[1] -  pout[4]) * (flow[2] - flow[1])) + (pin[4] - pin[1]) / (pout[2] - pout[1]);
               n_q[2]    = (pin[2] - pin[1]) / (flow[2] - flow[1]);
               rhs[2]    = pin[1] + pout[1] * n_pout[2] + flow[1] * n_q[2];
               if ( SCIPisFeasLE(scip, pin[0] + n_pout[2] * pout[0] + n_q[2] * flow[0], rhs[2]) && SCIPisFeasLE(scip, pin[3] + n_pout[2] * pout[3] + n_q[2] * flow[3], rhs[2]) )
               {
                  infeasible = FALSE;
               }
            }
         }
      }

      if ( infeasible )
      {
         i = 0;
         SCIPerrorMessage("In generateCaveNothingFixed: no overerstimator found for Case 5.\n");
      }
   }
   else
   {
      SCIPerrorMessage("Invalid shape for generateCaveNothingFixed.\n");
      return SCIP_INVALIDDATA;
   }

   if ( infeasible )
      return SCIP_OKAY;

   /* add offset to rhs */
   for ( j = 0; j < i; ++j )
   {
      rhs[j] += conshdlrdata->offset;
   }

   /* check if solution is cut off */
   if ( p_in_sol != NULL )
   {
      addcut = FALSE;
      for ( j = 0; j < i; ++j )
      {
         assert( p_out_sol != NULL );
         assert( q_sol != NULL );
         if ( dir == 1 )
         {
            if( SCIPisFeasGT(scip, n_pout[j] * (*p_out_sol) + (*p_in_sol) + n_q[j] * ABS(*q_sol), rhs[j]) )
            {
               addcut = TRUE;
               break;
            }
         }
         else
         {
            if( SCIPisFeasGT(scip, n_pout[j] * (*p_in_sol) + (*p_out_sol) + n_q[j] * ABS(*q_sol), rhs[j]) )
            {
               addcut = TRUE;
               break;
            }
         }
      }
   }

   if ( addcut )
   {
      while ( i > 0 )
      {
         --i;
         /* create constraint name */
         (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "cave_%s#%d", SCIPconsGetName(cons), consdata->nr_cave+1);

#if SCIP_VERSION >= 700
         SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &(cave[i]), conshdlr, s, -SCIPinfinity(scip), rhs[i], TRUE, TRUE, FALSE) );
#elif (SCIP_VERSION >= 602 && SCIP_SUBVERSION > 0)
         SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &(cave[i]), conshdlr, s, -SCIPinfinity(scip), rhs[i], TRUE, TRUE, FALSE) );
#else
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, &(cave[i]), conshdlr, s, -SCIPinfinity(scip), rhs[i], TRUE, TRUE, FALSE) );
#endif
         SCIP_CALL( SCIPcacheRowExtensions(scip, cave[i]) );

         SCIP_CALL( SCIPaddVarToRow(scip, cave[i], consdata->q_var, dir * n_q[i]) );
         if ( dir == 1 )
         {
            SCIP_CALL( SCIPaddVarToRow(scip, cave[i], consdata->p_out_var, n_pout[i]) );
            SCIP_CALL( SCIPaddVarToRow(scip, cave[i], consdata->p_in_var, 1.0) );
         }
         else
         {
            SCIP_CALL( SCIPaddVarToRow(scip, cave[i], consdata->p_out_var, 1.0) );
            SCIP_CALL( SCIPaddVarToRow(scip, cave[i], consdata->p_in_var, n_pout[i]) );
         }

         SCIP_CALL( SCIPflushRowExtensions(scip, cave[i]) );

#if SCIP_VERSION >= 500
         SCIP_CALL( SCIPaddRow(scip, cave[i], FALSE, &infeasible) );
#else
         SCIP_CALL( SCIPaddCut(scip, NULL, cave[i], FALSE, &infeasible) );
#endif

#ifdef OVERESTIMATOR_DEBUG
         if ( SCIPgetSubscipDepth(scip) == 0 )
         {
            printf("<%s>: p_out: LB %f \t UB %f\n", SCIPconsGetName(cons), p_var_lb, p_var_ub);
            printf("<%s>: q: LB %f \t UB %f\n", SCIPconsGetName(cons), q_var_lb, q_var_ub);
            SCIP_CALL( SCIPprintRow(scip, cave[i], NULL) );
         }
#endif
         assert( ! infeasible );
         SCIP_CALL( SCIPreleaseRow(scip, &cave[i]) );

         ++(consdata->nr_cave);
      }
   }

   return SCIP_OKAY;
}

/** function to generate the linear constraints for the concave overestimator of the pressure-input function
 *  cut gets added, if p_in_sol and p_out_sol are NULL or the (pseudo)solution gets cut off */
static
SCIP_RETCODE generateCaveFixedFlow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            p_in_sol,           /**< pointer to the (pseudo)solution of input pressure */
   SCIP_Real*            p_out_sol           /**< pointer to the (pseudo)solution of output pressure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_ROW* cave;
   SCIP_Real q_var_ub;
   SCIP_Real p_var_ub;
   SCIP_Real p_var_lb;
   SCIP_Real rhs;
   SCIP_Real pin_u;
   SCIP_Real pin_l;
   SCIP_Real n_out;
   SCIP_Real p_bound;
   int dir;
   SCIP_Bool addcut = TRUE;
   SCIP_Bool infeasible = FALSE;
   char s[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   if ( p_in_sol != NULL || p_out_sol != NULL )
   {
      assert( p_in_sol != NULL );
      assert( p_out_sol != NULL );
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   q_var_ub = SCIPcomputeVarUbLocal(scip, consdata->q_var);

   if ( SCIPisFeasGE(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var), 0.0) )
   {
      p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_out_var);
      p_var_lb = SCIPcomputeVarLbLocal(scip, consdata->p_out_var);
      dir = 1;
   }
   else if ( SCIPisFeasLE(scip, q_var_ub, 0.0) )
   {
      q_var_ub = - SCIPcomputeVarLbLocal(scip, consdata->q_var);
      p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_in_var);
      p_var_lb = SCIPcomputeVarLbLocal(scip, consdata->p_in_var);
      dir = -1;
   }
   else
   {
      SCIPdebugMessage("generateCaveFixedFlow works only, when the direction of the flow is fixed.\n");
      return SCIP_INVALIDDATA;
   }

   p_bound = (5.0 * consdata->c * q_var_ub)/(((SCIP_Real) conshdlrdata->machbound)  * 1e5 * consdata->A);
   if( SCIPisFeasGT(scip, p_bound, p_var_lb) )
   {
      p_var_lb = p_bound;
   }

   if ( SCIPisFeasGT(scip, p_var_lb, p_var_ub) )
   {
      SCIPerrorMessage("generateCaveFixedFlow: Upper pressure bound does not satisfy the mach bound.\n");
      return SCIP_ERROR;
   }

   pin_u = computeTrapezoidalScheme(scip, conshdlrdata, consdata, p_var_ub, q_var_ub, -1 );
   pin_l = computeTrapezoidalScheme(scip, conshdlrdata, consdata, p_var_lb, q_var_ub, -1 );
   if ( SCIPisNegative(scip, pin_u) || SCIPisNegative(scip, pin_l) )
   {
      SCIPerrorMessage("Something went wrong when performing computeTrapezoidalScheme.\n");
      return SCIP_ERROR;
   }
   rhs = (pin_l * p_var_ub - pin_u * p_var_lb) / (p_var_ub - p_var_lb);

   n_out = (pin_l - pin_u) / (p_var_ub - p_var_lb);

   /* add offset to rhs */
   rhs += conshdlrdata->offset;

   if ( (p_in_sol != NULL) && (p_out_sol != NULL) )
      if ( SCIPisFeasLE(scip, n_out * (*p_out_sol) + (*p_in_sol), rhs) )
         addcut = FALSE;

   if ( addcut )
   {
      /* create constraint name */
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "cave_%s#%d", SCIPconsGetName(cons), consdata->nr_cave+1);

#if SCIP_VERSION >= 700
      SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &cave, conshdlr, s, -SCIPinfinity(scip), rhs, TRUE, TRUE, FALSE) );
#elif (SCIP_VERSION >= 602 && SCIP_SUBVERSION > 0)
      SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &cave, conshdlr, s, -SCIPinfinity(scip), rhs, TRUE, TRUE, FALSE) );
#else
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cave, conshdlr, s, -SCIPinfinity(scip), rhs, TRUE, TRUE, FALSE) );
#endif
      SCIP_CALL( SCIPcacheRowExtensions(scip, cave) );

      if ( dir == 1 )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cave, consdata->p_out_var, n_out) );
         SCIP_CALL( SCIPaddVarToRow(scip, cave, consdata->p_in_var, 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cave, consdata->p_out_var, 1.0) );
         SCIP_CALL( SCIPaddVarToRow(scip, cave, consdata->p_in_var, n_out) );
      }

      SCIP_CALL( SCIPflushRowExtensions(scip, cave) );

#if SCIP_VERSION >= 500
      SCIP_CALL( SCIPaddRow(scip, cave, TRUE, &infeasible) );
#else
      SCIP_CALL( SCIPaddCut(scip, NULL, cave, TRUE, &infeasible) );
#endif

#ifdef OVERESTIMATOR_DEBUG
      if ( SCIPgetSubscipDepth(scip) == 0 )
      {
         printf("<%s>: p_out: LB %f \t UB %f\n", SCIPconsGetName(cons), p_var_lb, p_var_ub);
         printf("<%s>: q: LB %f \t UB %f\n", SCIPconsGetName(cons), q_var_lb, q_var_ub);
         SCIP_CALL( SCIPprintRow(scip, cave, NULL) );
      }
#endif

      assert( ! infeasible );
      SCIP_CALL( SCIPreleaseRow(scip, &cave) );
      ++(consdata->nr_cave);
   }

   return SCIP_OKAY;
}

/** function to generate the linear constraints for the concave overestimator of the pressure-input function
 *  cut gets added, if p_in_sol and p_out_sol are NULL or the (pseudo)solution gets cut off */
static
SCIP_RETCODE generateCaveFixedOutput(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            p_in_sol,           /**< pointer to the (pseudo)solution of input pressure */
   SCIP_Real*            p_out_sol,          /**< pointer to the (pseudo)solution of output pressure */
   SCIP_Real*            q_sol               /**< pointer to the (pseudo)solution of the flow */
   )
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_ROW*      cave;
   SCIP_Real      q_var_ub;
   SCIP_Real      q_var_lb;
   SCIP_Real      p_var_ub;
   SCIP_Real      rhs;
   SCIP_Real      pin_u;
   SCIP_Real      pin_l;
   SCIP_Real      n_q;
   SCIP_Real      q_bound;
   int            dir;
   SCIP_Bool      addcut     = TRUE;
   SCIP_Bool      infeasible = FALSE;
   char           s[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   if ( p_in_sol != NULL || q_sol != NULL )
   {
      assert( p_in_sol != NULL );
      assert( q_sol != NULL );
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   q_var_lb = SCIPcomputeVarLbLocal(scip, consdata->q_var);
   q_var_ub = SCIPcomputeVarUbLocal(scip, consdata->q_var);

   if ( SCIPisFeasGE(scip, q_var_lb, 0.0) )
   {
      p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_out_var);
      dir = 1;
   }
   else if ( SCIPisFeasLE(scip, q_var_ub, 0.0) )
   {
      q_var_ub = - SCIPcomputeVarLbLocal(scip, consdata->q_var);
      q_var_lb = - SCIPcomputeVarUbLocal(scip, consdata->q_var);
      p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_in_var);
      dir = -1;
   }
   else
   {
      SCIPdebugMessage("generateCaveFixedOutput works only, when the direction of the flow is fixed.\n");
      return SCIP_INVALIDDATA;
   }

   q_bound = (((SCIP_Real) conshdlrdata->machbound) * consdata->A * 1e5 * p_var_ub) / (5.0 * consdata->c);
   if ( SCIPisFeasLT(scip, q_bound, q_var_ub) )
   {
      q_var_ub = q_bound;
   }

   if ( SCIPisFeasGT(scip, q_var_lb, q_var_ub) )
   {
      SCIPerrorMessage("generateCaveFixedFlow: Upper flow bound does not satisfy the mach bound.\n");
      return SCIP_ERROR;
   }

   pin_u = computeTrapezoidalScheme(scip, conshdlrdata, consdata, p_var_ub, q_var_ub, -1 );
   pin_l = computeTrapezoidalScheme(scip, conshdlrdata, consdata, p_var_ub, q_var_lb, -1 );
   if ( SCIPisNegative(scip, pin_u) || SCIPisNegative(scip, pin_l) )
   {
      SCIPerrorMessage("Something went wrong when performing computeTrapezoidalScheme.\n");
      return SCIP_ERROR;
   }

   rhs = (- pin_u * q_var_lb + pin_l * q_var_ub) / (q_var_ub - q_var_lb);

   n_q = (pin_l - pin_u) * ((SCIP_Real) dir) / (q_var_ub - q_var_lb);

   /* add offset to rhs */
   rhs += conshdlrdata->offset;

   if ( (p_in_sol != NULL) && (q_sol != NULL) )
   {
      if ( SCIPisFeasLE(scip, n_q * (*q_sol) + (*p_in_sol), rhs) )
         addcut = FALSE;
   }

   if ( addcut )
   {
      /* create constraint name */
      (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "cave_%s#%d", SCIPconsGetName(cons), consdata->nr_cave+1);

#if SCIP_VERSION >= 700
      SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &cave, conshdlr, s, -SCIPinfinity(scip), rhs, TRUE, TRUE, FALSE) );
#elif (SCIP_VERSION >= 602 && SCIP_SUBVERSION > 0)
      SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &cave, conshdlr, s, -SCIPinfinity(scip), rhs, TRUE, TRUE, FALSE) );
#else
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &cave, conshdlr, s, -SCIPinfinity(scip), rhs, TRUE, TRUE, FALSE) );
#endif
      SCIP_CALL( SCIPcacheRowExtensions(scip, cave) );

      SCIP_CALL( SCIPaddVarToRow(scip, cave, consdata->q_var, n_q) );
      if ( dir == 1 )
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cave, consdata->p_in_var, 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPaddVarToRow(scip, cave, consdata->p_out_var, 1.0) );
      }

      SCIP_CALL( SCIPflushRowExtensions(scip, cave) );

#if SCIP_VERSION >= 500
      SCIP_CALL( SCIPaddRow(scip, cave, FALSE, &infeasible) );
#else
      SCIP_CALL( SCIPaddCut(scip, NULL, cave, FALSE, &infeasible) );
#endif

#ifdef OVERESTIMATOR_DEBUG
      if ( SCIPgetSubscipDepth(scip) == 0 )
      {
         printf("<%s>: p_out: LB %f \t UB %f\n", SCIPconsGetName(cons), p_var_lb, p_var_ub);
         printf("<%s>: q: LB %f \t UB %f\n", SCIPconsGetName(cons), q_var_lb, q_var_ub);
         SCIP_CALL( SCIPprintRow(scip, cave, NULL) );
      }
#endif

      assert( ! infeasible );
      SCIP_CALL( SCIPreleaseRow(scip, &cave) );
      ++(consdata->nr_cave);
   }

   return SCIP_OKAY;
}

/** function to decide which generateCave*-function to call
 *  returns the number of rows added to the lp */
static
SCIP_RETCODE generateCave(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            p_in_sol,           /**< pointer to the (pseudo)solution of input pressure */
   SCIP_Real*            p_out_sol,          /**< pointer to the (pseudo)solution of output pressure */
   SCIP_Real*            q_sol,              /**< pointer to the (pseudo)solution of the flow */
   int*                  counter             /**< counter for the number of generated rows or NULL */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_Real          q_var_ub;
   SCIP_Real          q_var_lb;
   SCIP_Real          p_var_ub;
   SCIP_Real          p_var_lb;
   int                nr_createdrows;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nr_createdrows = - consdata->nr_cave;

   q_var_lb = SCIPcomputeVarLbLocal(scip, consdata->q_var);
   q_var_ub = SCIPcomputeVarUbLocal(scip, consdata->q_var);

   if ( SCIPisFeasGE(scip, q_var_lb, 0.0) )
   {
      p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_out_var);
      p_var_lb = SCIPcomputeVarLbLocal(scip, consdata->p_out_var);
   }
   else if ( SCIPisFeasLE(scip, q_var_ub, 0.0)  )
   {
      p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_in_var);
      p_var_lb = SCIPcomputeVarLbLocal(scip, consdata->p_in_var);
   }
   else
   {
      SCIPdebugMessage("Cannot generate an overestimator if the orientation of the flow is no fixed.\n");
      return SCIP_OKAY;
   }

   if ( pressureIsEQ(p_var_lb, p_var_ub, conshdlrdata->pressure_feastol) && flowIsEQ(q_var_ub, q_var_lb, conshdlrdata->flow_feastol) )
   {
      SCIPdebugMessage("Cannot generate an overestimator if both output and flow are FeasEQ.\n");
      return SCIP_OKAY;
   }
   else if ( pressureIsEQ(p_var_lb, p_var_ub, conshdlrdata->pressure_feastol) )
   {
      SCIP_CALL( generateCaveFixedOutput(scip, conshdlr, cons, p_in_sol, p_out_sol, q_sol) );
   }
   else if ( flowIsEQ(q_var_ub, q_var_lb, conshdlrdata->flow_feastol) )
   {
      SCIP_CALL( generateCaveFixedFlow(scip, conshdlr, cons, p_in_sol, p_out_sol) );
   }
   else
   {
      SCIP_CALL( generateCaveNothingFixed(scip, conshdlr, cons, p_in_sol, p_out_sol, q_sol) );
   }

   nr_createdrows += consdata->nr_cave;
   if ( nr_createdrows > 0 )
   {
      SCIPdebugMessage("Added %d linear cuts to the pipe %s.\n", nr_createdrows, SCIPconsGetName(cons));
#ifdef OVERESTIMATOR_DEBUG_MIN_OUTPUT
      if ( SCIPgetSubscipDepth(scip) == 0 )
         printf("Added %d linear cuts to the pipe <%s>.\n", nr_createdrows, SCIPconsGetName(cons));
#elif defined OVERESTIMATOR_DEBUG
      if ( SCIPgetSubscipDepth(scip) == 0 )
         printf("Added %d linear cuts to the pipe <%s>.\n", nr_createdrows, SCIPconsGetName(cons));
#endif
   }

   if ( counter != NULL )
   {
      *counter = nr_createdrows;
   }

   return SCIP_OKAY;
}

/** function that creates a gradient cut for the given values of p_out and flow.
 *  If p_in_sol is NULL or separated by the generate cut, then the cut is added to the lp.
 *  @note: this function should only be called if the orientation of
 *  the flow is unique. */
static
SCIP_RETCODE addGradientCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_Real*            p_in_sol,           /**< solution value of input pressure or NULL */
   SCIP_Real             p_out,              /**< output pressure value to generate cut with, should be p_out_sol if p_in_sol != NULL  */
   SCIP_Real             flow,               /**< flow value, should be q_sol if p_in_sol != NULL */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_ROW* vex;
   SCIP_Real dp;
   SCIP_Real dq;
   SCIP_Real lhs;
   int dir;
   SCIP_Bool infeasible;
   char s[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( result != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( SCIPisFeasGE(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var), 0.0) )
      dir = 1;
   else
   {
      assert( SCIPisFeasLE(scip, SCIPcomputeVarUbLocal(scip, consdata->q_var), 0.0) );
      dir = -1;
   }

   assert ( SCIPisFeasGE(scip, (((SCIP_Real) conshdlrdata->machbound) * consdata->A * 1e5 * p_out), (5.0 * consdata->c * ABS(flow))) );
   lhs = computeModifiedEulerSchemeWithGradients(scip, conshdlrdata, consdata, p_out, ABS(flow), &dp, &dq);
   /* shift the cut by conshdlrdata->offset away from feasible set  */
   lhs -= conshdlrdata->offset;

   if ( p_in_sol != NULL )
   {
      if ( SCIPisFeasGE(scip, *p_in_sol, lhs) )
      {
         /* cut does not separate solution */
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }
   }

   lhs -= (dp * p_out + dq * ABS(flow));

   /* create constraint name */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "vex_%s%d", SCIPconsGetName(cons), consdata->nr_vex + 1);

#if SCIP_VERSION >= 700
   SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &vex, conshdlr, s, lhs, SCIPinfinity(scip), TRUE, TRUE, FALSE) );
#elif (SCIP_VERSION >= 602 && SCIP_SUBVERSION > 0)
   SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &vex, conshdlr, s, lhs, SCIPinfinity(scip), TRUE, TRUE, FALSE) );
#else
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &vex, conshdlr, s, lhs, SCIPinfinity(scip), TRUE, TRUE, FALSE) );
#endif
   SCIP_CALL( SCIPcacheRowExtensions(scip, vex) );

   if ( dir == 1 )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, vex, consdata->p_in_var, 1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, vex, consdata->q_var, - dq) );
      SCIP_CALL( SCIPaddVarToRow(scip, vex, consdata->p_out_var, - dp) );
   }
   else
   {
      SCIP_CALL( SCIPaddVarToRow(scip, vex, consdata->p_out_var, 1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, vex, consdata->q_var, dq) );
      SCIP_CALL( SCIPaddVarToRow(scip, vex, consdata->p_in_var, - dp) );
   }

   SCIP_CALL( SCIPflushRowExtensions(scip, vex) );

#if SCIP_VERSION >= 500
   SCIP_CALL( SCIPaddRow(scip, vex, FALSE, &infeasible) );
#else
   SCIP_CALL( SCIPaddCut(scip, NULL, vex, FALSE, &infeasible) );
#endif

#ifdef UNDERESTIMATOR_DEBUG_MIN_OUTPUT
   if ( SCIPgetSubscipDepth(scip) == 0 )
      printf("Add gradient cut for pipe %s\n", SCIPconsGetName(cons));
#elif defined UNDERESTIMATOR_DEBUG
   if ( SCIPgetSubscipDepth(scip) == 0 )
   {
      printf("Add gradient cut for pipe %s\n", SCIPconsGetName(cons));
      if ( dir == 1 )
      {
         printf("<%s>: p_out: LB %f \t UB %f\n", SCIPconsGetName(cons), SCIPcomputeVarLbLocal(scip, consdata->p_out_var),
            SCIPcomputeVarUbLocal(scip, consdata->p_out_var));
      }
      else
      {
         printf("<%s>: p_out: LB %f \t UB %f\n", SCIPconsGetName(cons), SCIPcomputeVarLbLocal(scip, consdata->p_in_var),
            SCIPcomputeVarUbLocal(scip, consdata->p_in_var));
      }
      printf("<%s>: q: LB %f \t UB %f\n", SCIPconsGetName(cons), SCIPcomputeVarLbLocal(scip, consdata->q_var),SCIPcomputeVarUbLocal(scip, consdata->q_var));
      SCIP_CALL( SCIPprintRow(scip, vex, NULL) );
   }
#endif

   *result = SCIP_SEPARATED;
   ++conshdlrdata->ngradientcuts;
   assert( SCIPgetRowLPFeasibility(scip, vex) < 0 );
   assert( ! infeasible );
   SCIP_CALL( SCIPreleaseRow(scip, &vex) );
   ++(consdata->nr_vex);

   return SCIP_OKAY;
}

/** this functions checks wheter a given solution is feasible for one ode constraint or not */
static
SCIP_RETCODE checkSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< pipe ode constraint handler */
   SCIP_CONS*            cons,               /**< ode constraint */
   SCIP_SOL*             sol,                /**< current pseudo, lp, etc solution or NULL */
   SCIP_RESULT*          result              /**< pointer to store the result, either feasible or infeasible */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real q_sol;
   SCIP_Real p_in_sol;
   SCIP_Real p_in_ub;
   SCIP_Real p_in_lb;
   SCIP_Real p_out_sol;
   SCIP_Real numerator;

   assert( cons != NULL );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   /* get data of constraint handler*/
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);

   /* read solution values */
   q_sol = SCIPgetSolVal(scip, sol, consdata->q_var);
   if ( SCIPisFeasGE(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var), 0.0) || SCIPisFeasPositive(scip, q_sol) )
   {
      p_in_sol = SCIPgetSolVal(scip, sol, consdata->p_in_var);
      p_out_sol = SCIPgetSolVal(scip, sol, consdata->p_out_var);
   }
   else
   {
      p_out_sol = SCIPgetSolVal(scip, sol, consdata->p_in_var);
      p_in_sol = SCIPgetSolVal(scip, sol, consdata->p_out_var);
   }

   numerator = (SCIP_Real) conshdlrdata->machbound;

   /* check "mach bound"
    * unnecessary for lp solutions, but has to be checked for pseudo solutions */
   if ( SCIPisFeasLT(scip, numerator * consdata->A * p_out_sol * 1e5, 5.0 * consdata->c * ABS(q_sol)) )
   {
      *result = SCIP_INFEASIBLE;
      consdata->feasible = FALSE;
      consdata->violation =(5.0 * consdata->c * ABS(q_sol)) / (numerator * consdata->A * 1e5) - p_out_sol;
      consdata->type_violation = MACHBOUND;
      return SCIP_OKAY;
   }

   /* compute relaxation */
   p_in_ub = computeTrapezoidalScheme(scip, conshdlrdata, consdata, p_out_sol, ABS(q_sol), -1);
   p_in_lb = computeModifiedEulerScheme(scip, conshdlrdata, consdata, p_out_sol, ABS(q_sol), -1);

   assert( pressureIsGE(p_in_ub, p_in_lb, conshdlrdata->pressure_feastol) );

   while ( ! pressureIsEQ(p_in_ub, p_in_lb, conshdlrdata->pressure_feastol) )
   {
      consdata->N *= 2;
#ifdef STEPSIZE_DEBUG
      printf("Increase N of pipe %s to %d.\n", SCIPconsGetName(cons), consdata->N);
#endif
      p_in_ub = computeTrapezoidalScheme(scip, conshdlrdata, consdata, p_out_sol, ABS(q_sol), -1);
      p_in_lb = computeModifiedEulerScheme(scip, conshdlrdata, consdata, p_out_sol, ABS(q_sol), -1);
   }

   assert( pressureIsGE(p_in_ub, p_in_lb, conshdlrdata->pressure_feastol) );

   if ( pressureIsGT(p_in_sol, p_in_ub, conshdlrdata->pressure_feastol) )
   {
      *result = SCIP_INFEASIBLE;
      consdata->feasible = FALSE;
      consdata->violation = p_in_sol - p_in_ub;
      consdata->type_violation = PIN_TOO_BIG;
      return SCIP_OKAY;
   }
   if ( pressureIsLT(p_in_sol, p_in_lb, conshdlrdata->pressure_feastol) )
   {
      *result = SCIP_INFEASIBLE;
      consdata->feasible = FALSE;
      consdata->violation = p_in_lb - p_in_sol;
      consdata->type_violation = PIN_TOO_LOW;
      return SCIP_OKAY;
   }

   consdata->feasible = TRUE;
   consdata->violation = 0.0;
   consdata->type_violation = FEASIBLE;
   return SCIP_OKAY;
}

/** function that returns the position of the constraint which the current LP solution violates the most
 *  returns -1 if the lp solution is feasible */
static
SCIP_RETCODE computeMostViolatedConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< pipe ode constraint handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< current LP, Pseudo, etc. solution or NULL */
   int*                  violatedcons        /**< pointer to store the most violated constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   SCIP_RESULT result;
   SCIP_Real maxviolation = 0.0;
   int i;
   int vcons = -1;

   assert( scip != NULL );

   for ( i = 0 ; i < nconss ; ++i )
   {
      cons = conss[i];
      assert( cons != NULL );

      SCIP_CALL( checkSolution(scip, conshdlr, cons, sol, &result) );

      if ( result == SCIP_INFEASIBLE )
      {
         /* get data of constraint */
         consdata = SCIPconsGetData(cons);
         assert( consdata != NULL);

         if ( SCIPisGT(scip, consdata->violation, maxviolation) )
         {
            vcons = i;
            maxviolation = consdata->violation;
         }
      }
   }

   *violatedcons = vcons;

   return SCIP_OKAY;
}

/** separates a solution  */
static
SCIP_RETCODE separateODE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< current solution */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real      q_var_ub;
   SCIP_Real      q_var_lb;
   int            counter = 0;
   SCIP_Real      p_in_sol;
   SCIP_Real      p_out_sol;
   SCIP_Real      q_sol;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);

   assert( ! consdata->feasible );
   *result = SCIP_DIDNOTFIND;

   q_var_ub = SCIPcomputeVarUbLocal(scip, consdata->q_var);
   q_var_lb = SCIPcomputeVarLbLocal(scip, consdata->q_var);

   /* our first priority is to fix the flow direction */
   if ( SCIPisFeasLT(scip, q_var_lb, 0.0) && SCIPisFeasGT(scip, q_var_ub, 0.0) )
      return SCIP_OKAY;

   /* read solution values */
   q_sol = SCIPgetSolVal(scip, sol, consdata->q_var);
   if ( SCIPisFeasGE(scip, q_var_lb, 0.0) )
   {
      p_in_sol = SCIPgetSolVal(scip, sol, consdata->p_in_var);
      p_out_sol = SCIPgetSolVal(scip, sol, consdata->p_out_var);
   }
   else if ( SCIPisFeasLE(scip, q_var_ub, 0.0) )
   {
      p_out_sol = SCIPgetSolVal(scip, sol, consdata->p_in_var);
      p_in_sol = SCIPgetSolVal(scip, sol, consdata->p_out_var);
   }
   else
   {
      SCIPerrorMessage("Flow direction should already be fixed if separateODE is called.\n");
      return SCIP_ERROR;
   }

   if ( consdata->type_violation == PIN_TOO_LOW )
   {
      SCIPdebugMessage("Add Gradient Cut.\n");
      /*lint -e{644}*/
      SCIP_CALL( addGradientCut(scip, conshdlr, cons, &p_in_sol, p_out_sol, q_sol, result) );
      return SCIP_OKAY;
   }

   if ( consdata->type_violation == PIN_TOO_BIG )
   {
      counter = 0;
      /* assert( counter == 0 ); */
      /* Try to cut off the solution via the concave overestimator, */
      SCIP_CALL( generateCave(scip, conshdlr, cons, &p_in_sol, &p_out_sol, &q_sol, &counter) );
      if ( counter > 0 )
      {
         SCIPdebugMessage("Add %d concave overestimator to the lp.\n",counter);
         *result = SCIP_SEPARATED;
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}

/** enforce ODE constraints
 *
 *  The function enforceODE assumes that all ode-constraints have recently been checked for feasibility.  enforceODE
 *  tries to resolve an infeasibility first by fixing the orientation of the flow, then bound propagation.  depending on
 *  the error type it tries to add an overestimator or a gradient cut. If nothing works, then extern branching
 *  candidates are added.
 */
static
SCIP_RETCODE enforceODE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< current solution */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real      q_var_ub;
   SCIP_Real      q_var_lb;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   /* get data of constraint */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL);

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);

   assert( ! consdata->feasible );
   *result = SCIP_INFEASIBLE;

   q_var_ub = SCIPcomputeVarUbLocal(scip, consdata->q_var);
   q_var_lb = SCIPcomputeVarLbLocal(scip, consdata->q_var);

   /* our first priority is to fix the flow direction */
   if ( SCIPisFeasNegative(scip, q_var_lb) && SCIPisFeasPositive(scip, q_var_ub) )
   {
      SCIPdebugMessage("Fix orientation of flow.\n");

      if ( probdata->noFlowBinvars )
      {
         SCIP_CALL( SCIPaddExternBranchCand(scip, consdata->q_var, MIN( - q_var_lb, q_var_ub), 0.0) );
         assert( SCIPcontainsExternBranchCand(scip, consdata->q_var) );
      }
      else
      {
         SCIP_Bool infeasible = FALSE;
         SCIP_Bool tightened = FALSE;

         /* If solinfeasible is true, the constraints might not be completely enforced. Then a binary flow direction
          * might be fixed, but the flow direction is not fixed. We then simply fix the flow direction, because it is
          * cheap and might pay of while resolving the LP. */
         if ( SCIPisFeasEQ(scip, SCIPcomputeVarLbLocal(scip, consdata->posFlowBinvar), 1.0) )
         {
            SCIP_CALL( SCIPinferVarLbCons(scip, consdata->q_var, 0.0, cons, POS_LB, TRUE, &infeasible, &tightened) );
         }
         else if ( SCIPisFeasEQ(scip, SCIPcomputeVarUbLocal(scip, consdata->posFlowBinvar), 0.0) )
         {
            SCIP_CALL( SCIPinferVarUbCons(scip, consdata->q_var, 0.0, cons, POS_UB, TRUE, &infeasible, &tightened) );
         }

         if ( infeasible )
            return SCIP_ERROR;

         if ( tightened )
         {
            *result = SCIP_REDUCEDDOM;
            return SCIP_OKAY;
         }

         if ( SCIPisFeasEQ(scip, SCIPcomputeVarLbLocal(scip, consdata->negFlowBinvar), 1.0) )
         {
            SCIP_CALL( SCIPinferVarUbCons(scip, consdata->q_var, 0.0, cons, NEG_LB, TRUE, &infeasible, &tightened) );
         }
         else if ( SCIPisFeasEQ(scip, SCIPcomputeVarUbLocal(scip, consdata->negFlowBinvar), 0.0) )
         {
            SCIP_CALL( SCIPinferVarLbCons(scip, consdata->q_var, 0.0, cons, NEG_UB, TRUE, &infeasible, &tightened) );
         }

         if ( infeasible )
            return SCIP_ERROR;

         if ( tightened )
         {
            *result = SCIP_REDUCEDDOM;
            return SCIP_OKAY;
         }

         /* Only add posFlowBinvar or negFlowBinvar!
          * Adding both can result in an error, when one is fixed to 1.0 and
          * then branching on the other occurrs although it is implicitly fixed.
          */
         /* SCIP_CALL( SCIPaddExternBranchCand(scip, consdata->negFlowBinvar, 1.0, 0.0) ); */
         SCIP_CALL( SCIPaddExternBranchCand(scip, consdata->posFlowBinvar, 1.0, 0.0) );
      }

      return SCIP_OKAY;
   }

   SCIP_CALL( separateODE(scip, conshdlr, cons, sol, result) );

   /* if the infeasibility could not be resolved, we try branching */
   /* we branch on the middle value of the variable bounds */
   if ( *result == SCIP_DIDNOTFIND )
   {
      SCIP_Real p_in_ub;
      SCIP_Real p_in_lb;
      SCIP_Real p_out_ub;
      SCIP_Real p_out_lb;

      *result = SCIP_INFEASIBLE;
      p_in_lb = SCIPcomputeVarLbLocal(scip, consdata->p_in_var);
      p_in_ub = SCIPcomputeVarUbLocal(scip, consdata->p_in_var);
      p_out_lb = SCIPcomputeVarLbLocal(scip, consdata->p_out_var);
      p_out_ub = SCIPcomputeVarUbLocal(scip, consdata->p_out_var);

      /* SCIP_CALL( SCIPaddExternBranchCand(scip, variable, score, value) ); */
      if ( !SCIPisEQ(scip, q_var_ub, q_var_lb) )
      {
         SCIP_CALL( SCIPaddExternBranchCand(scip, consdata->q_var, (q_var_ub - q_var_lb) / (2 * MAX(q_var_ub, -q_var_lb)), q_var_ub / 2 + q_var_lb / 2) );
      }
      if ( !SCIPisEQ(scip, p_in_ub, p_in_lb) )
      {
         SCIP_CALL( SCIPaddExternBranchCand(scip, consdata->p_in_var, (p_in_ub - p_in_lb) / (2 * p_in_ub), p_in_lb / 2 + p_in_ub / 2) );
      }
      if ( !SCIPisEQ(scip, p_out_ub, p_out_lb) )
      {
         SCIP_CALL( SCIPaddExternBranchCand(scip, consdata->p_out_var, (p_out_ub - p_out_lb) / (2 * p_out_ub), p_out_lb / 2 + p_out_ub / 2) );
      }

      assert( SCIPcontainsExternBranchCand(scip, consdata->p_in_var) || SCIPcontainsExternBranchCand(scip, consdata->p_out_var) || SCIPcontainsExternBranchCand(scip, consdata->q_var) );
   }

   return SCIP_OKAY;
}

/** enforcing pseudo solution */
static
SCIP_RETCODE enforcePS(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< current solution */
   SCIP_RESULT*          result              /**< pointer to store the result */
   )
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* q_var;
   SCIP_VAR* pout_var;
   SCIP_VAR* pin_var;
   SCIP_Real q_sol, pin_sol, pout_sol;
   SCIP_Real q_ub, pin_ub, pout_ub;
   SCIP_Real q_lb, pin_lb, pout_lb;
   SCIP_Bool q_lower, pin_lower, pout_lower;
   SCIP_Real new_p_bound;
   SCIP_Bool infeasible, tightened;
   SCIP_Bool preferFlow = FALSE;
   int dir;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   q_var = consdata->q_var;
   q_lb = SCIPcomputeVarLbLocal(scip, q_var);
   q_ub = SCIPcomputeVarUbLocal(scip, q_var);

   /* our first priority is to fix the flow direction */
   if ( flowIsNegative(q_lb, conshdlrdata->flow_feastol) && flowIsPositive(q_ub, conshdlrdata->flow_feastol) )
   {
      SCIPdebugMessage("Fix orientation of flow.\n");
      SCIP_CALL( SCIPaddExternBranchCand(scip, q_var, consdata->violation, 0.0) );
      assert( SCIPcontainsExternBranchCand(scip, q_var) );
      return SCIP_OKAY;
   }

   /* set the correct values _sol, _lb, _ub */
   if ( flowIsPositive(q_ub, conshdlrdata->flow_feastol) )
   {
      pin_var  = consdata->p_in_var;
      pout_var = consdata->p_out_var;
      dir = 1;
   }
   else if ( flowIsNegative(q_lb, conshdlrdata->flow_feastol) )
   {
      pin_var  = consdata->p_out_var;
      pout_var = consdata->p_in_var;
      q_ub     = - q_lb;
      q_lb     = - SCIPcomputeVarUbLocal(scip, q_var);
      dir = -1;
   }
   else
   {
      /* flow is (basically) zero, pin should be equal to pout */
      pin_var  = consdata->p_in_var;
      pout_var = consdata->p_out_var;
      dir = 0;
   }

   pin_lb  = SCIPcomputeVarLbLocal(scip, pin_var);
   pin_ub  = SCIPcomputeVarUbLocal(scip, pin_var);
   pout_lb = SCIPcomputeVarLbLocal(scip, pout_var);
   pout_ub = SCIPcomputeVarUbLocal(scip, pout_var);

   pin_sol  = SCIPgetSolVal(scip, sol, pin_var);
   q_sol    = SCIPgetSolVal(scip, sol, pout_var);
   pout_sol = ABS(SCIPgetSolVal(scip, sol, q_var)); /*lint !e666*/

   if ( dir == 0 )
   {
      /* check if problem is infeasible */
      if ( pressureIsGT(pin_lb, pout_ub, conshdlrdata->pressure_feastol) || pressureIsGT(pout_lb, pin_ub, conshdlrdata->pressure_feastol) )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      /* check if a pressure bound can be improved */
      if ( SCIPisGT(scip, pin_lb, pout_lb) )
      {
         SCIP_CALL( SCIPinferVarLbCons(scip, pout_var, pin_lb, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );

         assert( ! infeasible );
         if ( tightened )
         {
            SCIPdebugMessage("New lower bound on %s: <%f> \n", SCIPvarGetName(pout_var), pin_lb);
            *result = SCIP_REDUCEDDOM;
         }
      }
      else if ( SCIPisGT(scip, pout_lb, pin_lb) )
      {
         SCIP_CALL( SCIPinferVarLbCons(scip, pin_var, pout_lb, cons, BPROP_LOWER, FALSE, &infeasible, &tightened) );

         assert( ! infeasible );
         if ( tightened )
         {
            SCIPdebugMessage("New lower bound on %s: <%f> \n", SCIPvarGetName(pin_var), pout_lb);
            *result = SCIP_REDUCEDDOM;
         }
      }

      if ( SCIPisGT(scip, pin_ub, pout_ub) )
      {
         SCIP_CALL( SCIPinferVarUbCons(scip, pin_var, pout_ub, cons, BPROP_UPPER, FALSE, &infeasible, &tightened) );

         assert( ! infeasible );
         if ( tightened )
         {
            SCIPdebugMessage("New upper bound on %s: <%f> \n", SCIPvarGetName(pin_var), pout_ub);
            *result = SCIP_REDUCEDDOM;
         }
      }
      else if ( SCIPisGT(scip, pout_ub, pin_ub) )
      {
         SCIP_CALL( SCIPinferVarUbCons(scip, pout_var, pin_ub, cons, FPROP_UPPER, FALSE, &infeasible, &tightened) );

         assert( ! infeasible );
         if ( tightened )
         {
            SCIPdebugMessage("New upper bound on %s: <%f> \n", SCIPvarGetName(pout_var), pin_ub);
            *result = SCIP_REDUCEDDOM;
         }
      }

      if ( *result == SCIP_INFEASIBLE )
      {
         /* add branching candidates */
         if ( SCIPisGT(scip, pout_ub - pout_lb, pin_ub - pin_lb) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, pout_var, consdata->violation, pout_lb / 2 + pout_ub / 2) );
         }
         else
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, pin_var, consdata->violation, pin_lb / 2 + pin_ub / 2) );
         }
      }

      /* Eventually it would make sense to return SCIP_SOLVELP instead of SCIP_INFEASIBLE,
       * since at least one of the inequalities pin <= pout and pout <= pin should be part
       * of the LP, due to the flow direction variables.
       * If both binary variables are fixed to 0 then pin=pout should be there. If only one
       * is fixed to 0, we could add the remaining inequality instead.
       */
      return SCIP_OKAY;
   }

   /* pseudo solutions are set to the variable bounds best for the objective function; check which bound they are set to: */
   q_lower    = SCIPisEQ(scip, ABS(q_sol), q_lb);
   pin_lower  = SCIPisEQ(scip, pin_sol, pin_lb);
   pout_lower = SCIPisEQ(scip, pout_sol, pout_lb);

   if ( consdata->type_violation == PIN_TOO_LOW )
   {
      if ( pout_lower && pin_lower && q_lower )
      {
         /* reduce domain of pin by bound propagation */
         new_p_bound = computeModifiedEulerScheme(scip, conshdlrdata, consdata, pout_lb, q_lb, -1);
         if ( new_p_bound > -0.5 )
         {
            if ( SCIPisGT(scip, new_p_bound, pin_lb) )
            {
               if ( dir == -1 )
                  SCIP_CALL( SCIPinferVarLbCons(scip, pin_var, new_p_bound, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );
               else
                  SCIP_CALL( SCIPinferVarLbCons(scip, pin_var, new_p_bound, cons, BPROP_LOWER, FALSE, &infeasible, &tightened) );

               assert( ! infeasible );
               if ( tightened )
               {
                  SCIPdebugMessage("New lower bound on %s: <%f> \n", SCIPvarGetName(pin_var), new_p_bound);
                  *result = SCIP_REDUCEDDOM;
               }
            }
         }
      }
      else if ( pout_lower && !pin_lower && q_lower )
      {
         /* there cannot be a feasible solution */
         *result = SCIP_CUTOFF;
      }
      else if ( pout_lower && !pin_lower && !q_lower )
      {
         /* upper flow bound is too big, try to compute a new one, otherwise branching on flow is preferred */
         SCIP_Real new_flow;

         tightened = FALSE;

         if ( bisectedFlowBound(scip, conshdlrdata, cons, &new_flow, dir) )
         {
            if ( dir == 1 )
            {
               SCIP_CALL( SCIPinferVarUbCons(scip, q_var, new_flow, cons, BPROP_LOWER, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );
            }
            else if ( dir == -1)
            {
               SCIP_CALL( SCIPinferVarLbCons(scip, q_var, new_flow, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );
               assert( ! infeasible );
            }
         }

         if ( tightened )
         {
            *result = SCIP_REDUCEDDOM;
         }
         else
         {
            preferFlow = TRUE;
         }
      }
      else if ( !pout_lower && !pin_lower && q_lower )
      {
         /* pout_ub is too big, try to reduce the domain by forward propagation */
         new_p_bound = computeModifiedEulerScheme(scip, conshdlrdata, consdata, pin_ub, q_lb, 1);
         if ( new_p_bound > -0.5 )
         {
            if ( SCIPisLT(scip, new_p_bound, pout_ub) )
            {
               if ( dir == 1 )
               {
                  SCIP_CALL( SCIPinferVarUbCons(scip, pout_var, new_p_bound, cons, FPROP_UPPER, FALSE, &infeasible, &tightened) );
               }
               else
               {
                  SCIP_CALL( SCIPinferVarUbCons(scip, pout_var, new_p_bound, cons, BPROP_UPPER, FALSE, &infeasible, &tightened) );
               }

               assert( ! infeasible );
               if ( tightened )
               {
                  SCIPdebugMessage("New upper bound on %s: <%f> \n", SCIPvarGetName(pout_var), new_p_bound);
                  *result = SCIP_REDUCEDDOM;
               }
            }
         }
      }
   }
   else if ( consdata->type_violation == PIN_TOO_BIG )
   {
      if ( pout_lower && pin_lower && !q_lower )
      {
         /* pout_lb is too small, try to reduce the domain by forward propagation */
         new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata, consdata, pin_lb, q_ub, 1);
         if ( new_p_bound > -0.5 )
         {
            if ( SCIPisGT(scip, new_p_bound, pout_lb) )
            {
               if ( dir == 1 )
               {
                  SCIP_CALL( SCIPinferVarLbCons(scip, pout_var, new_p_bound, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );
               }
               else
               {
                  SCIP_CALL( SCIPinferVarLbCons(scip, pout_var, new_p_bound, cons, BPROP_LOWER, FALSE, &infeasible, &tightened) );
               }

               assert( ! infeasible );
               if ( tightened )
               {
                  SCIPdebugMessage("New lower bound on %s: <%f> \n",SCIPvarGetName(pout_var), new_p_bound);
                  *result = SCIP_REDUCEDDOM;
               }
            }
         }
      }
      else if ( !pout_lower && pin_lower && q_lower )
      {
         /* q_lb is too small, eventually this could be changed by flowTightening based on bisection.
          * So far this is only implemented for the upper flow bound */
         SCIPinfoMessage(scip, NULL, "EnforcePS: Case q_lb too small, eventually flowTightening could be useful.\n");
         preferFlow = TRUE;
      }
      else if ( !pout_lower && pin_lower && !q_lower )
      {
         /* there cannot be a feasible solution */
         *result = SCIP_CUTOFF;
      }
      else if ( !pout_lower && !pin_lower && !q_lower )
      {
         /* pin_ub is too big, try to reduce the domain by bound propagation */
         new_p_bound = computeTrapezoidalScheme(scip, conshdlrdata, consdata, pout_ub, q_ub, -1);
         if ( new_p_bound > -0.5 )
         {
            if ( SCIPisLT(scip, new_p_bound, pin_ub) )
            {
               if ( dir == 1 )
               {
                  SCIP_CALL( SCIPinferVarUbCons(scip, pin_var, new_p_bound, cons, BPROP_UPPER, FALSE, &infeasible, &tightened) );
               }
               else
               {
                  SCIP_CALL( SCIPinferVarUbCons(scip, pin_var, new_p_bound, cons, FPROP_LOWER, FALSE, &infeasible, &tightened) );
               }

               assert( ! infeasible );
               if ( tightened )
               {
                  SCIPdebugMessage("New upper bound on %s: <%f> \n",SCIPvarGetName(pin_var), new_p_bound);
                  *result = SCIP_REDUCEDDOM;
               }
            }
         }
      }
   }
   else if ( consdata->type_violation == MACHBOUND )
   {
      /* add branching candidates */
   }

   /* If the infeasibility is not resolved yet, we add extern branching candidates.
    * Thereby we prefer p_out and q, because we can get good bounds on p_in by bound propagation */
   if ( *result == SCIP_INFEASIBLE )
   {
      if ( !SCIPisEQ(scip, q_ub, q_lb) )
      {
         SCIP_CALL( SCIPaddExternBranchCand(scip, q_var, consdata->violation, q_ub / 2 + q_lb / 2) );
      }

      if ( !preferFlow || (SCIPgetNExternBranchCands(scip) == 0) )
      {
         if ( !SCIPisEQ(scip, pout_ub, pout_lb) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, pout_var, consdata->violation, pout_lb / 2 + pout_ub / 2) );
         }
      }

      if ( SCIPgetNExternBranchCands(scip) == 0 )
      {
         if ( !SCIPisEQ(scip, pin_ub, pin_lb) )
         {
            SCIP_CALL( SCIPaddExternBranchCand(scip, pin_var, consdata->violation, pin_lb / 2 + pin_ub / 2) );
         }
      }
   }

   return SCIP_OKAY;
}

/** optimization based bound tightening */
static
SCIP_RETCODE obbt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of the pipe_ode constraint handler */
   SCIP_STATUS*          status,             /**< status pointer */
   int*                  n_bdchgs            /**< pointer to store number of bound changes */
   )
{
   char           probname[SCIP_MAXSTRLEN];
   SCIP_HASHMAP*  varmapfw;
   SCIP*          subscip;
   SCIP_PROBDATA* origdata;
   SCIP_VAR**     subvars;
   SCIP_VAR*      flowvar;
   SCIP_VAR*      subflowvar;
   SCIP_CONS**    conss;
   SCIP_Bool      success;
   SCIP_Real      dualbound;                 /* dual solution value of the subproblem */
   /* SCIP_Real      primalbound; */               /* primal solution value of the subproblem */
   SCIP_Real      ub, lb;
   SCIP_Real      lb_improvement = 0.0;
   SCIP_Real      ub_improvement = 0.0;
   SCIP_Real      timelimit = 3600.0;
   SCIP_RETCODE   retcode;
   int            nconss;
   int            nsubvars;
   int            on_vars = 0;
   int            i;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( n_bdchgs != NULL );
   assert( *n_bdchgs == 0 );

   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "obbt#subscip#");

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   if ( ! success )
      return SCIP_OKAY;

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPgetNVars(scip)) );

   /* create a problem copy as sub SCIP */
#if SCIP_VERSION >= 700
   SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, probname, TRUE, FALSE, TRUE, FALSE, &success ) );
#else
   SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, probname, TRUE, FALSE, FALSE, &success ) );
#endif
   if ( ! success )
   {
      SCIPerrorMessage("Was not able to copy scip during Obbt.\n");
      return SCIP_ERROR;
   }

   /* get variables of the subscip */
   subvars = SCIPgetVars(subscip);
   nsubvars = SCIPgetNVars(subscip);

   /* change objective value to 0.0 */
   for (i = 0; i < nsubvars; i++)
   {
      SCIP_CALL( SCIPchgVarObj(subscip, subvars[i], 0.0) );
   }

#ifdef SCIP_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* copy time and memory limits */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );

   SCIP_CALL( SCIPgetRealParam(subscip, "limits/time", &timelimit) );
   if ( timelimit > 600.0 )
   {
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", 600.0) );
   }

   /* problem get data */
   origdata = SCIPgetProbData(scip);

   /* delete nonlinear constraints, such that problem can be solved fast */
   conss = SCIPgetConss(subscip);
   nconss = SCIPgetNConss(subscip);

   for (i = nconss-1; i >= 0; --i)
   {
#if SCIP_VERSION >= 900
      if ( strcmp("nonlinear", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 )
      {
         SCIP_CALL( SCIPdelCons(subscip, conss[i]) );
      }
#else
      if ( strcmp("abspower", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 )
      {
         SCIP_CALL( SCIPdelCons(subscip, conss[i]) );
      }
      else if ( strcmp("quadratic", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 )
      {
         SCIP_CALL( SCIPdelCons(subscip, conss[i]) );
      }
#endif
      else if ( strcmp("pipe_ode", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 )
      {
         SCIP_CALL( SCIPdelCons(subscip, conss[i]) );
      }
      /* else if ( strcmp("linear", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 ) */
      /* { */
      /*    /\* do nothing *\/ */
      /* } */
      /* else if ( strcmp("varbound", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 ) */
      /* { */
      /*    /\* do also nothing *\/ */
      /* } */
      /* else if ( strcmp("and", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 ) */
      /* { */
      /*    /\* do nothing *\/ */
      /* } */
      /* else if ( strcmp("bounddisjunction", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 ) */
      /* { */
      /*    /\* do also nothing *\/ */
      /* } */
      /* else */
      /* { */
      /*    printf("%s\n", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i]))); */
      /* } */
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) starting OBBT ...\n", SCIPgetSolvingTime(scip));

   /* loop over all arcs and try to minimize and maximize the flow variables */
   for (i = origdata->network->numarcs - 1; i >= 0; --i)
   {
      flowvar = SCIPvarGetTransVar(origdata->flowvars[i]);
      subflowvar = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, flowvar);

      /* possibly some variables have been deleted */
      if ( subflowvar == NULL )
         continue;

      /* can not change bounds of multi-aggregated variables */
      assert( SCIPvarGetStatus(flowvar) != SCIP_VARSTATUS_ORIGINAL );
      assert( SCIPvarGetStatus(flowvar) != SCIP_VARSTATUS_NEGATED );
      if ( SCIPvarGetStatus(flowvar) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(flowvar) == SCIP_VARSTATUS_AGGREGATED || SCIPvarGetStatus(flowvar) == SCIP_VARSTATUS_FIXED )
         continue;

      /* only try to compute new bounds if the current ones are not close */
      lb = SCIPcomputeVarLbGlobal(scip, flowvar);
      ub = SCIPcomputeVarUbGlobal(scip, flowvar);

      /* do not perform obbt, if bounds are close; obbt_mingap can be changed with the settings */
      if ( ub - lb < conshdlrdata->obbt_mingap )
         continue;

      /* count variables obbt is performed on */
      ++on_vars;

      /* First try to minimize the flow varible.
       * One could also try to maximize first, there is no reason for this order. */
      SCIP_CALL( SCIPchgVarObj(subscip, subflowvar, 1.0) );

      /* solve subscip */
      retcode = SCIPsolve(subscip);

      if ( retcode != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving subproblem in OBBT; sub-SCIP terminated with code <%d>.\n", retcode);
         SCIPABORT();

         break;
      }

      /* get solving status */
      *status = SCIPgetStatus(subscip);

      if ( *status == SCIP_STATUS_USERINTERRUPT )
      {
         break;
      }
      else if ( *status == SCIP_STATUS_INFEASIBLE )
      {
         break;
      }
      else if ( *status == SCIP_STATUS_TIMELIMIT )
      {
         break;
      }

      /* get dual bound of solution process -> provides variable bound in original problem */
      dualbound = SCIPgetDualbound(subscip);
      /* primalbound = SCIPgetPrimalbound(subscip); */

      /* free transformed problem such that we can optimize again */
      SCIP_CALL( SCIPfreeTransform(subscip) );

      /* have to change bounds of subflowvar after freeing the transformed problem */
      if ( *status == SCIP_STATUS_OPTIMAL )
      {
         if ( SCIPisFeasGT(scip, dualbound, lb) && SCIPisFeasLT(scip, dualbound, ub) )
         {
            /* SCIPinfoMessage(subscip, NULL,"Min: %-27s   Dual: %11.5f   Primal: %11.5f   LB: %11.5f   UB: %11.5f\n", SCIPvarGetName(subflowvar), dualbound, primalbound, lb, ub); */

            SCIP_CALL( SCIPchgVarLbGlobal(scip, flowvar, dualbound) );
            SCIP_CALL( SCIPchgVarLbGlobal(subscip, subflowvar, dualbound) );
            lb_improvement += dualbound - lb;
            lb = dualbound;
            ++(*n_bdchgs);
         }
         else if ( SCIPisFeasEQ(scip, ub, dualbound) )
         {
            SCIP_CALL( SCIPchgVarLbGlobal(scip, flowvar, ub) );
            SCIP_CALL( SCIPchgVarLbGlobal(subscip, subflowvar, ub) );
            lb_improvement += ub - lb;
            lb = ub;
            ++(*n_bdchgs);
         }
      }

      /* if new bounds are already tight, trying to maximize the flow is pointless */
      if ( SCIPisFeasEQ(scip, ub, lb) )
      {
         SCIP_CALL( SCIPchgVarObj(subscip, subflowvar, 0.0) );
         continue;
      }

      /* try to maximize flow */
      SCIP_CALL( SCIPchgVarObj(subscip, subflowvar, -1.0) );

      /* resolve subscip */
      retcode = SCIPsolve(subscip);

      if ( retcode != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving subproblem in OBBT; sub-SCIP terminated with code <%d>.\n", retcode);
         SCIPABORT();

         break;
      }

      /* get solving status */
      *status = SCIPgetStatus(subscip);

      if ( *status == SCIP_STATUS_USERINTERRUPT )
      {
         break;
      }
      else if ( *status == SCIP_STATUS_INFEASIBLE )
      {
         break;
      }
      else if ( *status == SCIP_STATUS_TIMELIMIT )
      {
         break;
      }

      /* get dual bound of solution process -> provides variable bound in original problem */
      dualbound = - SCIPgetDualbound(subscip);
      /* primalbound = - SCIPgetPrimalbound(subscip); */

      /* free transformed problem such that we can optimize again */
      SCIP_CALL( SCIPfreeTransform(subscip) );

      /* have to change bounds of subflowvar after freeing the transformed problem */
      if ( *status == SCIP_STATUS_OPTIMAL )
      {
         if ( SCIPisFeasGT(scip, dualbound, lb) && SCIPisFeasLT(scip, dualbound, ub) )
         {
            /* SCIPinfoMessage(subscip, NULL,"Max: %-27s   Dual: %11.5f   Primal: %11.5f   LB: %11.5f   UB: %11.5f\n",
               SCIPvarGetName(flowvar), dualbound, primalbound, lb, ub); */

            SCIP_CALL( SCIPchgVarUbGlobal(scip, flowvar, dualbound) );
            SCIP_CALL( SCIPchgVarUbGlobal(subscip, subflowvar, dualbound) );
            ub_improvement += ub - dualbound;
            ++(*n_bdchgs);
         }
         else if ( SCIPisFeasEQ(scip, lb, dualbound) )
         {
            SCIP_CALL( SCIPchgVarUbGlobal(scip, flowvar, lb) );
            SCIP_CALL( SCIPchgVarUbGlobal(subscip, subflowvar, lb) );
            ub_improvement += ub - lb;
            ++(*n_bdchgs);
         }
      }

      SCIP_CALL( SCIPchgVarObj(subscip, subflowvar, 0.0) );

      /* go to next variable */
   }

   /* some output */
   if ( *status != SCIP_STATUS_INFEASIBLE && *status != SCIP_STATUS_USERINTERRUPT )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) OBBT on %d variables produced %d bound changes in run %d.\n",
         SCIPgetSolvingTime(scip), on_vars, *n_bdchgs, SCIPgetNRuns(scip));
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "           Improved lower bounds by: %10.2f\n", lb_improvement);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "           Improved upper bounds by: %10.2f\n", ub_improvement);
   }

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* finally free subscipdata */
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}

/** function that tries to remove an ODE constraint if flow is (almost) fixed */
static
SCIP_RETCODE removeODECons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< pipe ode constraint handler */
   SCIP_CONS*            cons,               /**< constraint to be deleted/deactivated */
   SCIP_RESULT*          result,             /**< pointer to store result */
   int*                  nfixedvars,         /**< pointer to count the total number of fixed variables */
   int*                  naggrvars,          /**< pointer to count the total number of aggregated variables */
   int*                  ndelcons            /**< pointer to count the total number of deleted constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   SCIP_VARSTATUS     varstatus;
   SCIP_Bool          infeasible;
   SCIP_Bool          fixed;
   SCIP_Bool          aggregated;
   SCIP_Bool          redundant;
   SCIP_Bool          delete;
   SCIP_Real          qlb;
   SCIP_Real          qub;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( result != NULL );
   assert( nfixedvars != NULL );
   assert( naggrvars != NULL );
   assert( ndelcons != NULL );

   /* initialize data and variables */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   *result    = SCIP_DIDNOTFIND;
   infeasible = FALSE;
   fixed      = FALSE;
   delete     = FALSE;
   aggregated = FALSE;
   redundant  = FALSE;

   qlb = SCIPcomputeVarLbGlobal(scip, consdata->q_var);
   qub = SCIPcomputeVarUbGlobal(scip, consdata->q_var);

   varstatus = SCIPvarGetStatus(consdata->q_var);
   assert( varstatus != SCIP_VARSTATUS_ORIGINAL );

   /* try to fix flow to zero if not already */
   if ( varstatus != SCIP_VARSTATUS_FIXED )
   {
      SCIP_CALL( SCIPfixVar(scip, consdata->q_var, 0.0, &infeasible, &fixed) );

      if ( infeasible )
      {
         /* Could not fix flow variable to zero, i.e., either UB is FeasNegative or LB is FeasPositive.
          * Nevertheless, our custom feasibility tolerance for flows, accepts the flow as zero.
          * Therefore, we fix the variable to its "smaller" bound, i.e., the one closer to zero. */
         if ( SCIPisFeasPositive(scip, qlb) )
         {
            SCIP_CALL( SCIPfixVar(scip, consdata->q_var, qlb, &infeasible, &fixed) );
         }
         else if ( SCIPisFeasNegative(scip, qub) )
         {
            SCIP_CALL( SCIPfixVar(scip, consdata->q_var, qub, &infeasible, &fixed) );
         }
         else
         {
            SCIPerrorMessage("Error in removeODECons on pipe <%s>.\n", SCIPconsGetName(cons));
            SCIPerrorMessage("Could not fix flow variable <%s>:\n", SCIPvarGetName(consdata->q_var));
            SCIPerrorMessage("Lower Bound: %f\n", qlb);
            SCIPerrorMessage("Upper Bound: %f\n", qub);

            return SCIP_ERROR;
         }
      }

      if ( fixed )
      {
         SCIPdebugMessage("Fixed flow variable <%s> to %f.\n", SCIPvarGetName(consdata->q_var), SCIPcomputeVarLbGlobal(scip, consdata->q_var));
         ++(*nfixedvars);
      }
      else if ( varstatus == SCIP_VARSTATUS_AGGREGATED )
      {
         SCIP_VAR* aggrflow;
         SCIP_Real scalar;
         SCIP_Real constant;

         aggrflow = consdata->q_var;
         scalar = 1.0;
         constant = 0.0;
         SCIP_CALL( SCIPgetProbvarSum(scip, &aggrflow, &scalar, &constant) );

         if ( SCIPvarGetStatus(aggrflow) != SCIP_VARSTATUS_FIXED )
         {
            SCIPerrorMessage("Aggregated flow variable <%s> should be fixed at this point of removeODECons.\n", SCIPvarGetName(consdata->q_var));

            return SCIP_ERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Flow variable <%s> should be fixed at this point of removeODECons.\n", SCIPvarGetName(consdata->q_var));

         return SCIP_ERROR;
      }
   }
   else
   {
      /* check if fixedval is feasZero or only zero to custom feasibility tolerance? */

      SCIPdebugMessage("Flow variable <%s> already fixed to %f.\n", SCIPvarGetName(consdata->q_var), qlb);
   }

   /* updata flow bounds */
   if ( fixed )
   {
      qlb = SCIPcomputeVarLbGlobal(scip, consdata->q_var);
      qub = SCIPcomputeVarUbGlobal(scip, consdata->q_var);
   }

   /* fix binary flow variables */
   varstatus = SCIPvarGetStatus(consdata->posFlowBinvar);

   if ( varstatus == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("Cannot fix posFlowBinvar <%s> in removeODECons because it is multi-aggregated.\n", SCIPvarGetName(consdata->posFlowBinvar));

      return SCIP_ERROR;
   }
   else if ( varstatus != SCIP_VARSTATUS_FIXED )
   {
      SCIP_CALL( SCIPfixVar(scip, consdata->posFlowBinvar, 0.0, &infeasible, &fixed) );

      if ( infeasible )
      {
         if ( SCIPisFeasPositive(scip, qlb) )
         {
            SCIPdebugMessage("Cannot fix posFlowBinvar <%s> to 0 because flow is only zero for custom tolerance.\n", SCIPvarGetName(consdata->posFlowBinvar));
         }
         else if ( ! SCIPisFeasPositive(scip, qub) )
         {
            SCIPerrorMessage("Cannot fix posFlowBinvar <%s> to 0 although flow is not feas-positive.\n", SCIPvarGetName(consdata->posFlowBinvar));

            return SCIP_ERROR;
         }
      }

      if ( fixed )
      {
         SCIPdebugMessage("Fixed posFlowBinvar <%s> to 0.\n", SCIPvarGetName(consdata->posFlowBinvar));
         ++(*nfixedvars);
      }
   }
   else
   {
      SCIPdebugMessage("posFlowBinvar <%s> already fixed to %f.\n", SCIPvarGetName(consdata->posFlowBinvar), SCIPcomputeVarLbGlobal(scip, consdata->posFlowBinvar));
   }

   varstatus = SCIPvarGetStatus(consdata->negFlowBinvar);

   if ( varstatus == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIPerrorMessage("Cannot fix negFlowBinvar <%s> in removeODECons because it is multi-aggregated.\n", SCIPvarGetName(consdata->negFlowBinvar));

      return SCIP_ERROR;
   }
   else if ( varstatus != SCIP_VARSTATUS_FIXED )
   {
      SCIP_CALL( SCIPfixVar(scip, consdata->negFlowBinvar, 0.0, &infeasible, &fixed) );

      if ( infeasible )
      {
         if ( SCIPisFeasNegative(scip, qub) )
         {
            SCIPdebugMessage("Cannot fix negFlowBinvar <%s> to 0 because flow is only zero for custom tolerance.\n", SCIPvarGetName(consdata->negFlowBinvar));
         }
         else if ( ! SCIPisFeasNegative(scip, qlb) )
         {
            SCIPerrorMessage("Cannot fix negFlowBinvar <%s> to 0 although flow is not feas-negative.\n", SCIPvarGetName(consdata->negFlowBinvar));

            return SCIP_ERROR;
         }
      }

      if ( fixed )
      {
         SCIPdebugMessage("Fixed negFlowBinvar <%s> to 0.\n", SCIPvarGetName(consdata->negFlowBinvar));
         ++(*nfixedvars);
      }
   }
   else
   {
      SCIPdebugMessage("negFlowBinvar <%s> already fixed to %f.\n", SCIPvarGetName(consdata->negFlowBinvar), SCIPcomputeVarLbGlobal(scip, consdata->negFlowBinvar));
   }

   assert( SCIPvarGetStatus(consdata->p_in_var)  != SCIP_VARSTATUS_MULTAGGR );
   assert( SCIPvarGetStatus(consdata->p_out_var) != SCIP_VARSTATUS_MULTAGGR );

   SCIP_CALL( SCIPaggregateVars(scip, consdata->p_in_var, consdata->p_out_var, 1.0, -1.0, 0.0, &infeasible, &redundant, &aggregated) );

   /* if aggregation is infeasible, i.e., pin_lb > pout_ub or pin_ub > pout_lb, check if there is a feasible solution for custom comparison. */
   /* if node is infeasible, bound propagation probably detects this...  */
   if ( infeasible )
   {
      SCIP_Real max_lb;
      SCIP_Real min_ub;

      max_lb = MAX( SCIPcomputeVarLbGlobal(scip, consdata->p_in_var), SCIPcomputeVarLbGlobal(scip, consdata->p_out_var) ); /*lint !e666 */
      min_ub = MAX( SCIPcomputeVarUbGlobal(scip, consdata->p_in_var), SCIPcomputeVarUbGlobal(scip, consdata->p_out_var) ); /*lint !e666 */

      if ( pressureIsGT(max_lb, min_ub, conshdlrdata->pressure_feastol) )
      {
         *result = SCIP_CUTOFF;

         return SCIP_OKAY;
      }
   }

   /* if aggregation was successfull or condition is redundant we can delete the ODE constraint */
   if ( aggregated || redundant )
      delete = TRUE;

   if ( aggregated )
      ++(*naggrvars);

   /* delete ODE constraint if possible */
   if ( delete )
   {
      SCIPdebugMessage("Delete pipe_ode constraint <%s>.\n", SCIPconsGetName(cons));
      SCIP_CALL( SCIPdelCons(scip, cons) );

      *result = SCIP_SUCCESS;
      ++(*ndelcons);
   }

   return SCIP_OKAY;
}


/**
 *  Callback methods of subnlp heuristik
 */

/** user expression data */
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
struct SCIP_ExprData
#else
struct SCIP_UserExprData
#endif
{
   SCIP*                 scip;
   SCIP_CONSHDLRDATA*    conshdlrdata;       /**< data of constraint handler */
   SCIP_CONSDATA*        consdata;           /**< problem data structure */
   SCIP_EXPRCURV         curvature;          /**< curvature information */
};

#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )

/** check the curvature of the user expression */
static
SCIP_DECL_EXPRCURVATURE(SCIPexprCurvODE)
{
   SCIP_EXPRDATA* data;
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( expr != NULL );
   assert( success != NULL );

   *success = FALSE;

   data = SCIPexprGetData(expr);

   assert( data != NULL );
   assert( data->consdata != NULL );
   assert( data->scip == scip );

   consdata = data->consdata;
   conshdlrdata = data->conshdlrdata;

   /* if the flow is 0 */
   if ( SCIPisZero(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var)) && SCIPisZero(scip, SCIPcomputeVarUbLocal(data->scip, consdata->q_var)) )
   {
      data->curvature = SCIP_EXPRCURV_LINEAR;
      if ( exprcurvature == SCIP_EXPRCURV_LINEAR )
      {
         *success = TRUE;
         childcurv[0] = SCIP_EXPRCURV_LINEAR;
         childcurv[1] = SCIP_EXPRCURV_LINEAR;
      }
      else if ( exprcurvature == SCIP_EXPRCURV_CONVEX )
      {
         *success = TRUE;
         childcurv[0] = SCIP_EXPRCURV_LINEAR;
         childcurv[1] = SCIP_EXPRCURV_LINEAR;
      }
   }
   else if ( conshdlrdata->nlp_convex )
   {
      /* if the flow direction is fixed */
      if ( SCIPisGE(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var), 0.0) || SCIPisLE(scip, SCIPcomputeVarUbLocal(scip, consdata->q_var), 0.0) )
      {
         data->curvature = SCIP_EXPRCURV_CONVEX;
         if ( exprcurvature == SCIP_EXPRCURV_CONVEX )
         {
            *success = TRUE;
            childcurv[0] = SCIP_EXPRCURV_LINEAR;
            childcurv[1] = SCIP_EXPRCURV_LINEAR;
         }
      }
   }

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback for the expression "g(p_in, p_out, q) = p^l(p_out, q) - p_in <= 0.0"
 *
 *  The variables are p_in, p_out, q in this order. q^l is the lower bound obtained, e.g., my the modified Euler
 *  scheme. This function is convex, so the inequality is convex as well.
 *
 *  We have the following cases:
 *  - q = 0: In this case g(p_in, p_out, q) = p_out - p_in.
 *  - q > 0: Compute p^l(p_out,q) by the modified Euler scheme. Then g = p^l - p_in.
 *  - q < 0: Compute p^l(p_in,-q) by the modified Euler scheme with negative flow. Then g = p^l - p_out.
 */
static
SCIP_DECL_EXPREVAL(SCIPexprEvalODE)
{
   SCIP_EXPRDATA* data;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR** children;
   SCIP_Real q;
   SCIP_Real pin;
   SCIP_Real pout;
   SCIP_Real plb;

   assert( scip != NULL );
   assert( expr != NULL );
   assert( val != NULL );

   data = SCIPexprGetData(expr);
   assert( data != NULL );
   assert( scip == data->scip );
   assert( SCIPexprGetNChildren(expr) == 3 );

   consdata = data->consdata;
   assert( consdata != NULL );

   children = SCIPexprGetChildren(expr);
   assert( children != NULL );

   pin = SCIPexprGetEvalValue(children[0]);
   pout = SCIPexprGetEvalValue(children[1]);
   q = SCIPexprGetEvalValue(children[2]);

   /* determine value depending on flow direction */
   if ( SCIPisFeasZero(scip, q) )
   {
      /* flow is 0 */
      *val = pout - pin;
   }
   else if ( q > 0.0 )
   {
      /* compute lower bound in reverse direction, so pout is the input */
      plb = computeModifiedEulerScheme(scip, data->conshdlrdata, consdata, pout, q, -1);

      if ( SCIPisNegative(scip, plb) )
         *val = SCIP_INVALID;
      else
         *val = plb - pin;
   }
   else
   {
      assert( q < 0.0 );

      /* if we want a convex relaxation, we use the modified Euler scheme */
      if ( data->conshdlrdata->nlp_convex )
      {
         /* compute lower bound in reverse direction with negative flow */
         plb = computeModifiedEulerScheme(scip, data->conshdlrdata, consdata, pin, -q, -1);
      }
      else
      {
         /* otherwise we use the more exact RK4 */
         plb = computeRK4(scip, data->conshdlrdata, consdata, pin, -q, -1);
      }

      if ( SCIPisNegative(scip, plb) )
         *val = SCIP_INVALID;
      else
         *val = plb - pout;
   }

   return SCIP_OKAY;
}

/** backward derivative evaluation callback
 *
 *  We have the following cases:
 *  - q = 0: The gradient is (-1, 1, 0).
 *  - q > 0: Compute p^l(p_out,q). The gradient is (-1, \partial p^l/\partial p_out, \partial p^l/\partial q).
 *  - q < 0: Compute p^l(p_in,-q). The gradient is (\partial p^l/\partial p_in, -1, - \partial p^l/\partial q).
 */
static
SCIP_DECL_EXPRBWDIFF(SCIPexprBWODE)
{
   SCIP_EXPRDATA* data;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR** children;
   SCIP_Real q;
   SCIP_Real pin;
   SCIP_Real pout;

   assert( scip != NULL );
   assert( expr != NULL );
   assert( val != NULL );
   assert( childidx == 0 || childidx == 1 || childidx == 2 );

   data = SCIPexprGetData(expr);
   assert( data != NULL );
   assert( scip == data->scip );

   consdata = data->consdata;
   assert( consdata != NULL );

   children = SCIPexprGetChildren(expr);
   assert( children != NULL );

   pin = SCIPexprGetEvalValue(children[0]);
   pout = SCIPexprGetEvalValue(children[1]);
   q = SCIPexprGetEvalValue(children[2]);

   /* determine value depending on flow direction */
   if ( SCIPisFeasZero(scip, q) )
   {
      /* flow is 0 */
      if ( childidx == 0 )
         *val = -1.0;
      else if ( childidx == 1 )
         *val = 1.0;
      else
      {
         assert( childidx == 2 );
         *val = 0.0;
      }
   }
   else if ( q > 0.0 )
   {
      if ( childidx == 0 )
         *val = -1.0;
      else
      {
         SCIP_Real dpout;
         SCIP_Real dq;

         if ( data->conshdlrdata->nlp_convex )
            (void) computeModifiedEulerSchemeWithGradients(scip, data->conshdlrdata, consdata, pout, q, &dpout, &dq);
         else
            (void) computeRK4WithGradients(scip, data->conshdlrdata, consdata, pout, q, &dpout, &dq);

         if ( childidx == 1 )
            *val = dpout;
         else
         {
            assert( childidx == 2 );
            *val = dq;
         }
      }
   }
   else
   {
      assert( q < 0.0 );

      if ( childidx == 1 )
         *val = -1.0;
      else
      {
         SCIP_Real dpin;
         SCIP_Real dq;

         assert( childidx == 0 || childidx == 2 );
         if ( data->conshdlrdata->nlp_convex )
            (void) computeModifiedEulerSchemeWithGradients(scip, data->conshdlrdata, consdata, pin, REALABS(q), &dpin, &dq);
         else
            (void) computeRK4WithGradients(scip, data->conshdlrdata, consdata, pin, REALABS(q), &dpin, &dq);

         if ( childidx == 0 )
            *val = dpin;
         else
         {
            assert( childidx == 2 );
            *val = - dq; /* we need to reverse the direction */
         }
      }
   }

   return SCIP_OKAY;
}

/** forward derivative evaluation callback
 *
 *  Directional derivative callback for the constraint "p^l(p_out, q) - p_in <= 0.0".
 *
 *  The variables are p_in, p_out, q in this order.
 *
 *  We have the following cases:
 *  - q = 0: The gradient is (-1, 1, 0).
 *  - q > 0: Compute p^l(p_out,q). The gradient is (-1, \partial p^l/\partial p_out, \partial p^l/\partial q).
 *  - q < 0: Compute p^l(p_in,-q). The gradient is (\partial p^l/\partial p_in, -1, - \partial p^l/\partial q).
 */
static
SCIP_DECL_EXPRFWDIFF(SCIPexprFWdiffODE)
{
   SCIP_EXPRDATA* data;
   SCIP_CONSDATA* consdata;
   SCIP_EXPR** children;
   SCIP_Real q;
   SCIP_Real pin;
   SCIP_Real pout;

   assert( scip != NULL );
   assert( expr != NULL );
   assert( dot != NULL );

   assert( SCIPexprGetNChildren(expr) == 3 );

   data = SCIPexprGetData(expr);
   assert( data != NULL );
   assert( scip == data->scip );

   consdata = data->consdata;
   assert( consdata != NULL );

   children = SCIPexprGetChildren(expr);
   assert( children != NULL );

   pin = SCIPexprGetEvalValue(children[0]);
   pout = SCIPexprGetEvalValue(children[1]);
   q = SCIPexprGetEvalValue(children[2]);

   assert( SCIPexprGetDot(children[0]) != SCIP_INVALID );
   assert( SCIPexprGetDot(children[1]) != SCIP_INVALID );
   assert( SCIPexprGetDot(children[2]) != SCIP_INVALID );

   /* determine value depending on flow direction */
   *dot = 0.0;
   if ( SCIPisFeasZero(scip, q) )
   {
      /* flow is 0 */
      *dot += -1.0 * SCIPexprGetDot(children[0]);
      *dot += 1.0 * SCIPexprGetDot(children[1]);
   }
   else if ( q > 0.0 )
   {
      SCIP_Real dpout;
      SCIP_Real dq;

      if ( data->conshdlrdata->nlp_convex )
         (void) computeModifiedEulerSchemeWithGradients(scip, data->conshdlrdata, consdata, pout, q, &dpout, &dq);
      else
         (void) computeRK4WithGradients(scip, data->conshdlrdata, consdata, pout, q, &dpout, &dq);

      *dot += -1.0 * SCIPexprGetDot(children[0]);
      *dot += dpout * SCIPexprGetDot(children[1]);
      *dot += dq * SCIPexprGetDot(children[2]);
   }
   else
   {
      SCIP_Real dpin;
      SCIP_Real dq;

      assert( q < 0.0 );

      if ( data->conshdlrdata->nlp_convex )
         (void) computeModifiedEulerSchemeWithGradients(scip, data->conshdlrdata, consdata, pin, REALABS(q), &dpin, &dq);
      else
         (void) computeRK4WithGradients(scip, data->conshdlrdata, consdata, pin, REALABS(q), &dpin, &dq);

      *dot += dpin * SCIPexprGetDot(children[0]);
      *dot += -1.0 * SCIPexprGetDot(children[1]);
      *dot += - dq * SCIPexprGetDot(children[2]); /* reverse sign */
   }

   return SCIP_OKAY;
}


/** copy user data */
static
SCIP_DECL_EXPRCOPYDATA(SCIPexprCopyODE)
{
   SCIP_EXPRDATA* sourceexprdata;

   sourceexprdata = SCIPexprGetData(sourceexpr);
   assert( sourceexprdata != NULL );

   SCIP_CALL( SCIPduplicateBlockMemory(targetscip, targetexprdata, sourceexprdata) );

   return SCIP_OKAY;
}

/** free user expression data callback */
static
SCIP_DECL_EXPRFREEDATA(SCIPexprFreeODE)
{
   SCIP_EXPRDATA* exprdata;

   assert( expr != NULL );
   exprdata = SCIPexprGetData(expr);
   assert( exprdata != NULL );

   SCIPfreeBlockMemory(scip, &exprdata);

   return SCIP_OKAY;
}


/* ----------------------------------------------------*/
#else
/* ----------------------------------------------------*/


/** check the curvature of the user expression */
static
SCIP_DECL_USEREXPRCURV(SCIPuserexprCurvODE)
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP* scip;

   assert( data != NULL );
   assert( data->consdata != NULL );

   scip = data->scip;
   consdata = data->consdata;
   conshdlrdata = data->conshdlrdata;

   /* if the flow is 0 */
   /* if the flow is 0 */
   if ( SCIPisZero(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var)) && SCIPisZero(scip, SCIPcomputeVarUbLocal(scip, consdata->q_var)) )
   {
      data->curvature = SCIP_EXPRCURV_LINEAR;
   }
   else if ( conshdlrdata->nlp_convex )
   {
      /* if the flow direction is fixed */
      if ( SCIPisGE(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var), 0.0) || SCIPisLE(scip, SCIPcomputeVarUbLocal(scip, consdata->q_var), 0.0) )
      {
         data->curvature = SCIP_EXPRCURV_CONVEX;
      }
   }

   *result = data->curvature;

   return SCIP_OKAY;
}


/** function eval callback for the constraint "p^l(p_out, q) - p_in <= 0.0"
 *
 *  input:  SCIP_USEREXPRDATA* data,
 *          int nargs,                    << should be three
 *          SCIP_Real* argvals,           << values of the expression variables
 *          SCIP_Real* funcvalue,         << pointer to store the function evaluation, SCIP_INVALID if not defined
 *          SCIP_Real* gradient,          << if not null store the gradient here
 *          SCIP_Real* hessian            << should be NULL, since exprcapability is eval and gradient only
 *
 *  output: SCIP_RETCODE
 */
static
SCIP_DECL_USEREXPREVAL(SCIPuserexprEvalODE)
{
   SCIP_CONSDATA* consdata;
   SCIP* scip;
   SCIP_Real plb;

   assert( data != NULL );
   assert( nargs == 3 );
   assert( argvals != NULL );
   assert( funcvalue != NULL );
   assert( hessian == NULL );

   scip = data->scip;
   assert( scip != NULL );
   consdata = data->consdata;
   assert( consdata != NULL );

   /* determine value/gradients depending on flow direction */
   if ( SCIPisFeasZero(scip, argvals[2]) )
   {
      /* flow is 0 */
      *funcvalue = argvals[1] - argvals[0];
      if( gradient != NULL )
      {
         gradient[0] = -1.0;
         gradient[1] = 1.0;
         gradient[2] = 0.0;
      }
   }
   else if ( argvals[2] > 0.0 )
   {
      if ( gradient != NULL )
      {
         if ( data->conshdlrdata->nlp_convex )
            plb = computeModifiedEulerSchemeWithGradients(scip, data->conshdlrdata, consdata, argvals[1], argvals[2], &gradient[1], &gradient[2]);
         else
            plb = computeRK4WithGradients(scip, data->conshdlrdata, consdata, argvals[1], argvals[2], &gradient[1], &gradient[2]);

         gradient[0] = -1.0;
      }
      else
      {
         if ( data->conshdlrdata->nlp_convex )
            plb = computeModifiedEulerScheme(scip, data->conshdlrdata, consdata, argvals[1], argvals[2], -1);
         else
            plb = computeRK4(scip, data->conshdlrdata, consdata, argvals[1], argvals[2], -1);
      }

      if ( SCIPisNegative(scip, plb) )
         *funcvalue = SCIP_INVALID;
      else
         *funcvalue = plb - argvals[0];
   }
   else
   {
      assert( argvals[2] < 0.0 );
      if ( gradient != NULL )
      {
         if ( data->conshdlrdata->nlp_convex )
            plb = computeModifiedEulerSchemeWithGradients(scip, data->conshdlrdata, consdata, argvals[0], REALABS(argvals[2]), &gradient[0], &gradient[2]);
         else
            plb = computeRK4WithGradients(scip, data->conshdlrdata, consdata, argvals[0], REALABS(argvals[2]), &gradient[0], &gradient[2]);

         gradient[1] = -1.0;
         gradient[2] *= -1.0;  /* we need to reverse the direction */
      }
      else
      {
         if ( data->conshdlrdata->nlp_convex )
            plb = computeModifiedEulerScheme(scip, data->conshdlrdata, consdata, argvals[0], REALABS(argvals[2]), -1);
         else
            plb = computeRK4(scip, data->conshdlrdata, consdata, argvals[0], REALABS(argvals[2]), -1);
      }

      if ( SCIPisNegative(scip, plb) )
         *funcvalue = SCIP_INVALID;
      else
         *funcvalue = plb - argvals[1];
   }

   return SCIP_OKAY;
}

/** copy user data */
static
SCIP_DECL_USEREXPRCOPYDATA(SCIPuserexprCopyODE)
{
   assert( datasource != NULL );
   SCIP_ALLOC( BMSduplicateBlockMemory(blkmem, datatarget, datasource) );

   return SCIP_OKAY;
}

/** free user expression data callback */
static
SCIP_DECL_USEREXPRFREEDATA(SCIPuserexprFreeODE)
{
   assert( data != NULL );
   BMSfreeBlockMemory(blkmem, &data);
}
#endif

/**
 *  Callback methods of event handler
 */

/** exec the event handler */
static
SCIP_DECL_EVENTEXEC(eventExecPipeODE)
{  /*lint --e{715}*/
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   SCIP_VAR* eventvar;

   assert( eventhdlr != NULL );
   assert( eventdata != NULL );
   assert( strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0 );
   assert( event != NULL );

   cons = (SCIP_CONS*)eventdata;             /*lint !e740*/
   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);

   eventvar = SCIPeventGetVar(event);
   assert( eventvar != NULL );

   eventtype = SCIPeventGetType(event);
   switch ( eventtype )
   {
   case SCIP_EVENTTYPE_LBRELAXED:
      if ( eventvar == consdata->p_in_var )
      {
         if ( consdata->propvarind & PIPEODE_PIN_LB )
            consdata->propvarind ^= PIPEODE_PIN_LB;
      }
      else if ( eventvar == consdata->p_out_var )
      {
         if ( consdata->propvarind & PIPEODE_POUT_LB )
            consdata->propvarind ^= PIPEODE_POUT_LB;
      }
      else if ( eventvar == consdata->q_var )
      {
         if ( consdata->propvarind & PIPEODE_Q_LB )
            consdata->propvarind ^= PIPEODE_Q_LB;
      }
      else if ( eventvar == consdata->posFlowBinvar )
      {
         if ( consdata->propvarind & PIPEODE_POS_LB )
            consdata->propvarind ^= PIPEODE_POS_LB;
      }
      else if ( eventvar == consdata->negFlowBinvar )
      {
         if ( consdata->propvarind & PIPEODE_NEG_LB )
            consdata->propvarind ^= PIPEODE_NEG_LB;
      }
      break;

   case SCIP_EVENTTYPE_LBTIGHTENED:

      SCIPdebugMessage("<%s>: lower bound of p_in|p_out|q|posFlowBinvar|negFlowBinvar tightened changed.\n", SCIPconsGetName(cons));

      if ( eventvar == consdata->p_in_var )
         consdata->propvarind |= PIPEODE_PIN_LB;
      else if ( eventvar == consdata->p_out_var )
         consdata->propvarind |= PIPEODE_POUT_LB;
      else if ( eventvar == consdata->q_var )
         consdata->propvarind |= PIPEODE_Q_LB;
      else if ( eventvar == consdata->posFlowBinvar )
         consdata->propvarind |= PIPEODE_POS_LB;
      else if ( eventvar == consdata->negFlowBinvar )
         consdata->propvarind |= PIPEODE_NEG_LB;
      break;

   case SCIP_EVENTTYPE_UBRELAXED:
      if ( eventvar == consdata->p_in_var )
      {
         if ( consdata->propvarind & PIPEODE_PIN_UB )
            consdata->propvarind ^= PIPEODE_PIN_UB;
      }
      else if ( eventvar == consdata->p_out_var )
      {
         if ( consdata->propvarind & PIPEODE_POUT_UB )
            consdata->propvarind ^= PIPEODE_POUT_UB;
      }
      else if ( eventvar == consdata->q_var )
      {
         if ( consdata->propvarind & PIPEODE_Q_UB )
            consdata->propvarind ^= PIPEODE_Q_UB;
      }
      else if ( eventvar == consdata->posFlowBinvar )
      {
         if ( consdata->propvarind & PIPEODE_POS_UB )
            consdata->propvarind ^= PIPEODE_POS_UB;
      }
      else if ( eventvar == consdata->negFlowBinvar )
      {
         if ( consdata->propvarind & PIPEODE_NEG_UB )
            consdata->propvarind ^= PIPEODE_NEG_UB;
      }
      break;

   case SCIP_EVENTTYPE_UBTIGHTENED:

      SCIPdebugMessage("<%s>: upper bound of p_in|p_out|q|posFlowBinvar|negFlowBinvar tightened.\n", SCIPconsGetName(cons));

      if ( eventvar == consdata->p_in_var )
         consdata->propvarind |= PIPEODE_PIN_UB;
      else if ( eventvar == consdata->p_out_var )
         consdata->propvarind |= PIPEODE_POUT_UB;
      else if ( eventvar == consdata->q_var )
         consdata->propvarind |= PIPEODE_Q_UB;
      else if ( eventvar == consdata->posFlowBinvar )
         consdata->propvarind |= PIPEODE_POS_UB;
      else if ( eventvar == consdata->negFlowBinvar )
         consdata->propvarind |= PIPEODE_NEG_UB;
      break;

   default:
      SCIPerrorMessage("invalid event type.\n");
      SCIPABORT();
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}


/**
 *  Callback methods of constraint handler
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyPipeODE)
{  /*lint --e{715}*/
  assert(scip != NULL);
  assert(conshdlr != NULL);
  assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
  assert(valid != NULL);

  SCIP_CALL( SCIPincludeConshdlrPipeODE(scip) );

  *valid = TRUE;

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreePipeODE)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeMemory(scip, &conshdlrdata);

   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolPipeODE)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_NLROW* nlrow;
   SCIP_EXPR* odeexpr;
   SCIP_EXPR* varexprs[3];
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
   SCIP_EXPRDATA* data;
#else
   SCIP_EXPRTREE* exprtree;
   SCIP_USEREXPRDATA* data;
   SCIP_VAR* exprvars[5];
#endif
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   /* exit if no NLP is present or we do not want to use the representation */
   if ( !SCIPisNLPConstructed(scip) || ! conshdlrdata->nlp_representation )
      return SCIP_OKAY;

#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
   /* set up user expressions for all pipes that are used to evaluate the pipe equation */
   for (i = 0; i < nconss; ++i)
   {
      int nexpr = 3;

      consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );

      /* define varexpressions */
      SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[0], consdata->p_in_var, NULL, NULL) );
      SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[1], consdata->p_out_var, NULL, NULL) );
      SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[2], consdata->q_var, NULL, NULL) );

      /* define data */
      SCIP_CALL( SCIPallocBlockMemory(scip, &data) );
      data->scip = scip;
      data->conshdlrdata = conshdlrdata;
      data->consdata = consdata;
      data->curvature = SCIP_EXPRCURV_UNKNOWN;

      /* create user expression */
      SCIP_CALL( SCIPcreateExpr(scip, &odeexpr, conshdlrdata->exprhdlr, data, nexpr, varexprs, NULL, NULL) );

      /* add convex relaxation or equation */
      if ( conshdlrdata->nlp_convex )
      {
         /* possibly update curvature */
         if ( SCIPisGE(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var), 0.0) || SCIPisLE(scip, SCIPcomputeVarUbLocal(scip, consdata->q_var), 0.0) )
            data->curvature = SCIP_EXPRCURV_CONVEX;
         else if ( SCIPisZero(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var)) && SCIPisZero(scip, SCIPcomputeVarUbLocal(scip, consdata->q_var)) )
            data->curvature = SCIP_EXPRCURV_LINEAR;

         SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0, 0, NULL, NULL, odeexpr, - SCIPinfinity(scip), 0.0, data->curvature) );
      }
      else
      {
         /* add equation */
         SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0, 0, NULL, NULL, odeexpr, 0.0, 0.0, SCIP_EXPRCURV_UNKNOWN) );
      }
      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );

      SCIP_CALL( SCIPreleaseExpr(scip, &odeexpr) );
      SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[2]) );
      SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[1]) );
      SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[0]) );

      /* release row */
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   }

   /* ------------------------------------------------- */
#else
   /* ------------------------------------------------- */

   /* set up user expressions for all pipes that are used to evaluate the pipe equation */
   for (i = 0; i < nconss; ++i)
   {
      int nexpr = 3;

      consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );

      /* set variable array */
      exprvars[0] = consdata->p_in_var;
      exprvars[1] = consdata->p_out_var;
      exprvars[2] = consdata->q_var;

      /* define varexpressions */
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &varexprs[0], SCIP_EXPR_VARIDX, 0) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &varexprs[1], SCIP_EXPR_VARIDX, 1) );
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &varexprs[2], SCIP_EXPR_VARIDX, 2) );

      /* define data */
      SCIP_CALL( SCIPallocBlockMemory(scip, &data) );
      data->scip = scip;
      data->conshdlrdata = conshdlrdata;
      data->consdata = consdata;
      data->curvature = SCIP_EXPRCURV_UNKNOWN;

      /* create user expression */
      SCIP_CALL( SCIPexprCreateUser(SCIPblkmem(scip), &odeexpr, nexpr, varexprs, data,
            SCIP_EXPRINTCAPABILITY_FUNCVALUE | SCIP_EXPRINTCAPABILITY_GRADIENT,    /* capability */
            SCIPuserexprEvalODE,                                                   /* eval callback */
            NULL,                                                                  /* interval eval callback */
            SCIPuserexprCurvODE,                                                   /* curvature callback */
            NULL,                                                                  /* propagation callback */
            NULL,                                                                  /* estimate callback */
            SCIPuserexprCopyODE,                                                   /* copydata callback */
            SCIPuserexprFreeODE,                                                   /* freedata callback */
            NULL                                                                   /* print callback */
            ) );

      /**
       *  TODO: free/null varexpressions??
       */

      /* create expression tree */
      SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, odeexpr, nexpr, 0, NULL) );
      SCIP_CALL( SCIPexprtreeSetVars(exprtree, nexpr, exprvars) );

      /* add convex relaxation or equation */
      if ( conshdlrdata->nlp_convex )
      {
         /* possibly update curvature */
         if ( SCIPisGE(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var), 0.0) || SCIPisLE(scip, SCIPcomputeVarUbLocal(scip, consdata->q_var), 0.0) )
            data->curvature = SCIP_EXPRCURV_CONVEX;
         else if ( SCIPisZero(scip, SCIPcomputeVarLbLocal(scip, consdata->q_var)) && SCIPisZero(scip, SCIPcomputeVarUbLocal(scip, consdata->q_var)) )
            data->curvature = SCIP_EXPRCURV_LINEAR;

         SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0, 0, NULL, NULL, 0, NULL, 0, NULL, exprtree, - SCIPinfinity(scip), 0.0, data->curvature) );
      }
      else
      {
         /* add equation */
         SCIP_CALL( SCIPcreateNlRow(scip, &nlrow, SCIPconsGetName(conss[i]), 0.0, 0, NULL, NULL, 0, NULL, 0, NULL, exprtree, 0.0, 0.0, data->curvature) );
      }
      SCIP_CALL( SCIPaddNlRow(scip, nlrow) );

      SCIP_CALL( SCIPexprtreeFree(&exprtree) );

      /* release row */
      SCIP_CALL( SCIPreleaseNlRow(scip, &nlrow) );
   }
#endif
   return SCIP_OKAY;
}


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolPipeODE)
{  /*lint --e{715}*/
#ifdef SCIP_DEBUG
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   SCIP_SOL* sol;
   SCIP_Real Ma;
   SCIP_Real q;
   SCIP_Real p;
   int i;
   int pos_q = 0;
   int neg_q = 0;
   SCIP_Real maxma = 0;

   if ( restart || SCIPgetSubscipDepth(scip) != 0 )
      return SCIP_OKAY;

   sol = SCIPgetBestSol(scip);
   if ( sol == NULL )
      return SCIP_OKAY;

   for (i = 0; i < nconss; ++i)
   {
      cons = conss[i];
      consdata = SCIPconsGetData(cons);

      p = MIN(SCIPgetSolVal(scip, sol, consdata->p_in_var), SCIPgetSolVal(scip, sol, consdata->p_out_var) );
      q = SCIPgetSolVal(scip, sol, consdata->q_var);
      if ( SCIPisFeasPositive(scip, q) )
         ++pos_q;
      else if ( SCIPisFeasNegative(scip, q) )
         ++neg_q;

      q = ABS(q);

      Ma = (consdata->c * q) / (consdata->A * p) / 1e5;

      SCIPinfoMessage(scip, NULL, "Machnumber of pipe <%s>: \t %f\n", SCIPconsGetName(cons), Ma);
      maxma = MAX(maxma, Ma);
   }
   SCIPinfoMessage(scip, NULL, "Maximum Machnumber: \t %f\n", maxma);

   SCIPinfoMessage(scip, NULL, "\nTotal number of pipes: %d\n", nconss);
   SCIPinfoMessage(scip, NULL, "Number of pipes with positive flow: %d\n", pos_q);
   SCIPinfoMessage(scip, NULL, "Number of pipes with negative flow: %d\n", neg_q);
#endif

   return SCIP_OKAY;
}

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeletePipeODE)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_PROBDATA* probdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( consdata != NULL);
   assert( *consdata != NULL);

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIPdebugMessage("Deleting PipeODE constraint <%s>.\n", SCIPconsGetName(cons));

#ifdef STEPSIZE_DEBUG
   if ( ! SCIPconsIsOriginal(cons) && (SCIPgetSubscipDepth(scip) == 0) )
   {
      SCIP_Real stepsize = (*consdata)->L / (double)(*consdata)->N;
      SCIPinfoMessage(scip, NULL, "Final discretization on pipe <%s>:\t #%d \t %f m\n", SCIPconsGetName(cons), (*consdata)->N, stepsize);
   }
#endif

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* drop events on variable */
   if ( SCIPvarIsTransformed((*consdata)->q_var) )
   {
      if ( ! probdata->noFlowBinvars )
      {
         if ( !SCIPvarIsDeleted((*consdata)->negFlowBinvar) )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->negFlowBinvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)cons, -1) ); /*lint !e740*/
         }
         if ( !SCIPvarIsDeleted((*consdata)->posFlowBinvar) )
         {
            SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->posFlowBinvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)cons, -1) ); /*lint !e740*/
         }
      }
      if ( !SCIPvarIsDeleted((*consdata)->q_var) )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->q_var, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)cons, -1) ); /*lint !e740*/
      }
      if ( !SCIPvarIsDeleted((*consdata)->p_out_var) )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->p_out_var, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)cons, -1) ); /*lint !e740*/
      }
      if ( !SCIPvarIsDeleted((*consdata)->p_in_var) )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, (*consdata)->p_in_var, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*)cons, -1) ); /*lint !e740*/
      }
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransPipeODE)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_PROBDATA* probdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONSDATA* sourcedata;
   char s[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( sourcecons != NULL );
   assert( targetcons != NULL );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIPdebugMessage("transforming PipeODE constraint <%s>.\n", SCIPconsGetName(sourcecons) );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   /* transform variables */
   assert( sourcedata->p_in_var != NULL );
   assert( sourcedata->p_out_var != NULL );
   assert( sourcedata->q_var != NULL );

   SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->p_in_var, &(consdata->p_in_var)) );
   SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->p_out_var, &(consdata->p_out_var)) );
   SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->q_var, &(consdata->q_var)) );
   if ( probdata->noFlowBinvars )
   {
      consdata->negFlowBinvar = NULL;
      consdata->posFlowBinvar = NULL;
   }
   else
   {
      assert( sourcedata->negFlowBinvar != NULL );
      assert( sourcedata->posFlowBinvar != NULL );

      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->negFlowBinvar, &(consdata->negFlowBinvar)) );
      SCIP_CALL( SCIPgetTransformedVar(scip, sourcedata->posFlowBinvar, &(consdata->posFlowBinvar)) );
   }

   consdata->N = sourcedata->N;
   consdata->A = sourcedata->A;
   consdata->Asq = sourcedata->Asq;
   consdata->L = sourcedata->L;
   consdata->D = sourcedata->D;
   consdata->lambda = sourcedata->lambda;
   consdata->lcsqd = sourcedata->lcsqd;
   consdata->c = sourcedata->c;
   consdata->csq = sourcedata->csq;

   /* set default values to "new constraint" */
   consdata->nr_cave        = 0;
   consdata->nr_vex         = 0;
   consdata->feasible       = FALSE;
   consdata->violation      = SCIP_INVALID;
   consdata->type_violation = -1;
   consdata->propvarind     = 0;

   /* create constraint */
   (void) SCIPsnprintf(s, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));

   SCIP_CALL( SCIPcreateCons(scip, targetcons, s, conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* catch bound change events of variable */
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->p_in_var, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*) *targetcons, NULL) ); /*lint !e740*/
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->p_out_var, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*) *targetcons, NULL) ); /*lint !e740*/
   SCIP_CALL( SCIPcatchVarEvent(scip, consdata->q_var, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*) *targetcons, NULL) ); /*lint !e740*/
   if ( ! probdata->noFlowBinvars )
   {
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->posFlowBinvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*) *targetcons, NULL) ); /*lint !e740*/
      SCIP_CALL( SCIPcatchVarEvent(scip, consdata->negFlowBinvar, SCIP_EVENTTYPE_BOUNDCHANGED, conshdlrdata->eventhdlr, (SCIP_EVENTDATA*) *targetcons, NULL) ); /*lint !e740*/
   }

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpPipeODE)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   int i;
   SCIP_Real q_var_ub;
   SCIP_Real q_var_lb;
   SCIP_Real p_var_ub;
   SCIP_Real p_var_lb;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* loop through all constraints */
   for (i = nconss - 1; i >= 0; --i)
   {
      cons = conss[i];
      assert( cons != NULL );

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      q_var_lb = SCIPcomputeVarLbLocal(scip, consdata->q_var);
      q_var_ub = SCIPcomputeVarUbLocal(scip, consdata->q_var);

#ifdef FLOW_TIGHTENING_DEBUG
      if ( SCIPgetSubscipDepth(scip) == 0 )
      {
         SCIPinfoMessage(scip, NULL, "Flow Bounds on Pipe <%s>:\t LB %8.4f \t UB %8.4f\n", SCIPconsGetName(cons), q_var_lb, q_var_ub);
      }
#endif

      if ( ! SCIPisFeasNegative(scip, q_var_lb) )
      {
         p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_out_var);
         p_var_lb = SCIPcomputeVarLbLocal(scip, consdata->p_out_var);
      }
      else if ( !SCIPisFeasPositive(scip, q_var_ub) )
      {
         p_var_ub = SCIPcomputeVarUbLocal(scip, consdata->p_in_var);
         p_var_lb = SCIPcomputeVarLbLocal(scip, consdata->p_in_var);
      }
      else
      {
         /* Direction of the flow not uniquely determined.
          * The theory for this case has yet to be developed. */
         continue;
      }

      if ( SCIPisFeasEQ(scip, q_var_ub, q_var_lb) && SCIPisEQ(scip, p_var_lb, p_var_ub) )
      {
         /* TODO: implement some method
          * if the third variable is fixed too, delete the constraint
          * else the other pressure variable should get fixed too and the constraint deleted */
      }
      else
      {
         SCIP_CALL( generateCave(scip, conshdlr, cons, NULL, NULL, NULL, NULL) );
         /* TODO: add gradient cuts? */
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpPipeODE)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   int i = 0;
   int j = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeMostViolatedConstraint(scip, conshdlr, conss, nconss, NULL, &i) );

   if ( i == -1 )
      return SCIP_OKAY;

   SCIPdebugMessage("Most violated constraint is number %d\n",i+1);

   cons = conss[i];
   assert( cons != NULL );

   SCIP_CALL( separateODE(scip, conshdlr, cons, NULL, result) );

   while ( *result == SCIP_DIDNOTFIND )
   {
      if ( j == i )
         ++j;
      if ( j >= nconss )
         break;

      cons = conss[j];
      assert( cons != NULL );
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      if ( consdata->feasible == TRUE )
      {
         ++j;
         continue;
      }
      else
      {
         SCIP_CALL( separateODE(scip, conshdlr, cons, NULL, result) );
         ++j;
      }
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolPipeODE)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   int i = 0;
   int j = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( computeMostViolatedConstraint(scip, conshdlr, conss, nconss, sol, &i) );

   if ( i == -1 )
      return SCIP_OKAY;

   /* since i != -1 the lp solution is infeasible! */
   *result = SCIP_DIDNOTFIND;

   cons = conss[i];
   assert( cons != NULL );

   SCIP_CALL( separateODE(scip, conshdlr, cons, sol, result) );

   while ( *result == SCIP_DIDNOTFIND )
   {
      if ( j == i )
         ++j;
      if ( j >= nconss )
         break;

      cons = conss[j];
      assert( cons != NULL );
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      if ( consdata->feasible == TRUE )
      {
         ++j;
         continue;
      }
      else
      {
         SCIP_CALL( separateODE(scip, conshdlr, cons, sol, result) );
         ++j;
      }
   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpPipeODE)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   int i = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   *result = SCIP_FEASIBLE;

   /* we proceed even if solinfeasible is true, because the reductions might pay off (the node is not yet cut off) */

   /* check the solution, if it is infeasible, return the pipe with the biggest deviation of the computed bounds and the corresponding solution value */
   SCIP_CALL( computeMostViolatedConstraint(scip, conshdlr, conss, nconss, NULL, &i) );

   if ( i != -1 )
   {
      /* since i != -1 the lp solution is infeasible! */
      *result = SCIP_INFEASIBLE;

      cons = conss[i];
      assert( cons != NULL );

      SCIPclearExternBranchCands(scip);
      SCIP_CALL( enforceODE(scip, conshdlr, cons, NULL, result) );
   }

   if ( *result == SCIP_INFEASIBLE )
   {
      SCIPdebugMessage( "EnforceODE added %d extern branching candidates.\n",SCIPgetNExternBranchCands(scip));
   }

   /* TODO: Do something if no cut was found and no branching candidates are added! */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsPipeODE)
{  /*lint --e{715}*/
   SCIP_CONS* cons;
   SCIP_CONSDATA* consdata;
   int i = 0;
   int j = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   if ( objinfeasible )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   *result = SCIP_FEASIBLE;

   /* check the solution, if it is infeasible, return the pipe with the biggest deviation of the computed bounds and the corresponding solution value */
   SCIP_CALL( computeMostViolatedConstraint(scip, conshdlr, conss, nconss, NULL, &i) );

   if ( i != -1 )
   {
      /* since i != -1 the lp solution is infeasible! */
      *result = SCIP_INFEASIBLE;

      cons = conss[i];
      assert( cons != NULL );

      SCIPclearExternBranchCands(scip);
      SCIP_CALL( enforcePS(scip, conshdlr, cons, NULL, result) );
   }

   /* if we did not find a domain reduction or a branching candidate for the most violated pipe, we check the others */
   while ( (*result == SCIP_INFEASIBLE) && (SCIPgetNExternBranchCands(scip) == 0) && (j < nconss) )
   {
      if ( j == i )
      {
         ++j;
         continue;
      }

      cons = conss[j];
      assert( cons != NULL );

      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL );

      if ( ! consdata->feasible )
      {
         SCIP_CALL( enforcePS(scip, conshdlr, cons, NULL, result) );
      }
      ++j;
   }

   if ( *result == SCIP_INFEASIBLE )
   {
      if ( SCIPgetNExternBranchCands(scip) == 0 )
      {
         SCIPerrorMessage("EnforcePS could not find any domain reduction or branching candidate.\n");
      }
      else
      {
         SCIPdebugMessage("EnforcePS added %d extern branching candidates.\n", SCIPgetNExternBranchCands(scip));
      }
   }
   else if ( *result == SCIP_REDUCEDDOM )
   {
      SCIPdebugMessage("EnforcePS found a domain reduction.\n");
   }
   else if ( *result == SCIP_CUTOFF )
   {
      SCIPdebugMessage("EnforcePS found an infeasible node.\n");
   }

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckPipeODE)
{  /*lint --e{715}*/
   SCIP_CONS* cons = NULL;
   int i = 0;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   *result = SCIP_FEASIBLE;

   /* loop through all constraints */
   while ( i < nconss )
   {
      cons = conss[i];
      assert( cons != NULL );

      SCIP_CALL( checkSolution(scip, conshdlr, cons, sol, result) );

      if ( *result == SCIP_INFEASIBLE )
         break;

      /* check next pipe/constraint */
      ++i;
   }

   if ( printreason && (*result == SCIP_INFEASIBLE) )
   {
      SCIP_CONSDATA* consdata;

      assert( cons != NULL );
      consdata = SCIPconsGetData(cons);

      SCIPinfoMessage(scip, NULL, "[pipeODE] constraint <%s> is violated\n", SCIPconsGetName(cons));

      switch( consdata->type_violation )
      {
      case FEASIBLE:
         SCIPerrorMessage("Error in consCheckPipeODE: result is SCIP_INFEASIBLE, but constraint is not violated!\n");
         break;
      case MACHBOUND:
         SCIPinfoMessage(scip, NULL, "   violation: Output pressure is %f <bar> too small for the mach-bound.\n", consdata->violation);
         break;
      case PIN_TOO_BIG:
         SCIPinfoMessage(scip, NULL, "   violation: Input pressure is %f <bar> too big.\n", consdata->violation);
         break;
      case PIN_TOO_LOW:
         SCIPinfoMessage(scip, NULL, "   violation: Input pressure is %f <bar> too small.\n", consdata->violation);
         break;
      default:
         SCIPerrorMessage("Unknown violation type!\n");
         return SCIP_INVALIDDATA;
      }
   }

   /* FEASIBLE = -1,                                 /\**< solution is feasible *\/ */
   /* MACHBOUND = 1,                                 /\**< machbound is not fulfilled at the output *\/ */
   /* PIN_TOO_BIG = 2,                               /\**< the input pressure is too big *\/ */
   /* PIN_TOO_LOW = 3     */

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropPipeODE)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   int i;
   int counter = 0;
   SCIP_VARSTATUS varstatus;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   conshdlrdata =  SCIPconshdlrGetData(conshdlr);

   *result = SCIP_DIDNOTFIND;

   for (i = nconss -1; i >= 0 ; --i)
   {
      cons = conss[i];
      assert( cons != NULL );

      /* get data of constraint */
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL);

      varstatus = SCIPvarGetStatus(consdata->q_var);

      /* do not propagate, if no variable bound has changed */
      if ( consdata->propvarind == PIPEODE_NONE )
         continue;

      SCIP_CALL( machBound(scip, conshdlrdata, cons, result, &counter) );
      if ( *result == SCIP_CUTOFF )
         return SCIP_OKAY;

      if ( conshdlrdata->with_flowTightening )
      {
         if ( (varstatus != SCIP_VARSTATUS_MULTAGGR) && (varstatus != SCIP_VARSTATUS_FIXED) )
         {
            SCIP_CALL( flowTightening(scip, conshdlrdata, cons, &counter) );
         }
      }

      SCIP_CALL( boundPropagation(scip, conshdlrdata, cons, result, &counter) );

      /* reset marker for changed variable bounds */
      consdata->propvarind = 0;

      if ( *result == SCIP_CUTOFF )
         return SCIP_OKAY;
   }

   if ( counter > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolPipeODE)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_PROBDATA* probdata;
   SCIP_CONSDATA* consdata;
   SCIP_VARSTATUS varstatus;
   SCIP_CONS* cons;
   int bdchgs = 0;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   *result = SCIP_DIDNOTFIND;

   /* loop through all constraints */
   for (i = nconss-1; i >= 0; --i)
   {
      cons = conss[i];
      assert( cons != NULL );

      /* get data of constraint */
      consdata = SCIPconsGetData(cons);
      assert( consdata != NULL);

      /* if nothing changed then skip this constraint */
      if ( consdata->propvarind != PIPEODE_NONE )
      {
         varstatus = SCIPvarGetStatus(consdata->q_var);

         /* CALL machBound */
         SCIP_CALL( machBound(scip, conshdlrdata, cons, result, &bdchgs) );

         if ( (*result) == SCIP_CUTOFF )
         {
            *nchgbds += bdchgs;
            return SCIP_OKAY;
         }

         if ( conshdlrdata->with_flowTightening )
         {
            if ( (varstatus != SCIP_VARSTATUS_FIXED) && (varstatus != SCIP_VARSTATUS_MULTAGGR) )
            {
               /* CALL flowTightening  */
               SCIP_CALL( flowTightening(scip, conshdlrdata, cons, &bdchgs) );
            }
         }

         /* CALL boundPropagation */
         SCIP_CALL( boundPropagation(scip, conshdlrdata, cons, result, &bdchgs) );

         if ( (*result) == SCIP_CUTOFF )
         {
            *nchgbds += bdchgs;
            if ( SCIPgetSubscipDepth(scip) == 0 )
               SCIPinfoMessage(scip, NULL, "Bound propagation detected infeasibility while propagating bounds on pipe <%s> in presolving.\n", SCIPconsGetName(cons));
            return SCIP_OKAY;
         }
      }
   }

   /* Perform OBBT only once per presolving, since it is time expensive.
    * Furthermore, sometimes an error occurs when performing obbt after restarts.
    */
   if ( conshdlrdata->with_obbt && SCIPgetNRuns(scip) == 1 )
   {
      if ( presoltiming == SCIP_PRESOLTIMING_FINAL && (SCIPgetSubscipDepth(scip) == 0) )
      {
         SCIP_STATUS status      = SCIP_STATUS_UNKNOWN;
         int         obbt_bdchgs = 0;

         SCIP_CALL( obbt( scip, conshdlrdata, &status, &obbt_bdchgs) );

         bdchgs += obbt_bdchgs;

         if ( status == SCIP_STATUS_USERINTERRUPT )
         {
            if ( bdchgs > 0 )
            {
               *nchgbds += bdchgs;
               *result = SCIP_SUCCESS;
            }
            SCIPtryTerminate();
            return SCIP_OKAY;
         }
         else if ( status == SCIP_STATUS_INFEASIBLE )
         {
            *result = SCIP_CUTOFF;
            SCIPinfoMessage(scip, NULL, "OBBT detected infeasibility in run %d.\n", SCIPgetNRuns(scip));
            *nchgbds += bdchgs;
            return SCIP_OKAY;
         }
      }
   }

   *nchgbds += bdchgs;

   if ( presoltiming == SCIP_PRESOLTIMING_FINAL && (SCIPgetSubscipDepth(scip) == 0) )
   {
      SCIP_Real qlb, qub;
      SCIP_RESULT outcome;

      for (i = nconss - 1; i >= 0; --i)
      {
         cons = conss[i];
         assert( cons != NULL );

         outcome = SCIP_DIDNOTRUN;

         /* get data of constraint */
         consdata = SCIPconsGetData(cons);
         assert( consdata != NULL);

         qlb = SCIPcomputeVarLbLocal(scip, consdata->q_var);
         qub = SCIPcomputeVarUbLocal(scip, consdata->q_var);

         if ( ! probdata->noFlowBinvars && flowIsZero(qlb, conshdlrdata->flow_feastol) && flowIsZero( qub, conshdlrdata->flow_feastol) )
         {
            SCIP_CALL( removeODECons(scip, conshdlr, cons, &outcome, nfixedvars, naggrvars, ndelconss) );
         }
         if ( outcome == SCIP_SUCCESS )
         {
            *result = SCIP_SUCCESS;
         }
         else if ( outcome == SCIP_CUTOFF )
         {
            *result = SCIP_CUTOFF;
            SCIPdebugMessage("removeODECons returned SCIP_CUTOFF for pipe <%s>.\n", SCIPconsGetName(cons));
            return SCIP_OKAY;
         }
      }
   }

   if ( presoltiming == SCIP_PRESOLTIMING_FINAL )
   {
      if ( SCIPisIpoptAvailableIpopt() && conshdlrdata->nlp_representation )
      {
         /* enable nlp */
         SCIPenableNLP(scip);
      }
   }

   if ( bdchgs > 0 )
      *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropPipeODE)
{  /*lint --e{715}*/
    SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infervar != NULL );

   consdata = SCIPconsGetData(cons);
   *result = SCIP_DIDNOTFIND;

   switch ( inferinfo )
   {
   case FPROP_UPPER:
      if ( infervar == consdata->p_out_var )
      {
         /* printf("Bound change reason FPROP_UPPER. Added flow and pressure-in variable to conflict.\n"); */
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->q_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_in_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   case FPROP_LOWER:
      if ( infervar == consdata->p_out_var )
      {
         /* printf("Bound change reason FPROP_LOWER. Added flow and pressure-in variable to conflict.\n"); */
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->q_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->p_in_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else if ( infervar == consdata->q_var )
      {
         /* printf("Bound change reason FPROP_LOWER. Added the pressure variables to conflict.\n"); */
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_out_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->p_in_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   case BPROP_UPPER:
      if ( infervar == consdata->p_in_var )
      {
         /* printf("Bound change reason BPROP_UPPER. Added flow and pressure-out variable to conflict.\n"); */
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->q_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_out_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   case BPROP_LOWER:
      if ( infervar == consdata->p_in_var )
      {
         /* printf("Bound change reason BPROP_LOWER. Added flow and pressure-out variable to conflict.\n"); */
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->q_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->p_out_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else if ( infervar == consdata->q_var )
      {
         /* printf("Bound change reason BPROP_LOWER. Added the pressure variables to conflict.\n"); */
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_in_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->p_out_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   case PIN_GT_POUT:
      if ( infervar == consdata->q_var )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_out_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->p_in_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else if ( infervar == consdata->negFlowBinvar )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_out_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->p_in_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   case PIN_LT_POUT:
      if ( infervar == consdata->q_var )
      {
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->p_out_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_in_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else if ( infervar == consdata->posFlowBinvar )
      {
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->p_out_var, bdchgidx) );
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_in_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   case MACH_BOUND:
      if ( infervar == consdata->q_var )
      {
         if ( SCIPgetVarUbAtIndex(scip, consdata->p_in_var, bdchgidx, FALSE) <= SCIPgetVarUbAtIndex(scip, consdata->p_out_var, bdchgidx, FALSE) )
         {
            SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_in_var, bdchgidx) );
         }
         else
         {
            SCIP_CALL( SCIPaddConflictUb(scip, consdata->p_out_var, bdchgidx) );
         }
         *result = SCIP_SUCCESS;
      }
      else if ( infervar == consdata->p_in_var )
      {
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->q_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else if ( infervar == consdata->p_out_var )
      {
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->q_var, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable does not belong to the given constraint.\n");
      }
      break;

   case POS_LB:
      if ( infervar == consdata->q_var )
      {
         /* this case should only occure if the noFlowBinvars option is not used */
         assert( consdata->posFlowBinvar != NULL );

         SCIP_CALL( SCIPaddConflictLb(scip, consdata->posFlowBinvar, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   case POS_UB:
      if ( infervar == consdata->q_var )
      {
         /* this case should only occur if the noFlowBinvars option is not used */
         assert( consdata->posFlowBinvar != NULL );

         SCIP_CALL( SCIPaddConflictUb(scip, consdata->posFlowBinvar, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   case NEG_LB:
      if ( infervar == consdata->q_var )
      {
         /* this case should only occur if the noFlowBinvars option is not used */
         assert( consdata->negFlowBinvar != NULL );

         SCIP_CALL( SCIPaddConflictLb(scip, consdata->negFlowBinvar, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   case NEG_UB:
      if ( infervar == consdata->q_var )
      {
         /* this case should only occur if the noFlowBinvars option is not used */
         assert( consdata->negFlowBinvar != NULL );

         SCIP_CALL( SCIPaddConflictUb(scip, consdata->negFlowBinvar, bdchgidx) );
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPerrorMessage("Changed variable and reason do not add up.\n");
      }
      break;

   default:
      SCIPerrorMessage("Not implemented conflict type:\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockPipeODE)
{  /*lint --e{715} */
   SCIP_PROBDATA* probdata;
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   SCIPdebugMessage("Locking PipeODE constraint <%s>.\n", SCIPconsGetName(cons));

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);

   /* the constaint may be violated in any way */
   SCIP_CALL( SCIPaddVarLocks(scip, consdata->q_var, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   SCIP_CALL( SCIPaddVarLocks(scip, consdata->p_in_var, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   SCIP_CALL( SCIPaddVarLocks(scip, consdata->p_out_var, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   if ( ! probdata->noFlowBinvars )
   {
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->posFlowBinvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );
      SCIP_CALL( SCIPaddVarLocks(scip, consdata->negFlowBinvar, nlockspos + nlocksneg, nlockspos + nlocksneg) );
   }

   return SCIP_OKAY;
}

/** constraint display method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRINT(consPrintPipeODE)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of PipeODE constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPrintPipeODE NULL
#endif


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyPipeODE)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SCIP_CONSDATA* sourcedata;
   SCIP_Bool success;

   assert( scip != NULL );
   assert( sourcescip != NULL );
   assert( valid != NULL );

   SCIPdebugMessage("copying PipeODE constraint <%s>.\n", SCIPconsGetName(sourcecons) );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   /* get data of original constraint */
   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL);

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   /* transform variables */
   assert( sourcedata->p_in_var != NULL );
   assert( sourcedata->p_out_var != NULL );
   assert( sourcedata->q_var != NULL );

   SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->p_in_var, &(consdata->p_in_var), varmap, consmap, global, &success) );
   assert( success );

   SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->p_out_var, &(consdata->p_out_var), varmap, consmap, global, &success) );
   assert( success );

   SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->q_var, &(consdata->q_var), varmap, consmap, global, &success) );
   assert( success );

   if ( probdata->noFlowBinvars )
   {
      consdata->posFlowBinvar = NULL;
      consdata->negFlowBinvar = NULL;
   }
   else
   {
      assert( sourcedata->negFlowBinvar != NULL );
      assert( sourcedata->posFlowBinvar != NULL );

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->negFlowBinvar, &(consdata->negFlowBinvar), varmap, consmap, global, &success) );
      assert( success );

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->posFlowBinvar, &(consdata->posFlowBinvar), varmap, consmap, global, &success) );
      assert( success );
   }

   consdata->N = sourcedata->N;
   consdata->A = sourcedata->A;
   consdata->Asq = sourcedata->Asq;
   consdata->L = sourcedata->L;
   consdata->D = sourcedata->D;
   consdata->lambda = sourcedata->lambda;
   consdata->lcsqd = sourcedata->lcsqd;
   consdata->c = sourcedata->c;
   consdata->csq = sourcedata->csq;

   /* set default values to "new constraint" */
   consdata->nr_cave = 0;
   consdata->nr_vex = 0;
   consdata->feasible = FALSE;
   consdata->violation = SCIP_INVALID;
   consdata->type_violation = -1;
   consdata->propvarind = 0;

   /* find the PipeODE constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("PipeOde constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, SCIPconsGetName(sourcecons), conshdlr, consdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons),
         SCIPconsIsEnforced(sourcecons), SCIPconsIsChecked(sourcecons),
         SCIPconsIsPropagated(sourcecons), SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons),
         SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
#if 0
static
SCIP_DECL_CONSPARSE(consParsePipeODE)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of PipeODE constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consParsePipeODE NULL
#endif


/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsPipeODE)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( cons != NULL );
   assert( success != NULL );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL);

   /* get data of constraint */
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL);

   SCIPdebugMessage("Locking PipeODE constraint <%s>.\n", SCIPconsGetName(cons));

   *success = FALSE;

   if ( probdata->noFlowBinvars )
   {
      if ( varssize >= 3 )
      {
         vars[0] = consdata->p_in_var;
         vars[1] = consdata->p_out_var;
         vars[2] = consdata->q_var;
         *success = TRUE;
      }
   }
   else
   {
      if ( varssize >= 5 )
      {
         vars[0] = consdata->p_in_var;
         vars[1] = consdata->p_out_var;
         vars[2] = consdata->q_var;
         vars[3] = consdata->posFlowBinvar;
         vars[4] = consdata->negFlowBinvar;
         *success = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsPipeODE)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;

   assert( success != NULL );
   assert( nvars != NULL );

   probdata = SCIPgetProbData(scip);

   if ( probdata->noFlowBinvars )
      *nvars = 3;
   else
      *nvars = 5;

   *success = TRUE;

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for PipeODE constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrPipeODE(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->ngradientcuts       = 0;
   conshdlrdata->machbound           = 2;
   conshdlrdata->pressure_feastol    = 0.01;
   conshdlrdata->flow_feastol        = 0.01;
   conshdlrdata->offset              = 0.001;
   conshdlrdata->eventhdlr           = NULL;
   conshdlrdata->with_flowTightening = TRUE;
   conshdlrdata->with_obbt           = TRUE;
   conshdlrdata->obbt_mingap         = 10.0;

   /* create event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecPipeODE, NULL) );
   if ( conshdlrdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for pipe ODE constraints not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpPipeODE, consEnfopsPipeODE, consCheckPipeODE, consLockPipeODE,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyPipeODE, consCopyPipeODE) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeletePipeODE) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolPipeODE) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreePipeODE) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsPipeODE) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsPipeODE) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolPipeODE) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpPipeODE) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParsePipeODE) );
#if SCIP_VERSION >= 320
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolPipeODE, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
#else
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolPipeODE, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
#endif
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintPipeODE) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropPipeODE, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropPipeODE) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpPipeODE, consSepasolPipeODE, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransPipeODE) );

   /* include expression handler */
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
   {
      SCIP_EXPRHDLR* exprhdlr;

      SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
            EXPRHDLR_PRECEDENCE, SCIPexprEvalODE, NULL) );
      assert( exprhdlr != NULL );

      SCIPexprhdlrSetCurvature(exprhdlr, SCIPexprCurvODE);
      SCIPexprhdlrSetDiff(exprhdlr, SCIPexprBWODE, SCIPexprFWdiffODE, NULL);
      SCIPexprhdlrSetCopyFreeData(exprhdlr, SCIPexprCopyODE, SCIPexprFreeODE);
      conshdlrdata->exprhdlr = exprhdlr;
   }
#endif

   /* define new feasibility parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "pressure_feastol", "tolerance used for comparison of pressure values",
         &(conshdlrdata->pressure_feastol), FALSE, 1e-2, 1e-4, 1.0, NULL, NULL) );
   /* default value for comparing pressure: 1e-2 bar
    * maximal tolerance: 1 bar
    * minimal tolerance: 1e-4 bar
    */
   SCIP_CALL( SCIPaddRealParam(scip, "flow_feastol", "tolerance used for comparison of flow values",
         &(conshdlrdata->flow_feastol), FALSE, 1e-2, 1e-4, 1.0, NULL, NULL) );
   /* default value for comparing pressure: 1e-2 kg/s
    * maximal tolerance: 1 kg/s
    * minimal tolerance: 1e-4 kg/s
    */
   SCIP_CALL( SCIPaddIntParam(scip, "machbound", "numerator used for bound on mach number: v/c = machbound/5",
         &(conshdlrdata->machbound), FALSE, 2, 1, 4, NULL, NULL) );
   /* default value for the mach bound: 0.4 kg/s
    * maximal mach bound: 0.8
    * minimal mach bound: 0.2
    */
   SCIP_CALL( SCIPaddRealParam(scip, "cut_offset", "used to relax cuts of the form ax <= b by adding cut_offset to b",
         &(conshdlrdata->offset), FALSE, 1e-3, 1e-5, 1e-1, NULL, NULL) );
   /* default value for the cut_offset: 1e-3 bar
    * maximal value: 0.1
    * minimal value: 1e-5
    * @note: Should be smaller than pressure_feastol!
    */

   /* add parameter for flowTightening */
   SCIP_CALL( SCIPaddBoolParam(scip, "with_flowTightening", "if propagation based flowTightening should be used",
         &(conshdlrdata->with_flowTightening), FALSE, DEFAULT_WITH_FLOWTIGHTENING, NULL, NULL) );

   /* add parameters for obbt */
   SCIP_CALL( SCIPaddBoolParam(scip, "with_obbt", "if OBBT should be used",
         &(conshdlrdata->with_obbt), FALSE, DEFAULT_WITH_OBBT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "obbt_mingap", "only perform obbt for flow variable if bounds differ by at least obbt_mingap",
         &(conshdlrdata->obbt_mingap), FALSE, DEFAULT_OBBT_MINGAP, 1.0, SCIPinfinity(scip), NULL, NULL) );

   /* add parameter for NLP expressions */
   SCIP_CALL( SCIPaddBoolParam(scip, "nlp_representation", "whether the pipes should be represented in the NLP",
         &(conshdlrdata->nlp_representation), FALSE, DEFAULT_NLP_REPRESENTATION, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlp_convex", "whether the pipes are represented by a convex relaxation",
         &(conshdlrdata->nlp_convex), FALSE, DEFAULT_NLP_CONVEX, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a PipeODE constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsPipeODE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             p_in_var,           /**< variable for incoming pressure  */
   SCIP_VAR*             p_out_var,          /**< variable for leaving pressure  */
   SCIP_VAR*             q_var,              /**< variable for flow on pipe  */
   SCIP_VAR*             posFlowBinvar,      /**< binary variable if flow is nonnegative */
   SCIP_VAR*             negFlowBinvar,      /**< binary variable if flow is nonpositive */
   int                   N,                  /**< number of discretization steps  */
   SCIP_Real             A,                  /**< cross sectional area of pipe */
   SCIP_Real             L,                  /**< length of pipe */
   SCIP_Real             D,                  /**< diameter of  pipe */
   SCIP_Real             lambda,             /**< friction coefficient */
   SCIP_Real             c,                  /**< speed of sound */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generaion)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsPipeODE() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the PipeODE constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("PipeOde constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->p_in_var = p_in_var;
   consdata->p_out_var = p_out_var;
   consdata->q_var = q_var;
   consdata->posFlowBinvar = posFlowBinvar;
   consdata->negFlowBinvar = negFlowBinvar;
   consdata->N = N;
   consdata->A = A;
   consdata->Asq = A * A;
   consdata->L = L;
   consdata->D = D;
   consdata->lambda = lambda;
   consdata->c = c;
   consdata->csq = c * c;
   consdata->lcsqd = lambda * consdata->csq / D;

#ifdef STEPSIZE_DEBUG
   SCIPinfoMessage(scip, NULL, "Initial discretization on pipe <%s>:\t #%d \t %f m\n", name, N, L/N);
#endif

   consdata->nr_cave = 0;
   consdata->nr_vex = 0;
   consdata->feasible = FALSE;
   consdata->violation = 0.0;
   consdata->type_violation = -1;
   consdata->propvarind = 63;                /* if constraint is new try to propagate all bounds */

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a PipeODE constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicPipeODE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             p_in_var,           /**< variable for incoming pressure  */
   SCIP_VAR*             p_out_var,          /**< variable for leaving pressure  */
   SCIP_VAR*             q_var,              /**< variable for flow on pipe  */
   SCIP_VAR*             posFlowBinvar,      /**< binary variable if flow is nonnegative */
   SCIP_VAR*             negFlowBinvar,      /**< binary variable if flow is nonpositive */
   int                   N,                  /**< number of discretization steps  */
   SCIP_Real             A,                  /**< cross sectional area of pipe */
   SCIP_Real             L,                  /**< length of pipe */
   SCIP_Real             D,                  /**< diameter of  pipe */
   SCIP_Real             lambda,             /**< friction coefficient */
   SCIP_Real             c                   /**< speed of sound */
   )
{
   SCIP_CALL( SCIPcreateConsPipeODE(scip, cons, name, p_in_var, p_out_var, q_var, posFlowBinvar, negFlowBinvar, N, A, L, D, lambda, c, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** returns the number of gradient cuts */
int SCIPconsPipeODEGetNGradientCuts(
   SCIP_CONSHDLR*        conshdlr            /**< pointer to hold the created constraint */
   )
{
   assert( conshdlr != NULL );
   assert( SCIPconshdlrGetData(conshdlr) != NULL );

   return ( SCIPconshdlrGetData(conshdlr)->ngradientcuts );
}
