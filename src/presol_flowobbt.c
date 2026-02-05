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

/**@file    presol_flowobbt.c
 * @ingroup DEFPLUGINS_PRESOL
 * @brief   OBBT for flow bounds
 * @author  Marc Pfetsch
 *
 * OBBT within SCIP seems to be fruitless, so there is a specialized OBBT for gas problems here. If the ODE constraint
 * handler is applied, a specific OBBT is used. Thus, this presolver is only useful for the algebraic model.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "presol_flowobbt.h"
#include "scip/scipdefplugins.h"       /* needed for the secondary SCIP instance */
#include "elements_gas.h"

#define PRESOL_NAME             "flowobbt"
#define PRESOL_DESC             "apply OBBT for tightening flow bounds"
#define PRESOL_PRIORITY                -100000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS                0      /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING                   SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define OBBT_MINGAP                     10.0  /**< only perform OBBT for flow variable if bounds differ by at least OBBT_MINGAP */
#define OBBT_REMOVEBINVAR               TRUE  /**< remove all constraints that include binary variables */

/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   SCIP_Bool             ran;                /**< whether we ran already */
};


/*
 * Local methods
 */

#if OBBT_REMOVEBINVAR

/** remove all nonlinear constraints and constraints that include binary variables */
static
SCIP_RETCODE simplifyBinary(
   SCIP*                 scip,               /**< original scip */
   SCIP*                 subscip,            /**< subscip */
   SCIP_HASHMAP*         varmapfw,           /**< hashmap between original and sub-scip */
   SCIP_PRESOLDATA*      presoldata          /**< presolving data */
   )
{
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_Bool success;
   int nconss;
   int nvars;
   int i;

   assert( subscip != NULL );
   assert( varmapfw != NULL );
   assert( presoldata != NULL );

   /* get subscip data */
   conss = SCIPgetConss(subscip);
   nconss = SCIPgetNConss(subscip);
   nvars = SCIPgetNVars(subscip);

   SCIP_CALL( SCIPallocBufferArray(subscip, &vars, nvars) );

   for (i = nconss - 1; i >= 0; --i)
   {
      SCIP_Bool removecons = FALSE;
      int nconsvars;
      int j;

      if ( strcmp("nonlinear", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 )
      {
         SCIP_CALL( SCIPdelCons(subscip, conss[i]) );
      }
      else
      {
         SCIP_CALL( SCIPgetConsNVars(subscip, conss[i], &nconsvars, &success) );

         if ( success )
         {
            SCIP_CALL( SCIPgetConsVars(subscip, conss[i], vars, nvars, &success) );
            assert( success );

            for (j = 0; j < nconsvars; ++j)
            {
               if ( SCIPvarIsBinary(vars[j]) )
               {
                  removecons = TRUE;
                  break;
               }
            }

            if ( removecons )
            {
               SCIP_CALL( SCIPdelCons(subscip, conss[i]) );
            }
         }
      }
   }

   /* Do not delete variables, since they will be removed in presolving and removal does not seem to work, because of
    * the negated variables which do not seem to be deleted. */

   SCIPfreeBufferArray(subscip, &vars);

   return SCIP_OKAY;
}

#else

/** remove nonlinear constraints */
static
SCIP_RETCODE simplifyNLTV(
   SCIP*                 scip,               /**< original scip */
   SCIP*                 subscip,            /**< subscip */
   SCIP_HASHMAP*         varmapfw,           /**< hashmap between original and sub-scip */
   SCIP_PRESOLDATA*      presoldata          /**< presolving data */
   )
{
   SCIP_CONS** conss;
   int nconss;
   int i;

   assert( subscip != NULL );
   assert( varmapfw != NULL );
   assert( presoldata != NULL );

   /* get subscip constraint data */
   conss = SCIPgetConss(subscip);
   nconss = SCIPgetNConss(subscip);

   for (i = nconss-1; i >= 0; --i)
   {
      if ( strcmp("nonlinear", SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])) ) == 0 )
      {
         SCIP_CALL( SCIPdelCons(subscip, conss[i]) );
      }
   }

   return SCIP_OKAY;
}
#endif


/*
 * Callback methods of presolver
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeFlowOBBT)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecFlowOBBT)
{  /*lint --e{715}*/
   char probname[SCIP_MAXSTRLEN];
   SCIP_PRESOLDATA* presoldata;
   SCIP_PROBDATA* origdata;
   SCIP_HASHMAP* varmapfw;
   SCIP* subscip;
   SCIP_VAR** subvars;
   SCIP_VAR* flowvar;
   SCIP_VAR* subflowvar;
   SCIP_Bool success;
   SCIP_Real dualbound;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real lb_improvement = 0.0;
   SCIP_Real ub_improvement = 0.0;
   SCIP_Real timelimit = 3600.0;
   SCIP_RETCODE retcode;
   SCIP_STATUS status;
   int noldfixedvars;
   int noldchgbds;
   int nsubvars;
   int napplied = 0;
   int i;

   assert( scip != NULL );
   assert( presol != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* do not run for sub-SCIPs */
   if ( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   /* get presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert( presoldata != NULL );

   /* exit if we already ran */
   if ( presoldata->ran )
      return SCIP_OKAY;

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );
   if ( ! success )
      return SCIP_OKAY;

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), SCIPgetNVars(scip)) );

   /* create a problem copy as sub SCIP */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "obbt#subscip#");
   SCIP_CALL( SCIPcopy(scip, subscip, varmapfw, NULL, probname, TRUE, FALSE, TRUE, FALSE, &success ) );

   if ( ! success )
   {
      SCIPerrorMessage("Was not able to copy SCIP during OBBT.\n");
      SCIPhashmapFree(&varmapfw);
      SCIP_CALL( SCIPfree(&subscip) );
      return SCIP_OKAY; /* do not return with an error */
   }

   /* get variables of the subscip */
   subvars = SCIPgetVars(subscip);
   nsubvars = SCIPgetNVars(subscip);

   /* change objective value of subscip to 0.0 */
   for (i = 0; i < nsubvars; i++)
   {
      SCIP_CALL( SCIPchgVarObj(subscip, subvars[i], 0.0) );
   }

   /* get problem data */
   origdata = SCIPgetProbData(scip);
   noldfixedvars = *nfixedvars;
   noldchgbds = *nchgbds;

   /* simplify problem */
#if OBBT_REMOVEBINVAR
   SCIP_CALL( simplifyBinary(scip, subscip, varmapfw, presoldata) );
#else
   SCIP_CALL( simplifyNLTV(scip, subscip, varmapfw, presoldata) );
#endif

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) starting OBBT ...\n", SCIPgetSolvingTime(scip));

#ifdef SCIP_MORE_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 100000000) );

   SCIP_CALL( SCIPpresolve(subscip) );
   SCIP_CALL( SCIPwriteTransProblem(subscip, "flowobbt.cip", "cip", FALSE) );
   SCIP_CALL( SCIPfreeTransform(subscip) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* copy time and memory limits */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );

   /* update time limit */
   SCIP_CALL( SCIPgetRealParam(subscip, "limits/time", &timelimit) );
   if ( timelimit > 600.0 )
   {
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", 600.0) );
   }

   /* loop over all arcs and try to minimize and maximize the flow variables */
   for (i = origdata->network->numarcs - 1; i >= 0; --i)
   {
      flowvar = SCIPvarGetTransVar(origdata->flowvars[i]);
      subflowvar = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, flowvar);

      /* possibly some variables have been deleted */
      if ( subflowvar == NULL )
         continue;

      /* cannot change bounds of (multi-)aggregated variables, irgnore already fixed variables */
      assert( SCIPvarGetStatus(flowvar) != SCIP_VARSTATUS_ORIGINAL );
      assert( SCIPvarGetStatus(flowvar) != SCIP_VARSTATUS_NEGATED );
      if ( SCIPvarGetStatus(flowvar) == SCIP_VARSTATUS_MULTAGGR || SCIPvarGetStatus(flowvar) == SCIP_VARSTATUS_AGGREGATED || SCIPvarGetStatus(flowvar) == SCIP_VARSTATUS_FIXED )
         continue;

      /* only try to compute new bounds if the current ones are not too close */
      lb = SCIPcomputeVarLbGlobal(scip, flowvar);
      ub = SCIPcomputeVarUbGlobal(scip, flowvar);

      /* do not perform obbt, if bounds are close; obbt_mingap can be changed with the settings */
      if ( ub - lb < OBBT_MINGAP )
         continue;

      /* count variables obbt is performed on */
      ++napplied;

      /* First try to minimize the flow varible. One could also try to maximize first, there is no reason for this
       * order. */
      SCIP_CALL( SCIPchgVarObj(subscip, subflowvar, 1.0) );

#if 0
      SCIP_CALL( SCIPpresolve(subscip) );
      SCIP_CALL( SCIPwriteTransProblem(subscip, "flowobbttrans.cip", "cip", FALSE) );
#endif

      /* solve subscip */
      retcode = SCIPsolve(subscip);
      if ( retcode != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving subproblem in OBBT; sub-SCIP terminated with code <%d>.\n", retcode);
         SCIPABORT();

         break;
      }

      /* check solving status */
      status = SCIPgetStatus(subscip);
      if ( status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_INFEASIBLE || status == SCIP_STATUS_TIMELIMIT )
         break;

      /* get dual bound of solution process -> provides variable bound in original problem */
      dualbound = SCIPgetDualbound(subscip);

      /* free transformed problem such that we can optimize again */
      SCIP_CALL( SCIPfreeTransform(subscip) );

      /* have to change bounds of subflowvar after freeing the transformed problem */
      if ( status == SCIP_STATUS_OPTIMAL )  /* other status could also lead to a valid bound */
      {
         assert( SCIPisFeasGE(scip, dualbound, lb) );
         assert( SCIPisFeasLE(scip, dualbound, ub) );

         /* if the bound has been improved */
         if ( SCIPisGT(scip, dualbound, lb) )
         {
            /* if lower bound is equal to upper bound */
            if ( SCIPisFeasEQ(scip, dualbound, ub) )
            {
               dualbound = ub;
               ++(*nfixedvars);
            }
            else
               ++(*nchgbds);

            /* change bounds */
            SCIP_CALL( SCIPchgVarLbGlobal(scip, flowvar, dualbound) );
            SCIP_CALL( SCIPchgVarLbGlobal(subscip, subflowvar, dualbound) );
            lb_improvement += dualbound - lb;
            lb = dualbound;
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

      /* solve subscip */
      retcode = SCIPsolve(subscip);
      if ( retcode != SCIP_OKAY )
      {
         SCIPwarningMessage(scip, "Error while solving subproblem in OBBT; sub-SCIP terminated with code <%d>.\n", retcode);
         SCIPABORT();

         break;
      }

      /* check solving status */
      status = SCIPgetStatus(subscip);
      if ( status == SCIP_STATUS_USERINTERRUPT || status == SCIP_STATUS_INFEASIBLE || status == SCIP_STATUS_TIMELIMIT )
         break;

      /* get dual bound of solution process -> provides variable bound in original problem */
      dualbound = - SCIPgetDualbound(subscip);

      /* free transformed problem such that we can optimize again */
      SCIP_CALL( SCIPfreeTransform(subscip) );

      /* have to change bounds of subflowvar after freeing the transformed problem */
      if ( status == SCIP_STATUS_OPTIMAL )  /* other status could also lead to a valid bound */
      {
         assert( SCIPisFeasGE(scip, dualbound, lb) );
         assert( SCIPisFeasLE(scip, dualbound, ub) );

         /* if the bound has been improved */
         if ( SCIPisLT(scip, dualbound, ub) )
         {
            /* if upper bound bound is equal to lower bound */
            if ( SCIPisFeasEQ(scip, dualbound, lb) )
            {
               dualbound = lb;
               ++(*nfixedvars);
            }
            else
               ++(*nchgbds);

            /* change bounds */
            SCIP_CALL( SCIPchgVarUbGlobal(scip, flowvar, dualbound) );
            SCIP_CALL( SCIPchgVarUbGlobal(subscip, subflowvar, dualbound) );
            ub_improvement += ub - dualbound;
         }
      }

      SCIP_CALL( SCIPchgVarObj(subscip, subflowvar, 0.0) );

      /* go to next variable */
   }

   /* some output */
   if ( status != SCIP_STATUS_INFEASIBLE && status != SCIP_STATUS_USERINTERRUPT )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) OBBT on %d variables produced %d bound changes and %d fixed variables.\n",
         SCIPgetSolvingTime(scip), napplied, *nchgbds - noldchgbds, *nfixedvars - noldfixedvars);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "           Improved lower bounds by: %10.2f\n", lb_improvement);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "           Improved upper bounds by: %10.2f\n", ub_improvement);
   }

   if ( (*nchgbds - noldchgbds) + (*nfixedvars - noldfixedvars) > 0 )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* finally free subscipdata */
   SCIP_CALL( SCIPfree(&subscip) );

   presoldata->ran = TRUE;

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the flow OBBT presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolFlowOBBT(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );
   presoldata->ran = FALSE;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecFlowOBBT, presoldata) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeFlowOBBT) );

   return SCIP_OKAY;
}
