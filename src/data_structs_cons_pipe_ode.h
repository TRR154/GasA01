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

/**@file   data_structs_cons_pipe_ode.h
 * @brief  provides typedef for conshdlrdata and consdata
 * @author Oliver Habeck
 */

#include <scip/type_cons.h>

/*
 * Data structures
 */

#define PIPEODE_NONE     0              /**< no bound has changed */
#define PIPEODE_PIN_LB   0x00000001u    /**< the lower bound of p_in has changed */
#define PIPEODE_PIN_UB   0x00000002u    /**< the upper bound of p_in has changed */
#define PIPEODE_POUT_LB  0x00000004u    /**< the lower bound of p_out has changed */
#define PIPEODE_POUT_UB  0x00000008u    /**< the upper bound of p_out has changed */
#define PIPEODE_Q_LB     0x00000010u    /**< the lower bound of q has changed */
#define PIPEODE_Q_UB     0x00000020u    /**< the upper bound of q has changed */
#define PIPEODE_POS_LB   0x00000040u    /**< the lower bound of positive flow binvar has changed */
#define PIPEODE_POS_UB   0x00000080u    /**< the upper bound of positive flow binvar has changed */
#define PIPEODE_NEG_LB   0x00000100u    /**< the lower bound of negative flow binvar has changed */
#define PIPEODE_NEG_UB   0x00000200u    /**< the upper bound of negative flow binvar has changed */

typedef unsigned int PIPEODEVARIND;     /**< indicator for changes in variable bounds in pipe ODE constraint */

/** constraint data for pipe_ode constraints */
struct SCIP_ConsData
{
   SCIP_VAR*             p_in_var;       /**< variable for pressure entering pipe */
   SCIP_VAR*             p_out_var;      /**< variable for pressure leaving pipe */
   SCIP_VAR*             q_var;          /**< variable for flow on pipe */
   SCIP_VAR*             posFlowBinvar;  /**< binary variable if flow is nonnegative */
   SCIP_VAR*             negFlowBinvar;  /**< binary variable if flow is nonpositive */
   int                   N;              /**< number of discretization steps  */
   SCIP_Real             A;              /**< cross sectional area of pipe */
   SCIP_Real             Asq;            /**< squared cross sectional area of pipe */
   SCIP_Real             L;              /**< length of pipe */
   SCIP_Real             D;              /**< diameter of  pipe */
   SCIP_Real             lambda;         /**< friction coefficient */
   SCIP_Real             c;              /**< speed of sound */
   SCIP_Real             csq;            /**< speed of sound squared */
   SCIP_Real             lcsqd;          /**< lambda*c^2/D */
   int                   nr_cave;        /**< number of all linear overestimating constraints added */
   int                   nr_vex;         /**< number of all gradient cuts added */
   SCIP_Bool             feasible;       /**< boolean whether or not a solution was feasible when checked last */
   SCIP_Real             violation;      /**< absolute value of the violation */
   int                   type_violation; /**< type of the violation */
   PIPEODEVARIND         propvarind;     /**< indicator which variable bound has changed (if any), used for propagation */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   int                   ngradientcuts;      /**< number of gradient cuts */
   int                   machbound;          /**< numerator used for bound on mach number: v/c = machbound/5 */
   SCIP_Real             pressure_feastol;   /**< tolerance used for comparisons of pressure values */
   SCIP_Real             flow_feastol;       /**< tolerance used for comparisons of flow values */
   SCIP_Real             offset;             /**< used to relax cuts of the form ax <= b by adding cut_offset to b */
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for bound change events */
   SCIP_Bool             with_flowTightening;/**< if propagation based flowTightening should be used */
   SCIP_Bool             with_obbt;          /**< if obbt should be used */
   SCIP_Real             obbt_mingap;        /**< only perform obbt for variable if bounds differ by at least mingap */
   SCIP_Bool             nlp_representation; /**< whether the pipes should be represented in the NLP */
   SCIP_Bool             nlp_convex;         /**< whether the pipes are represented by a convex relaxation */
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
   SCIP_EXPRHDLR*        exprhdlr;           /**< expression handler for pipes */
#endif
};
