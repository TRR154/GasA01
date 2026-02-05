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

/**@file cons_pipe_ode.h
 * @ingroup CONSHDLRS
 * @brief  constraint handler for PipeODE constraints
 * @author Oliver Habeck
 * @author Alexandra Stern
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_PIPE_ODE_H__
#define __SCIP_CONS_PIPE_ODE_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for PipeOde constraints and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrPipeODE(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates and captures a PipeODE constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
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
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
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
   );

/** creates and captures a PipeODE constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsBasicPipeODE(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_VAR*             p_in_var,           /**< variable for incoming pressure  */
   SCIP_VAR*             p_out_var,          /**< variable for leaving pressure  */
   SCIP_VAR*             q_var,              /**< variable for flow on pipe  */
   SCIP_VAR*             posFlowBinvar,      /**< binary variable if flow is nonnegative */
   SCIP_VAR*             negFlowBinvar,      /**< binary variable if flow is nonpositive */
   int                   N,                  /**< number of discretization steps */
   SCIP_Real             A,                  /**< cross sectional area of pipe */
   SCIP_Real             L,                  /**< length of pipe */
   SCIP_Real             D,                  /**< diameter of  pipe */
   SCIP_Real             lambda,             /**< friction coefficient */
   SCIP_Real             c                   /**< speed of sound */
   );

/** returns the number of gradient cuts */
SCIP_EXPORT
int SCIPconsPipeODEGetNGradientCuts(
   SCIP_CONSHDLR*        conshdlr            /**< pointer to hold the created constraint */
   );

#ifdef __cplusplus
}
#endif

#endif
