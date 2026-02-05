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

/**@file   ode_functions_gas.h
 * @brief  functions used for the relaxation of the stationary, isothermal Euler equations
 * @author Alexandra Stern
 */

#ifndef __ODE_FUNCTIONS_GAS__
#define __ODE_FUNCTIONS_GAS__

#include <scip/scip.h>
#include <scip/misc.h>

#ifdef __cplusplus
extern "C" {
#endif

/** computes the modified euler scheme either in direction of flow (dir==1)
 *  or against (dir==-1).
 *
 *  Returns -1 if either the pressure falls below the
 *  critical bound 5/4 * c * q/A (dir == 1) or if something goes terribly wrong and
 *  the pressure decreases instead of increases when dir == -1.
 */
SCIP_EXPORT
SCIP_Real computeModifiedEulerScheme(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of the constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   int                   dir                 /**< direction of computation */
   );

/** computes the modified euler scheme against the direction of the flow
 *  and the gradient of the "incoming" pressure with respect to the flow
 *  and the "outgoing" pressure
 */
SCIP_EXPORT
SCIP_Real computeModifiedEulerSchemeWithGradients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   SCIP_Real*            dp,                 /**< pointer to the partial derivative dp with no unit */
   SCIP_Real*            dq                  /**< pointer to the partial derivative dq in  bar * (kg / s)^{-1} */
   );

/** compute the pressure bound via the trapezoidal scheme */
SCIP_EXPORT
SCIP_Real computeTrapezoidalScheme(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   int                   dir                 /**< direction of computation */
   );

/** computes RK4 solutions either in direction of flow (dir==1) or against (dir==-1). */
SCIP_EXPORT
SCIP_Real computeRK4(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of the constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   int                   dir                 /**< direction of computation */
   );

/** computes RK4 solutions against the direction of the flow and the gradient of the "incoming" pressure with respect to
 *  the flow and the "outgoing" pressure.
 */
SCIP_EXPORT
SCIP_Real computeRK4WithGradients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   SCIP_Real*            dp,                 /**< pointer to the partial derivative dp with no unit */
   SCIP_Real*            dq                  /**< pointer to the partial derivative dq in  bar * (kg / s)^{-1} */
   );

#ifdef __cplusplus
}
#endif

#endif
