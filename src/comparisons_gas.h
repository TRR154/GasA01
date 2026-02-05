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

/**@file   comparisons_gas.h
 * @brief  functions to compare pressure and flow values with specific feasibility tolerances
 * @author Oliver Habeck
 */

#ifndef __COMPARISONS_GAS__
#define __COMPARISONS_GAS__

#include <scip/scip.h>
#include <scip/misc.h>

#ifdef __cplusplus
extern "C" {
#endif

/** checks if pressure values are in range of eps,
 *  i.e., ABS(pval1 - pval2) <= eps
 */
SCIP_EXPORT
SCIP_Bool pressureIsEQ(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   );

/** checks if pval1 is eps less than pval2,
 *  i.e., pval1 - pval2 < -eps
 */
SCIP_EXPORT
SCIP_Bool pressureIsLT(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   );

/** checks if pval1 is less or equal to pval2 with tolerance eps,
 *  i.e., pval1 - pval2 <= eps
 */
SCIP_EXPORT
SCIP_Bool pressureIsLE(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   );

/** checks if pval1 is eps greater than pval2,
 *  i.e., pval1 - pval2 > eps
 */
SCIP_EXPORT
SCIP_Bool pressureIsGT(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   );

/** checks if pval1 is greater or equal to pval2 with tolerance eps,
 *  i.e., pval1 - pval2 >= -eps
 */
SCIP_EXPORT
SCIP_Bool pressureIsGE(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   );

/** checks if flow values are in range of eps,
 *  i.e., REALABS(fval1 - fval2) <= eps
 */
SCIP_EXPORT
SCIP_Bool flowIsEQ(
   SCIP_Real             fval1,              /**< first flow value to be compared */
   SCIP_Real             fval2,              /**< second flow value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   );

/** checks if flow value is zero,
 *  i.e., REALABS(fval) <= eps
 */
SCIP_EXPORT
SCIP_Bool flowIsZero(
   SCIP_Real             fval,               /**< flow value */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   );

/** checks if flow value is positive,
 *  i.e., fval > eps
 */
SCIP_EXPORT
SCIP_Bool flowIsPositive(
   SCIP_Real             fval,               /**< flow value */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   );

/** checks if flow value is negative,
 *  i.e., fval < - eps
 */
SCIP_EXPORT
SCIP_Bool flowIsNegative(
   SCIP_Real             fval,               /**< flow value */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   );


#ifdef __cplusplus
}
#endif

#endif
