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

/**@file   comparisons_gas.c
 * @brief  functions to compare pressure and flow values with specific feasibility tolerances
 * @author Oliver Habeck
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "comparisons_gas.h"
#include <scip/scip.h>
#include <scip/misc.h>


/** checks if pressure values are in range of eps,
 *  i.e., REALABS(pval1 - pval2) <= eps
 */
SCIP_Bool pressureIsEQ(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   )
{
   return EPSEQ(pval1, pval2, eps);
}

/** checks if pval1 is eps less than pval2,
 *  i.e., pval1 - pval2 < -eps
 */
SCIP_Bool pressureIsLT(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   )
{
   return EPSLT(pval1, pval2, eps);
}

/** checks if pval1 is less or equal to pval2 with tolerance eps,
 *  i.e., pval1 - pval2 <= eps
 */
SCIP_Bool pressureIsLE(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   )
{
   return EPSLE(pval1, pval2, eps);
}


/** checks if pval1 is eps greater than pval2,
 *  i.e., pval1 - pval2 > eps
 */
SCIP_Bool pressureIsGT(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   )
{
   return EPSGT(pval1, pval2, eps);
}


/** checks if pval1 is greater or equal to pval2 with tolerance eps,
 *  i.e., pval1 - pval2 >= -eps
 */
SCIP_Bool pressureIsGE(
   SCIP_Real             pval1,              /**< first pressure value to be compared */
   SCIP_Real             pval2,              /**< second pressure value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   )
{
   return EPSGE(pval1, pval2, eps);
}


/** checks if flow values are in range of eps,
 *  i.e., REALABS(fval1 - fval2) <= eps
 */
SCIP_Bool flowIsEQ(
   SCIP_Real             fval1,              /**< first flow value to be compared */
   SCIP_Real             fval2,              /**< second flow value to be compared */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   )
{
   return EPSEQ(fval1, fval2, eps);
}


/** checks if flow value is zero,
 *  i.e., REALABS(fval) <= eps
 */
SCIP_Bool flowIsZero(
   SCIP_Real             fval,               /**< flow value */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   )
{
   return EPSZ(fval, eps);
}


/** checks if flow value is positive,
 *  i.e., fval > eps
 */
SCIP_Bool flowIsPositive(
   SCIP_Real             fval,               /**< flow value */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   )
{
   return EPSP(fval, eps);
}


/** checks if flow value is negative,
 *  i.e., fval < - eps
 */
SCIP_Bool flowIsNegative(
   SCIP_Real             fval,               /**< flow value */
   SCIP_Real             eps                 /**< feasibility tolerance for comparison */
   )
{
   return EPSN(fval, eps);
}
