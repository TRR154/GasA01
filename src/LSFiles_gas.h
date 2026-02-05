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

/**@file   LSFiles_gas.h
 * @brief  Creates and reads ls-file containing information about the network and the solution values of the variables
 * @author Alexandra Stern
 * @author Oliver Habeck
 */

#include <scip/misc.h>
#include <scip/scip.h>

#ifndef __LSFILES_GAS__
#define __LSFILES_GAS__

#include <scip/scip.h>
#include <scip/misc.h>


#ifdef __cplusplus
extern "C" {
#endif

/** write LS file */
SCIP_EXPORT
SCIP_RETCODE writeLSF(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_SOL*             sol                 /**< solution to be written */
   );

/** reads LS file and fixes variables to the given values */
SCIP_EXPORT
SCIP_RETCODE readLSF(
   SCIP*                 scip,               /**< SCIP instance */
   const char*           pathToLsf           /**< path to LSF file */
   );

#ifdef __cplusplus
}
#endif

#endif
