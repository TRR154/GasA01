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

/**@file   probdata_gas.h
 * @brief  handling of data needed for solving stationary gas transport problems
 * @author Marc Pfetsch
 */

#ifndef __PROBDATA_GAS__
#define __PROBDATA_GAS__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/* forward declaration */
struct GAS_network;

/** initialize the probdata structure */
SCIP_EXPORT
SCIP_RETCODE GASinitProb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< problem data structure */
   );

/** create stationary gas transport problem instance */
SCIP_EXPORT
SCIP_RETCODE GAScreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   const char*           netfilename,        /**< name of network file to read */
   const char*           scnfilename,        /**< name of scenario file to read */
   const char*           csfilename          /**< name of compressor file to read */
   );

/** create problem instance with existing network and determine neighbors */
SCIP_EXPORT
SCIP_RETCODE GAScreateProbNetwork(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   struct GAS_network*   network             /**< network */
   );

/** write network in xml file */
SCIP_EXPORT
SCIP_RETCODE GASwriteNetwork(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   const char*           filename,           /**< name of network file to write */
   SCIP_Bool             withframework       /**< use key "framework"? */
   );

#ifdef __cplusplus
}
#endif

#endif
