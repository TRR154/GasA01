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

/**@file   generateModel_gas.h
 * @brief  generate constraints for describing stationary gas transport networks
 * @author Imke Joormann
 */

#ifndef __GENERATEMODEL_GAS__
#define __GENERATEMODEL_GAS__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif



/** create stationary gas transport problem model */
SCIP_EXPORT
SCIP_RETCODE GASgenerateModel(
   SCIP*                 scip                /**< SCIP data structure */
   );


#ifdef __cplusplus
}
#endif

#endif

