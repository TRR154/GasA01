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

/**@file   unitConversion_gas.h
 * @brief  conversion of units used in stationary gas transport problems
 * @author Imke Joormann
 */

#ifndef __UNITCONVERSION_GAS__
#define __UNITCONVERSION_GAS__

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

SCIP_EXPORT
SCIP_RETCODE change_unit_meter(
   const char*           unit,
   SCIP_Real*            value
   );

SCIP_EXPORT
SCIP_RETCODE change_unit_bar(
   const char*           unit,
   SCIP_Real*            value
   );

SCIP_EXPORT
SCIP_RETCODE change_flow_unit(
   const char*           unit,
   SCIP_Real*            value,
   SCIP_Real            normDensity
   );

SCIP_EXPORT
SCIP_RETCODE change_density_unit(
   const char*           unit,
   SCIP_Real*               value
   );

SCIP_EXPORT
SCIP_RETCODE change_temperature_unit(
   const char*           unit,
   SCIP_Real*               value
   );

SCIP_EXPORT
SCIP_RETCODE change_molarmass_unit(
   const char*           unit,
   SCIP_Real*            value
   );

SCIP_EXPORT
SCIP_RETCODE change_enthalpy_unit(
   const char*           unit,
   SCIP_Real*            value
   );

SCIP_EXPORT
SCIP_RETCODE check_constant_unit(
   const char*           unit,
   SCIP_Real*            value
   );

SCIP_EXPORT
SCIP_RETCODE change_gasconstant_unit(
   const char*           unit,
   SCIP_Real*            value
   );

#ifdef __cplusplus
}
#endif

#endif
