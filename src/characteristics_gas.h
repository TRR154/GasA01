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

/**@file   characteristics_gas.h
 * @brief  handling of data needed for solving stationary gas transport problems
 * @author Alexandra Stern
 * @author Imke Joormann
 */

#ifndef __CHARACTERISTICS_GAS__
#define __CHARACTERISTICS_GAS__

#include <scip/scip.h>
#include <scip/misc.h>
#include "elements_gas.h"

#ifdef __cplusplus
extern "C" {
#endif

/* general constants */
extern SCIP_Real pi;                         /**< value of Pi */
extern SCIP_Real universalGasConstant;       /**< universal gas constant */
extern SCIP_Real p_norm;                     /**< constant value for pressure, unit: bar */
extern SCIP_Real T_0;                        /**< constant value for temperature, unit: K */


/** compute the constant in the Weymouth constraint */
SCIP_EXPORT
SCIP_Real computeWeymouthConstant(
   GAS_Network*          network,            /**< network data structure */
   SCIP_Real             length,             /**< length of pipe */
   SCIP_Real             diameter,           /**< diameter of pipe */
   SCIP_Real             roughness,          /**< roughness of pipe */
   SCIP_Real             meanpressure,       /**< mean pressure value */
   SCIP_Bool             papay               /**< if formula of papay is used */
   );

/** compute the constant in the MIXING Weymouth constraint */
SCIP_EXPORT
SCIP_Real computeWeymouthConstantMixing(
   GAS_Network*          network,            /**< network data structure */
   SCIP_Real             length,             /**< length of pipe */
   SCIP_Real             diameter,           /**< diameter of pipe */
   SCIP_Real             roughness,          /**< roughness of pipe */
   SCIP_Real             meanpressure,       /**< mean pressure value */
   SCIP_Bool             papay               /**< if formula of papay is used */
   );

/** compute mean compressibility factor using mean pressure and constant gas temperature */
SCIP_EXPORT
SCIP_Real MeanCompressibilityFactor(
   GAS_Network*          network,            /**< gas network data structure */
   SCIP_Real             meanpressure,       /**< mean pressure value */
   SCIP_Bool             papay               /**< if formula of papay is used */
   );

/** compute the speed of sound */
SCIP_EXPORT
SCIP_Real computeSpeedOfSound(
   SCIP_Real             T_m,                /**< gas temperature */
   SCIP_Real             z_m,                /**< mean compressibility factor */
   SCIP_Real             molarMass           /**< molar mass of the gas */
   );

/** compute the MIXED speed of sound */
SCIP_EXPORT
SCIP_Real computeLambda(
   SCIP_Real             T_m,                /**< gas temperature */
   SCIP_Real             z_m,                /**< mean compressibility factor */
   SCIP_Real             mixingRatio         /**< mixing ratio of the gases (molar fractions) */
   );

/** calculate mol % from mass % */
SCIP_EXPORT
SCIP_Real massToMolePercent(
   SCIP_Real             mixingRatio,        /**< the mixing ratio in Mass % */
   SCIP_Real             massFlow            /**< value of the mass flow in kg/s */
);

/** compute Nikuradze value */
SCIP_EXPORT
SCIP_Real NikuradzesEquation (
   SCIP_Real             diameter,           /**< diameter of pipe */
   SCIP_Real             roughness           /**< roughness of pipe */
   );

/** mean pressure to use for compressibility factor */
SCIP_EXPORT
SCIP_Real ComputeMeanPressure(
   SCIP_Real             lowerpressureV,     /**< lower pressure of GAS_Node V */
   SCIP_Real             upperpressureV,     /**< upper pressure of GAS_Node V */
   SCIP_Real             lowerpressureW,     /**< lower pressure of GAS_Node W */
   SCIP_Real             upperpressureW      /**< upper pressure of GAS_Node V */
   );

/** AGA approximation for compressibility factor */
SCIP_EXPORT
SCIP_Real AGA(
   SCIP_Real             T_c,                /**< pseudocriticalTemperature */
   SCIP_Real             T_m,                /**< gasTemperature */
   SCIP_Real             p_c,                /**< pseudocriticalpressure */
   SCIP_Real             p_m                 /**< mean pressure */
   );

/** formula of Papay for compressibility factor */
SCIP_EXPORT
SCIP_Real Papay(
   SCIP_Real             T_c,                /**< pseudocriticalTemperature */
   SCIP_Real             T_m,                /**< gasTemperature */
   SCIP_Real             p_c,                /**< pseudocriticalpressure */
   SCIP_Real             p_m                 /**< mean pressure */
   );

/** compute constant for nonlinear resistor */
SCIP_EXPORT
SCIP_Real ConstantNonlinearResistor(
   GAS_Network*          network,            /**< network data structure */
   SCIP_Real             meanpressure,       /**< mean pressure value */
   SCIP_Real             dragfactor,         /**< dragfactor of pipe */
   SCIP_Real             diameter,           /**< diameter of pipe */
   SCIP_Bool             papay               /**< if formula of papay is used */
   );

/** calculate slope
 *
 * Use length and heightdiff to calculate the adjacent side. Then divide heightdiff by adjacent to get the slope.
 */
SCIP_EXPORT
SCIP_Real calculateSlope(
   SCIP_Real             heightdiff,         /**< height difference of innode to outnode */
   SCIP_Real             length              /**< pipe length */
   );

#ifdef __cplusplus
}
#endif

#endif
