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

/**@file   characteristics_gas.c
 * @brief  handling of data needed for creating constraints
 * @author Alexandra Stern
 */

#include "characteristics_gas.h"
#include <scip/misc.h>
#include <math.h>

/* #define SCIP_DEBUG */

/* general constants */
SCIP_Real pi = 3.14159265359;
SCIP_Real universalGasConstant = 8.3144621;  /**< unit: universal gas constant */
SCIP_Real p_norm = 1.01325;                  /**< constant value for pressure, unit: bar */
SCIP_Real T_0 = 273.15;                      /**< constant value for temperature, unit: K */


/** compute Nikuradze value */
SCIP_Real NikuradzesEquation(
   SCIP_Real             diameter,           /**< diameter of pipe */
   SCIP_Real             roughness           /**< roughness of pipe */
   )
{
   SCIP_Real lambda;

   lambda = pow(2 * log10(diameter / roughness) + 1.138, -2.0);

   return lambda;
}

/** compute mean compressibility factor using mean pressure and constant gas temperature */
SCIP_Real MeanCompressibilityFactor(
   GAS_Network*          network,            /**< gas network data structure */
   SCIP_Real             meanpressure,       /**< mean pressure value */
   SCIP_Bool             papay               /**< if formula of papay is used */
   )
{
   SCIP_Real z_m;

   if ( papay )
   {
      z_m = Papay(network->pseudocriticalTemperature, network->gasTemperature, network->pseudocriticalPressure, meanpressure);
   }
   else
   {
      z_m = AGA(network->pseudocriticalTemperature, network->gasTemperature, network->pseudocriticalPressure, meanpressure);
   }

   return z_m;
}

/** compute the speed of sound */
SCIP_Real computeSpeedOfSound(
   SCIP_Real             T_m,                /**< gas temperature */
   SCIP_Real             z_m,                /**< mean compressibility factor */
   SCIP_Real             molarMass           /**< molar mass of the gas */
   )
{
   SCIP_Real c;

   c = sqrt( universalGasConstant * T_m * z_m / molarMass );

   return c;
}

/** mean pressure to use for compressibility factor */
SCIP_Real ComputeMeanPressure(
   SCIP_Real             lowerpressureV,     /**< lower pressure of GAS_Node V */
   SCIP_Real             upperpressureV,     /**< upper pressure of GAS_Node V */
   SCIP_Real             lowerpressureW,     /**< lower pressure of GAS_Node W */
   SCIP_Real             upperpressureW      /**< upper pressure of GAS_Node W */
   )
{
   SCIP_Real p;

   p = 0.5 * (fmax(lowerpressureV, lowerpressureW) + fmin(upperpressureV, upperpressureW));

   return p;
}

/** AGA approximation for compressibility factor */
SCIP_Real AGA(
   SCIP_Real             T_c,                /**< pseudocriticalTemperature */
   SCIP_Real             T_m,                /**< gasTemperature */
   SCIP_Real             p_c,                /**< pseudocriticalpressure */
   SCIP_Real             p_m                 /**< pressure */
   )
{
   SCIP_Real z_m;

   z_m = 1 + 0.257 * (p_m / p_c) - 0.533 * (p_m / p_c) * (T_c / T_m);

   return z_m;
}

/** formula of Papay for compressibility factor */
SCIP_Real Papay(
   SCIP_Real             T_c,                /**< pseudocriticalTemperature */
   SCIP_Real             T_m,                /**< gasTemperature */
   SCIP_Real             p_c,                /**< pseudocriticalpressure */
   SCIP_Real             p_m                 /**< pressure */
   )
{
   SCIP_Real z_m;

   z_m = 1 - 3.52 * (p_m / p_c) * exp(-2.26 * T_m / T_c) + 0.247 * (p_m / p_c) * (p_m / p_c) * exp(-1.878 * T_m / T_c);

   return z_m;
}

/** compute the constant in the Weymouth constraint */
SCIP_Real computeWeymouthConstant(
   GAS_Network*          network,            /**< network data structure */
   SCIP_Real             length,             /**< length of pipe */
   SCIP_Real             diameter,           /**< diameter of pipe */
   SCIP_Real             roughness,          /**< roughness of pipe */
   SCIP_Real             meanpressure,       /**< mean pressure value */
   SCIP_Bool             papay               /**< if formula of papay is used */
   )
{
   SCIP_Real lambda;
   SCIP_Real z_m;
   SCIP_Real LAMBDA;

   /*  friction coefficient with Nikuradze */
   lambda = NikuradzesEquation(diameter, roughness);

   z_m = MeanCompressibilityFactor(network, meanpressure, papay);

   /* To ensure that units are right in the weymouth equality we have to divide with 10^10 such that the
    * unit of LAMBDA is bar^2 * (kg/s)^-2. */

   LAMBDA = (16.0 * lambda * length * z_m * network->gasTemperature * universalGasConstant) / (pow( pi, 2.0 ) * pow( diameter, 5.0 ) * network->molarMass1 * 1e10);

   return LAMBDA;
}

/** compute the constant in the MIXING Weymouth constraint */
SCIP_Real computeWeymouthConstantMixing(
   GAS_Network*          network,            /**< network data structure */
   SCIP_Real             length,             /**< length of pipe */
   SCIP_Real             diameter,           /**< diameter of pipe */
   SCIP_Real             roughness,          /**< roughness of pipe */
   SCIP_Real             meanpressure,       /**< mean pressure value */
   SCIP_Bool             papay               /**< if formula of papay is used */
   )
{
   SCIP_Real lambda;
   SCIP_Real z_m;
   SCIP_Real LAMBDA;

   /*  friction coefficient with Nikuradze */
   lambda = NikuradzesEquation( diameter, roughness ); /* mult by 0.65 ensures that lambda is more realisitc i.e closer to colebrook-white and other good approx*/

   z_m = MeanCompressibilityFactor(network, meanpressure, papay);

   /* To ensure that units are right in the weymouth equality we have to divide with 10^10 such that the
    * unit of LAMBDA is bar^2 * (kg/s)^-2. */

   /* missing values are added in the constraints  */
   LAMBDA = (16.0 * lambda * length * z_m * network->gasTemperature * universalGasConstant) / (pow( pi, 2.0 ) * pow( diameter, 5.0 ) * 1e10);

   return LAMBDA;
}


/** compute constant for nonlinear resistor */
SCIP_Real ConstantNonlinearResistor(
   GAS_Network*          network,            /**< network data structure */
   SCIP_Real             meanpressure,       /**< mean pressure value */
   SCIP_Real             dragfactor,         /**< dragfactor of pipe */
   SCIP_Real             diameter,           /**< diameter of pipe */
   SCIP_Bool             papay               /**< if formula of papay is used */
   )
{
   SCIP_Real z_m;
   SCIP_Real c_a;

   z_m = MeanCompressibilityFactor(network, meanpressure, papay);

   c_a = (8 * dragfactor * network->gasTemperature * z_m * universalGasConstant) / (pow( pi, 2.0 ) * pow( diameter, 4.0 ) * network->molarMass1 * 1e10);

   return c_a;
}

/** calculate slope
 *
 * Use length and heightdiff to calculate the adjacent side. Then divide heightdiff by adjacent to get the slope.
 */
SCIP_Real calculateSlope(
   SCIP_Real             heightdiff,         /**< height difference of innode to outnode */
   SCIP_Real             length              /**< pipe length */
   )
{
   SCIP_Real slope;
   SCIP_Real adjacent;

   adjacent = sqrt(pow(length, 2.0) - pow(heightdiff, 2.0));
   slope = heightdiff/adjacent;

   return slope;
}
