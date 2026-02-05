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

/**@file   unitConversion_gas.c
 * @brief  conversion of units used in stationary gas transport problems
 * @author Marc Pfetsch
 * @author Alexandra Stern
 * @author Imke Joormann
 */

#include "unitConversion_gas.h"
#include <string.h>
#include <scip/misc.h>

SCIP_RETCODE change_unit_meter(
   const char*           unit,
   SCIP_Real*            value
   )
{
   if ( strcmp(unit, "mm") == 0 )
   {
      *value= (*value) * 0.001;
   }
   else if ( strcmp(unit, "cm") == 0 )
   {
      *value= (*value) * 0.01;
#ifdef DEBUG
      SCIPdebugMessage( "Unit change from centimeter to meter.\n");
#endif
   }
   else if ( strcmp(unit, "km") == 0 )
   {
      *value= (*value) * 1000;
#ifdef DEBUG
      SCIPdebugMessage( "Unit change from kilometer to meter.\n");
#endif
   }
   else if ( strcmp(unit, "m") == 0 )
   {
#ifdef DEBUG
      SCIPdebugMessage( "The unit is meter.\n");
#endif
   }
   else
   {
      SCIPerrorMessage("Conversion of unit %s to m not implemented.\n", unit);
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

SCIP_RETCODE change_unit_bar(
   const char*           unit,
   double*               value
   )
{
   if ( strcmp(unit, "barg") == 0 )
   {
      *value = (*value) + 1.01325;
   }
   else if ( strcmp(unit, "bar") == 0 )
   {
      /* *value= (*value); */
   }
   else if ( strcmp(unit, "1") == 0 )
   {
      /* do nothing */
   }
   else
   {
      SCIPerrorMessage("Conversion of unit %s to bar not implemented.\n", unit);
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}

SCIP_RETCODE change_flow_unit(
   const char*           unit,
   SCIP_Real*            value,
   SCIP_Real             normDensity
   )
{
   if ( strcmp(unit, "1000m_cube_per_hour") == 0 )
   {
      *value= ((*value) / 3.6) * normDensity;
#ifdef DEBUG
      SCIPdebugMessage("Unit change from volumetric norm flow in 1000m_cube_per_hour to mass flow in kg_per_s.\n");
#endif
   }
   else if ( strcmp(unit, "1000000m_cube_per_day") == 0 )
   {
      *value= ((*value) / 0.864) * normDensity;
   }
   else if ( strcmp(unit, "kg_per_s") == 0  )
   {
      /* do nothing */
#ifdef DEBUG
      SCIPdebugMessage("The flow unit is kg_per_sec.\n");
#endif
   }
   else
   {
      SCIPerrorMessage("Conversion of unit %s to kg per second not implemented.\n", unit);
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}


SCIP_RETCODE change_density_unit(
   const char*           unit,
   SCIP_Real*            value
   )
{  /*lint --e{715} */
   assert( value != NULL );

   if ( strcmp(unit, "kg_per_m_cube") == 0 )
   {
      /* do nothing */
#ifdef DEBUG
      SCIPdebugMessage("The density unit is kg_per_m_cube.\n");
#endif
   }
   else
   {
      SCIPerrorMessage("Conversion of unit %s to kg per m cube not implemented.\n", unit);
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}

SCIP_RETCODE change_temperature_unit(
   const char*           unit,
   SCIP_Real*            value
   )
{
   if ( strcmp(unit, "K") == 0 )
   {
      /* do nothing */
#ifdef DEBUG
      SCIPdebugMessage("The temperature unit is Kelvin.\n");
#endif
   }
   else if ( strcmp(unit, "Celsius") == 0  )
   {
      *value = (*value) + 273.15;
#ifdef DEBUG
      SCIPdebugMessage("Unit change from Celsius to Kelvin.\n");
#endif
   }
   else
   {
      SCIPerrorMessage("Conversion of unit %s to Kelvin is not implemented.\n", unit);
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}

SCIP_RETCODE change_molarmass_unit(
   const char*           unit,
   SCIP_Real*               value
   )
{
   if ( strcmp(unit, "kg_per_mol") == 0 )
   {
      /* do nothing */
#ifdef DEBUG
      SCIPdebugMessage("The molarmass unit is kg_per_mol.\n");
#endif
   }
   else if ( strcmp(unit, "kg_per_kmol") == 0  )
   {
      *value = (*value) / 1000;
#ifdef DEBUG
      SCIPdebugMessage("Unit change from kg_per_kmol to kg_per_mol.\n");
#endif
   }
   else
   {
      SCIPerrorMessage("Conversion of unit %s to  kg_per_mol is not implemented.\n", unit);
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}

SCIP_RETCODE change_enthalpy_unit(
   const char*           unit,
   SCIP_Real*            value
   )
{
   if ( strcmp(unit, "J_per_kg") == 0 )
   {
      /* do nothing */
#ifdef DEBUG
      SCIPdebugMessage("The enthalpy unit is J_per_kg.\n");
#endif
   }
   else if ( strcmp(unit, "kJ_per_kg") == 0  )
   {
      *value= (*value) / 1000;
#ifdef DEBUG
      SCIPdebugMessage("Unit change from kJ_per_kg to J_per_kg.\n");
#endif
   }
   else
   {
      SCIPerrorMessage("Conversion of unit %s to  J_per_kg is not implemented.\n", unit);
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}

SCIP_RETCODE check_constant_unit(
   const char*           unit,
   SCIP_Real*            value
   )
{  /*lint --e{715} */
   assert( value != NULL );

   if ( strcmp(unit, "1") == 0 )
   {
      /* do nothing */
#ifdef DEBUG
      SCIPdebugMessage("The isentropic exponent unit is 1.\n");
#endif
   }
   else
   {
      SCIPerrorMessage("Conversion of unit %s to  1 is not implemented.\n", unit);
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}

SCIP_RETCODE change_gasconstant_unit(
   const char*           unit,
   SCIP_Real*            value
   )
{
   if ( strcmp(unit, "J_per_kg_per_K") == 0 )
   {
      /* do nothing */
#ifdef DEBUG
      SCIPdebugMessage("The enthalpy unit is J_per_kg.\n");
#endif
   }
   else if ( strcmp(unit, "kJ_per_kg_per_K") == 0  )
   {
      *value= (*value) / 1000;
#ifdef DEBUG
      SCIPdebugMessage("Unit change from kJ_per_kg to J_per_kg.\n");
#endif
   }
   else
   {
      SCIPerrorMessage("Conversion of unit %s to  J_per_kg_per_K is not implemented.\n", unit);
      return SCIP_INVALIDDATA;
   }
   return SCIP_OKAY;
}
