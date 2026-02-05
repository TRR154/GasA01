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

/**@file   LSFiles_gas.c
 * @brief  Creates and reads ls-file containing information about the network and the solution values of the variables
 * @author Alexandra Stern
 * @author Oliver Habeck
 * @author Dennis Gabriel
 */

#include <time.h>
#include <scip/misc.h>
#include <scip/scip.h>
#include <LSFiles_gas.h>
#include <elements_gas.h>
#include "unitConversion_gas.h"
#include "myxmldefs.h"
#include <string.h>
#include <scip/pub_message.h>
#include <xml/xml.h>

/** function to print information about the flow on an arc */
static
void printFlowOnArc(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_SOL*             sol,                /**< solution to be written */
   FILE*                 solfile,            /**< solution file */
   SCIP_VAR*             flowvar             /**< flow variable that should be printed */
   )
{
   assert( flowvar != NULL );
   fprintf(solfile, "        <flow unit=\"kg_per_s\" value=\"%.15f\"></flow>\n", SCIPgetSolVal(scip, sol, flowvar));
}

/** main function printing information */
SCIP_RETCODE writeLSF(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_SOL*             sol                 /**< solution to be written */
   )
{  /*lint --e{668,747}*/
   SCIP_PROBDATA* probdata;
   GAS_Network* network;
   FILE *solfile;
   GAS_Node* node;
   GAS_Arc* arc;
   GAS_CS* cs;
   GAS_Valve* valve;
   GAS_Controlvalve* cv;
   SCIP_Real pressure;
   SCIP_Real pdiff;
   SCIP_Real muNode;
   SCIP_Real muArc;
   SCIP_Real muN1;
   SCIP_Real muN2;
   SCIP_Real Z2;
   SCIP_Real Z1;
   SCIP_Real muMassArc;
   int binvar;
   char e1[3] = "  ";
   char e2[5] = "    ";
   char e3[7] = "      ";
   char e4[9] = "        ";
   char name[SCIP_MAXSTRLEN] = "";
   time_t t = time(NULL);
   struct tm tm = *localtime(&t);
   int i;
   int j;

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   network = probdata->network;
   assert( network != NULL);

   if ( sol == NULL )
   {
      SCIPdebugMessage( "No solution file created.\n");
      return SCIP_OKAY;
   }

   /* set the name of the solution file */
   strcat(name, probdata->netname);
   strcat(name, "-");
   strcat(name, probdata->scnname);
   strcat(name, ".lsf");

   solfile = fopen(name, "w");
   if ( solfile == NULL )
   {
      SCIPerrorMessage("Could not open file <%s>.\n", name);
      return SCIP_OKAY;
   }

   fprintf(solfile, "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\" ?>\n");

   /* print information about solution, model, objective, solver... */
   fprintf(solfile, "<solution");
   fprintf(solfile," network=\"%s\" nomination=\"%s\"", probdata->netname, probdata->scnname);
   if ( probdata->noObjective )
   {
      fprintf(solfile, " objective=\"no objective\"");
   }
   else
   {
      if ( probdata->minCompSum )
         fprintf(solfile, " objective=\"minimize number of running compressors\" value=\"%g\"", SCIPgetPrimalbound(scip));
      else if ( probdata->minPressSum )
         fprintf(solfile, " objective=\"minimize sum of pressures\" value=\"%f\"", SCIPgetPrimalbound(scip));
      else if ( probdata->powerLoss )
         fprintf(solfile, " objective=\"minimize the power loss\" value=\"%f\"", SCIPgetPrimalbound(scip));
      else if ( probdata->minSlack )
         fprintf(solfile, " objective=\"minimize maximal violation of pressure bounds\" value=\"%f\"", SCIPgetPrimalbound(scip));
      else if ( probdata->minSlackPerBound )
         fprintf(solfile, " objective=\"minimize sum of pressure bound violation\" value=\"%f\"", SCIPgetPrimalbound(scip));
      else if (probdata->maxFlow )
         fprintf(solfile, " objective=\"maximize the flow\" value=\"%f\" ", SCIPgetPrimalbound(scip));
      else
         fprintf(solfile, " objective=\"maximize sum of pressures\" value=\"%f\"", - SCIPgetPrimalbound(scip));

      if ( SCIPgetStatus(scip) == SCIP_STATUS_OPTIMAL )
      {
         fprintf(solfile, " status=\"optimal\"");
      }
      else
      {
         fprintf(solfile, " status=\"feasible\" gap=\"%f\"", SCIPgetGap(scip));
      }
   }

   /* print model options */
   if ( probdata->algebraic )
      fprintf(solfile, " pressureLossModel=\"algebraic\"");
   else
      fprintf(solfile, " pressureLossModel=\"ODE\"");

   if ( probdata->relaxLowerBounds )
      fprintf(solfile, " pressureBounds=\"relax lower pressure bounds by %3.1f bar\"", probdata->relaxLBvalue);
   if ( probdata->relaxUpperBounds )
      fprintf(solfile, " pressureBounds=\"relax upper pressure bounds by %3.1f bar\"", probdata->relaxUBvalue);
   if ( probdata->relaxCVbounds )
      fprintf(solfile, " pressureBounds=\"relax pressure bounds of controlvalves by %3.1f bar\"", probdata->relaxCVvalue);

   if ( probdata->network->numcompressor > 0 )
   {
      if ( probdata->boxConstraintModel )
         fprintf(solfile, " cs-model=\"box constraint\"");
      else
         fprintf(solfile, " cs-model=\"idealized\"");
   }

   if ( probdata->network->firstCycle != NULL )
   {
      if ( probdata->allCycles )
         fprintf(solfile, " noCircularFlowConstraints=\"use all cycles\"");
      else
         fprintf(solfile, " noCircularFlowConstraints=\"use cycle basis\"");
   }

   /* scip and solver information */
   fprintf(solfile, " SCIP-version=\"%d.%d.%d\"", SCIPmajorVersion(), SCIPminorVersion(), SCIPtechVersion());
#ifndef NDEBUG
   fprintf(solfile, " mode=\"debug\"");
#else
   fprintf(solfile, " mode=\"optimized\"");
#endif
   fprintf(solfile, " LP-solver=\"%s\"", SCIPlpiGetSolverName());
   /*  tm.tm_year starts with year 1900 as year 0, for some reason january is month 0 */
   fprintf(solfile, " date=\"%d.%d.%d\">\n", tm.tm_mday, tm.tm_mon +1 , tm.tm_year + 1900);

   /* print nodes (sources, sinks, innodes) */
   fprintf(solfile, "%s<nodes>\n", e1);
   fprintf(solfile, "%s<sources>\n", e2);

   /* entries/sources */
   node = &(network->nodes_ptr[0]);
   for (i = 0; i < network->numnodes; i++)
   {
      if ( node->type == ENTRY )
      {
         fprintf(solfile, "%s<source id=\"%s\">\n", e3, node->id);

         /* pressure is displayed in bar */
         pressure = SCIPgetSolVal(scip, sol, node->pressurevar);
         fprintf(solfile, "%s<pressure unit=\"bar\" value=\"%.15f\"></pressure>\n",e4,pressure);

         /* print pressure square variable */
         if ( node->needPIvar )
         {
            pressure = SCIPgetSolVal(scip, sol, node->PIvar);
            fprintf(solfile, "%s<pressureSquare unit=\"bar_square\" value=\"%.15f\"></pressureSquare>\n", e4, pressure);
         }

         /* Flow is displayed in kg_per_s */
         fprintf( solfile, "%s<flow unit=\"kg_per_s\" value=\"%.15f\"></flow>\n", e4, node->flow );

         if ( probdata->mixing || probdata->linearMixing )
         {
            /*Node mixing ratio in mass percent*/
            muNode = SCIPgetSolVal(scip, sol, node->nodeMixingRatio);
            fprintf(solfile, "%s<mixingratio unit=\"mass-percent\" value=\"%.15f\"></mixingratio>\n", e4, muNode);
         }

         if ( probdata->nodeMu && ! probdata->linearMixing )
         {
            /*Node mixing ratio in mol percent*/
            muNode = SCIPgetSolVal( scip, sol, node->nodeMixingRatioMol );
            fprintf( solfile, "%s<mixingratio unit=\"mol-percent\" value=\"%.15f\"></mixingratio>\n", e4, muNode );
         }
         fprintf(solfile, "%s</source>\n", e3);
      }
      ++node;
   }
   fprintf(solfile, "%s</sources>\n", e2);

   /* sinks */
   fprintf(solfile, "%s<sinks>\n", e2);
   node = &(network->nodes_ptr[0]);
   for (i=0; i<network->numnodes; i++)
   {
      if ( node->type == EXIT )
      {
         fprintf(solfile, "%s<sink id=\"%s\">\n", e3, node->id);

         /* pressure is displayed in bar */
         pressure = SCIPgetSolVal(scip, sol, node->pressurevar);
         fprintf(solfile,  "%s<pressure unit=\"bar\" value=\"%.15f\"></pressure>\n", e4, pressure);

         /* print pressure square variable */
         if ( node->needPIvar )
         {
            pressure = SCIPgetSolVal(scip, sol, node->PIvar);
            fprintf(solfile, "%s<pressureSquare unit=\"bar_square\" value=\"%.15f\"></pressureSquare>\n", e4, pressure);
         }

         /* Flow is displayed in kg_per_s */
         fprintf( solfile, "%s<flow unit=\"kg_per_s\" value=\"%.15f\"></flow>\n", e4, ABS( node->flow ) );

         if ( probdata->mixing || probdata->linearMixing )
         {
            /*Node mixing ratio in mass percent*/
            muNode = SCIPgetSolVal(scip, sol, node->nodeMixingRatio);
            fprintf(solfile, "%s<mixingratio unit=\"mass-percent\" value=\"%.15f\"></mixingratio>\n", e4, muNode);
         }

         if ( probdata->nodeMu && ! probdata->linearMixing )
         {
            /* Node mixing ratio in mol percent */
            muNode = SCIPgetSolVal( scip, sol, node->nodeMixingRatioMol );
            fprintf( solfile, "%s<mixingratio unit=\"mol-percent\" value=\"%.15f\"></mixingratio>\n", e4, muNode );
         }

         fprintf(solfile, "%s</sink>\n", e3);
      }
      ++node;
   }
   fprintf(solfile, "%s</sinks>\n", e2);

   /* innodes */
   fprintf(solfile, "%s<innodes>\n", e2);
   node = &(network->nodes_ptr[0]);
   for (i=0; i<network->numnodes; i++)
   {
      if ( node->type == INNODE )
      {
         fprintf(solfile, "%s<innode id=\"%s\">\n", e3, node->id);
         /* pressure is displayed in bar*/
         pressure = SCIPgetSolVal(scip, sol, node->pressurevar);
         fprintf(solfile, "%s<pressure unit=\"bar\" value=\"%.15f\"></pressure>\n", e4, pressure);

         /* print pressure square variable */
         if ( node->needPIvar )
         {
            pressure = SCIPgetSolVal(scip, sol, node->PIvar);
            fprintf(solfile, "%s<pressureSquare unit=\"bar_square\" value=\"%.15f\"></pressureSquare>\n", e4, pressure);
         }

         if ( probdata->mixing || probdata->linearMixing )
         {
            /*Node mixing ratio in mass percent*/
            muNode = SCIPgetSolVal(scip, sol, node->nodeMixingRatio);
            fprintf(solfile, "%s<mixingratio unit=\"mass-percent\" value=\"%.15f\"></mixingratio>\n", e4, muNode);
         }

         if (probdata->nodeMu && ! probdata->linearMixing)
         {
            /* Node mixing ratio in mol percent */
            muNode = SCIPgetSolVal( scip, sol, node->nodeMixingRatioMol );
            fprintf( solfile, "%s<mixingratio unit=\"mol-percent\" value=\"%.15f\"></mixingratio>\n", e4, muNode );
         }
         fprintf(solfile, "%s</innode>\n", e3);
      }
      ++node;
   }
   fprintf(solfile, "%s</innodes>\n", e2);
   fprintf(solfile, "%s</nodes>\n", e1);

   /* Print connections */
   fprintf(solfile, "%s<connections>\n", e1);

   /* print pipes with sourcenode, targetnode and flow */
   fprintf(solfile, "%s<pipes>\n", e2);
   arc = &(network->arcs_ptr[0]);
   for (i = 0; i < network->numarcs; i++)
   {
      if ( arc->type == PIPE )
      {
         fprintf(solfile, "%s<pipe id=\"%s\" from=\"%s\" to=\"%s\">\n", e3, arc->id, arc->sourcenode->id, arc->targetnode->id);
         printFlowOnArc(scip, sol, solfile, arc->flowvar);

         pdiff = SCIPgetSolVal(scip, sol, arc->sourcenode->pressurevar);
         pdiff -= SCIPgetSolVal(scip, sol, arc->targetnode->pressurevar);
         fprintf(solfile, "%s<pressureDifference unit=\"bar\" value=\"%.15f\"></pressureDifference>\n", e4, pdiff);

         if ( probdata->noFlowBinvars == FALSE )
         {
            binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, arc->negativeFlowBinvar) );
            fprintf(solfile, "%s<negFlowBinvar value=\"%d\"></negFlowBinvar>\n", e4, binvar);
            binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, arc->positiveFlowBinvar) );
            fprintf(solfile, "%s<posFlowBinvar value=\"%d\"></posFlowBinvar>\n", e4, binvar);
         }

         if ( probdata->nodeMu )
         {
            /* Arc mixing ratio in mole and mass percent */
            if ( ! probdata->linearMixing )
            {
               muN1 = SCIPgetSolVal(scip, sol, arc->sourcenode->nodeMixingRatioMol);
               muN2 = SCIPgetSolVal(scip, sol, arc->targetnode->nodeMixingRatioMol);
               Z1 = SCIPgetSolVal(scip, sol, arc->positiveFlowBinvar);
               Z2 = SCIPgetSolVal(scip, sol, arc->negativeFlowBinvar);
               muArc = Z1 * muN1 + Z2 * muN2;
               fprintf(solfile, "%s<mixingratioMol unit=\"mol_percent\" value=\"%.15f\"></mixingratioMol>\n", e4, muArc);
            }

            muN1 = SCIPgetSolVal( scip, sol, arc->sourcenode->nodeMixingRatio );
            muN2 = SCIPgetSolVal( scip, sol, arc->targetnode->nodeMixingRatio );
            Z1 = SCIPgetSolVal( scip, sol, arc->positiveFlowBinvar );
            Z2 = SCIPgetSolVal( scip, sol, arc->negativeFlowBinvar );
            muMassArc = Z1 * muN1 + Z2 * muN2;
            fprintf( solfile, "%s<mixingratioMass unit=\"mass_percent\" value=\"%.15f\"></mixingratioMass>\n", e4, muMassArc );

            /* LAMBDA */
            muArc = SCIPgetSolVal( scip, sol, arc->Lambda );
            fprintf( solfile, "%s<LAMBDA unit=\"resistance\" value=\"%.15f\"></LAMBDA>\n", e4, muArc );
         }
         else if ( probdata->mixing || probdata->linearMixing )
         {
            /* Arc mixing ratio in mole and mass percent*/
            if ( ! probdata->linearMixing )
            {
               muArc = SCIPgetSolVal(scip, sol, arc->mixingRatio);
               fprintf(solfile, "%s<mixingratioMol unit=\"mol_percent\" value=\"%.15f\"></mixingratioMol>\n", e4, muArc);
            }
            muMassArc = SCIPgetSolVal( scip, sol, arc->mixingRatioMass );
            fprintf( solfile, "%s<mixingratioMass unit=\"mass_percent\" value=\"%.15f\"></mixingratioMass>\n", e4, muMassArc );

            /*LAMBDA*/
            muArc = SCIPgetSolVal( scip, sol, arc->Lambda );
            fprintf( solfile, "%s<LAMBDA unit=\"resistance\" value=\"%.15f\"></LAMBDA>\n", e4, muArc );
         }

         fprintf(solfile, "%s</pipe>\n", e3);
      }
      ++arc;
   }
   fprintf(solfile, "%s</pipes>\n", e2);

    /* print shortpipes with sourcenode, targetnode and flow */
   if ( network->numshortpipes!= 0 )
   {
      fprintf(solfile, "%s<shortPipes>\n", e2);
      for (i=0; i<network->numarcs; i++)
      {
         arc = &(network->arcs_ptr[i]);
         if ( arc->type == SHORTPIPE )
         {
            fprintf(solfile, "%s<shortPipe id=\"%s\" from=\"%s\" to=\"%s\">\n", e3, arc->id, arc->sourcenode->id, arc->targetnode->id);
            printFlowOnArc(scip, sol, solfile, arc->flowvar);

            if ( probdata->noFlowBinvars == FALSE )
            {
               binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, arc->negativeFlowBinvar) );
               fprintf(solfile, "%s<negFlowBinvar value=\"%d\"></negFlowBinvar>\n", e4, binvar);
               binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, arc->positiveFlowBinvar) );
               fprintf(solfile, "%s<posFlowBinvar value=\"%d\"></posFlowBinvar>\n", e4, binvar);
            }
            fprintf(solfile, "%s</shortPipe>\n",e3);
         }
      }
      fprintf(solfile, "%s</shortPipes>\n", e2);
   }

   /* print valves with sourcenode, targetnode and flow */
   if ( network->numvalves != 0 )
   {
      SCIP_Bool printvalves = FALSE;

      for (i = 0; i < network->numarcs; i++)
      {
         arc = &(network->arcs_ptr[i]);

         if ( arc->type == VALVE )
         {
            valve = arc->detailed_info;
            /* skip bypasses */
            if (valve->bypassed_cs != NULL || valve->bypassed_cv != NULL )
               continue;

            /* only print valve section if there is one */
            if ( ! printvalves )
            {
               printvalves = TRUE;
               fprintf(solfile, "%s<valves>\n", e2);
            }

            fprintf(solfile, "%s<valve id=\"%s\" from=\"%s\" to=\"%s\">\n", e3, arc->id, arc->sourcenode->id, arc->targetnode->id);
            printFlowOnArc(scip, sol, solfile, arc->flowvar);
            binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, valve->valve_binvar) );
            fprintf(solfile, "%s<open value=\"%d\"></open>\n", e4, binvar);
            pdiff = SCIPgetSolVal(scip, sol, arc->sourcenode->pressurevar);
            pdiff -= SCIPgetSolVal(scip, sol, arc->targetnode->pressurevar);
            fprintf(solfile, "%s<pressureDifferential unit=\"bar\" value=\"%.15f\"></pressureDifferential>\n", e4, pdiff);
            if ( probdata->noFlowBinvars == FALSE )
            {
               binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, arc->negativeFlowBinvar) );
               fprintf(solfile, "%s<negFlowBinvar value=\"%d\"></negFlowBinvar>\n", e4, binvar);
               binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, arc->positiveFlowBinvar) );
               fprintf(solfile, "%s<posFlowBinvar value=\"%d\"></posFlowBinvar>\n", e4, binvar);
            }
            fprintf(solfile, "%s</valve>\n", e3);
         }
      }

      if ( printvalves )
         fprintf(solfile, "%s</valves>\n", e2);
   }

   /* print controlvalves with sourcenode, targetnode and flow */
   if ( network->numcontrolvalves != 0 )
   {
      fprintf(solfile, "%s<controlValves>\n", e2);

      for (i = 0; i < network->numarcs; i++)
      {
         arc = &(network->arcs_ptr[i]);
         if ( arc->type == CONTROLVALVE )
         {
            cv = arc->detailed_info;
            fprintf(solfile, "%s<controlValve id=\"%s\" from=\"%s\" to=\"%s\">\n", e3, arc->id, arc->sourcenode->id, arc->targetnode->id);
            binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, cv->binvar));

            if ( binvar == 1 )
            {
               printFlowOnArc(scip, sol, solfile, arc->flowvar);
               pdiff = SCIPgetSolVal(scip, sol, arc->sourcenode->pressurevar);
               pdiff -= SCIPgetSolVal(scip, sol, arc->targetnode->pressurevar);
               fprintf(solfile, "%s<active value=\"1\"></active>\n", e4);
               fprintf(solfile, "%s<open value=\"1\"></open>\n", e4);
               fprintf(solfile, "%s<pressureDifferential unit=\"bar\" value=\"%.15f\"></pressureDifferential>\n", e4, pdiff);
            }
            else if ( cv->internalBypassRequired == 1 )
            {
               printFlowOnArc(scip, sol, solfile, cv->bypass->flowvar);
               valve = cv->bypass->detailed_info;
               binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, valve->valve_binvar) );
               fprintf(solfile, "%s<active value=\"0\"></active>\n", e4);
               fprintf(solfile, "%s<open value=\"%d\"></open>\n", e4, binvar);
            }
            else
            {
               printFlowOnArc(scip, sol, solfile, arc->flowvar);
               fprintf(solfile, "%s<active value=\"0\"></active>\n", e4);
               fprintf(solfile, "%s<open value=\"0\"></open>\n", e4);
            }
            fprintf(solfile, "%s</controlValve>\n", e3);
         }
      }
      fprintf(solfile, "%s</controlValves>\n", e2);
   }

   /* print compressors */
   if ( network->numcompressor != 0 )
   {
      fprintf(solfile, "%s<compressorStations>\n", e2);

      for (i = 0; i < network->numarcs; i++)
      {
         arc = &(network->arcs_ptr[i]);
         if ( arc->type == CS )
         {
            cs = arc->detailed_info;
            fprintf(solfile,"%s<compressorStation id=\"%s\" from=\"%s\" to=\"%s\">\n", e3, arc->id, arc->sourcenode->id, arc->targetnode->id);

            binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, cs->compressor_binvar) );

            if ( binvar == 1 )
            {
               printFlowOnArc(scip, sol, solfile, arc->flowvar);
               if ( probdata->boxConstraintModel )
               {
                  for (j = 0; j < cs->numconfigurations; ++j)
                  {
                     binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, cs->boxcons[j].config_binvar) );
                     if ( binvar == 1 )
                        break;
                  }
                  fprintf(solfile, "%s<configuration id=\"%s\"></configuration>\n", e4, cs->boxcons[j].id);
               }
               else
                  fprintf(solfile, "%s<configuration id=\"active\"></configuration>\n", e4);

               pdiff = SCIPgetSolVal(scip, sol, arc->targetnode->pressurevar);
               pdiff -= SCIPgetSolVal(scip, sol, arc->sourcenode->pressurevar);
               fprintf(solfile,"%s<pressureIncrease value =\"%.15f\"></pressureIncrease>\n", e4, pdiff);
            }
            else
            {
               if ( cs->internalBypassRequired == 1 )
               {
                  SCIP_VAR* valve_binvar = probdata->VALVE_binvars[cs->bypassPosition];
                  binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, valve_binvar) );
               }

               if ( binvar == 1 )
               {
                  printFlowOnArc(scip, sol, solfile, network->arcs_ptr[i+1].flowvar); /*lint !e679 */
                  fprintf(solfile, "%s<configuration id=\"bypass\"></configuration>\n", e4);
               }
               else
               {
                  printFlowOnArc(scip, sol, solfile, arc->flowvar);
                  fprintf(solfile, "%s<configuration id=\"closed\"></configuration>\n", e4);
               }
            }

            /* print resistor data if existing */
            if ( cs->NLRin_pressure != NULL )
            {
               fprintf(solfile, "%s<resistorIn>\n", e4);

               pressure = SCIPgetSolVal(scip, sol, cs->NLRin_pressure);
               fprintf(solfile, "%s  <pressureOut unit=\"bar\" value=\"%.15f\"></pressureOut>\n", e4, pressure);
               pdiff = SCIPgetSolVal(scip, sol, cs->NLRin_pressureDiff);
               fprintf(solfile, "%s  <pressureDifferential unit=\"bar\" value=\"%.15f\"></pressureDifferential>\n", e4, pdiff);

               fprintf(solfile, "%s</resistorIn>\n", e4);
            }

            if ( cs->NLRout_pressure != NULL )
            {
               fprintf(solfile, "%s<resistorOut>\n", e4);

               pressure = SCIPgetSolVal(scip, sol, cs->NLRout_pressure);
               fprintf(solfile, "%s  <pressureIn unit=\"bar\" value=\"%.15f\"></pressureIn>\n", e4, pressure);
               pdiff = SCIPgetSolVal(scip, sol, cs->NLRout_pressureDiff);
               fprintf(solfile, "%s  <pressureDifferential unit=\"bar\" value=\"%.15f\"></pressureDifferential>\n", e4, pdiff);

               fprintf(solfile, "%s</resistorOut>\n", e4);
            }

            fprintf(solfile, "%s</compressorStation>\n", e3);
         }
      }

      fprintf(solfile, "%s</compressorStations>\n", e2);
   }

   /* print resistor information*/
   if ( network->numlinresistor != 0 || network->numnonlinresistor != 0 )
   {
      fprintf(solfile, "%s<resistors>\n", e2);

      for (i = 0; i < network->numarcs; i++)
      {
         arc = &(network->arcs_ptr[i]);
         if ( arc->type == RESISTOR )
         {
            fprintf(solfile,"%s<resistor id=\"%s\" from=\"%s\" to=\"%s\">\n", e3, arc->id, arc->sourcenode->id, arc->targetnode->id);
            printFlowOnArc(scip, sol, solfile, arc->flowvar);
            pdiff = SCIPgetSolVal(scip, sol, arc->sourcenode->pressurevar) - SCIPgetSolVal(scip, sol, arc->targetnode->pressurevar);
            fprintf(solfile, "%s<pressureDifferential unit=\"bar\" value=\"%.15f\"></pressureDifferential>\n", e4, pdiff);
            if ( probdata->noFlowBinvars == FALSE )
            {
               binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, arc->negativeFlowBinvar) );
               fprintf(solfile, "%s<negFlowBinvar value=\"%d\"></negFlowBinvar>\n", e4, binvar);
               binvar = (int) SCIPround(scip, SCIPgetSolVal(scip, sol, arc->positiveFlowBinvar) );
               fprintf(solfile, "%s<posFlowBinvar value=\"%d\"></posFlowBinvar>\n", e4, binvar);
            }
            fprintf(solfile, "%s</resistor>\n", e3);
         }
      }

      fprintf(solfile, "%s</resistors>\n", e2);
   }
   fprintf(solfile, "%s</connections>\n", e1);
   fprintf(solfile, "</solution>\n");
   fclose(solfile);

   SCIPinfoMessage(scip, NULL, "\nWritten solution to LS file <%s>.\n", name);

   return SCIP_OKAY;
}


/** reads solution data for a pipe given by the xml_node */
static
SCIP_RETCODE readPipeSolution(
   SCIP*                 scip,               /**< SCIP instance */
   const XML_NODE*       xml_arc             /**< XML block for arc */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_HASHTABLE* ArcsTable;
   GAS_Arc* arc;
   const XML_NODE* xml_node;
   const char* id;
   const char* value;
   const char* unit;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool fixed = FALSE;
   SCIP_Real q;
   SCIP_Real qlb;
   SCIP_Real qub;
   int pbin = 0;
   int nbin = 0;
   int lb;
   int ub;

   assert( scip != NULL );
   assert( xml_arc != NULL );

   probdata = SCIPgetProbData(scip);
   ArcsTable = probdata->network->Arcshashtable;

   /* find variable id */
   id = SCIPxmlGetAttrval(xml_arc, "id");
   if ( id == NULL )
   {
      SCIPerrorMessage("Attribute \"id\" of pipe not found.\n");
      return SCIP_READERROR;
   }

   /* retrieve GAS_Node from hashtable */
   arc = (GAS_Arc*) SCIPhashtableRetrieve(ArcsTable, (void*) id);
   if ( arc == NULL )
   {
      SCIPerrorMessage("Could not find GAS_Arc <%s> given in the solution file. Maybe solution file does not match network?\n", id);
      return SCIP_READERROR;
   }

   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "flow", 0, 1);
   if ( xml_node == NULL )
   {
      SCIPerrorMessage("Could not find the flow on pipe <%s>.\n", id);
      return SCIP_READERROR;
   }
   else
   {
      unit = SCIPxmlGetAttrval(xml_node, "unit");
      if ( unit == NULL )
      {
         SCIPerrorMessage("Flow unit of pipe <%s> not found.\n", id);
         return SCIP_READERROR;
      }
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Flow value of pipe <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      q = atof(value);
      SCIP_CALL( change_flow_unit(unit, &q, probdata->network->normDensity) );

      qlb = SCIPcomputeVarLbGlobal(scip, arc->flowvar);
      qub = SCIPcomputeVarUbGlobal(scip, arc->flowvar);

      /* fix flow */
      if ( SCIPisFeasGT(scip, qlb, q) || SCIPisFeasLT(scip, qub, q) )
      {
         SCIPerrorMessage("Flow on pipe <%s> given in solution file does not satisfy the variable bounds!\n", id);
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPfixVar(scip, arc->flowvar, q, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
   }

   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "negFlowBinvar", 0, 1);
   if ( xml_node != NULL )
   {
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Negative flow binvar of pipe <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      nbin = atoi(value);

      lb = (int) SCIPcomputeVarLbGlobal(scip, arc->negativeFlowBinvar);
      ub = (int) SCIPcomputeVarUbGlobal(scip, arc->negativeFlowBinvar);

      if ( lb <= nbin && nbin <= ub )
      {
         SCIP_CALL( SCIPfixVar(scip, arc->negativeFlowBinvar, (SCIP_Real) nbin, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
      else
      {
         SCIPerrorMessage("Cannot fix flowBinvar <%s> to %d, since this contradicts the variable bounds.\n", SCIPvarGetName(arc->negativeFlowBinvar), nbin);
         return SCIP_INVALIDDATA;
      }
   }

   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "posFlowBinvar", 0, 1);
   if ( xml_node != NULL )
   {
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Positive flow binvar of pipe <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      pbin = atoi(value);

      lb = (int) SCIPcomputeVarLbGlobal(scip, arc->positiveFlowBinvar);
      ub = (int) SCIPcomputeVarUbGlobal(scip, arc->positiveFlowBinvar);

      if ( lb <= pbin && pbin <= ub )
      {
         SCIP_CALL( SCIPfixVar(scip, arc->positiveFlowBinvar, (SCIP_Real) pbin, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
      else
      {
         SCIPerrorMessage("Cannot fix flowBinvar <%s> to %d, since this contradicts the variable bounds.\n", SCIPvarGetName(arc->positiveFlowBinvar), pbin);
         return SCIP_INVALIDDATA;
      }
   }

   if ( nbin+pbin == 2 )
   {
      SCIPerrorMessage("Positive and Negative flow binvars of pipe <%s> cannot be both 1.\n", id);
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** reads solution data of a valve given by the xml_node */
static
SCIP_RETCODE readValveSolution(
   SCIP*                 scip,               /**< SCIP instance */
   const XML_NODE*       xml_arc             /**< XML block for arc */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_HASHTABLE* ArcsTable;
   GAS_Arc* arc;
   GAS_Valve* valve;
   const XML_NODE* xml_node;
   const char* id;
   const char* value;
   const char* unit;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool fixed = FALSE;
   SCIP_Bool fix_pbin = FALSE;
   SCIP_Bool fix_nbin = FALSE;
   SCIP_Real q;
   SCIP_Real pdiff;
   int valve_binvar;
   int pbin = 0;
   int nbin = 0;
   int lb;
   int ub;

   assert( scip != NULL );
   assert( xml_arc != NULL );

   probdata = SCIPgetProbData(scip);
   ArcsTable = probdata->network->Arcshashtable;

   /* find variable id */
   id = SCIPxmlGetAttrval(xml_arc, "id");
   if ( id == NULL )
   {
      SCIPerrorMessage("Attribute \"id\" of valve not found.\n");
      return SCIP_READERROR;
   }

   /* retrieve GAS_Node from hashtable */
   arc = (GAS_Arc*) SCIPhashtableRetrieve(ArcsTable, (void*) id);
   if ( arc == NULL )
   {
      SCIPerrorMessage("Could not find GAS_Arc <%s> given in the solution file. Maybe solution file does not match network?\n", id);
      return SCIP_READERROR;
   }
   valve = arc->detailed_info;

   /* find status of valve */
   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "open", 0, 1);
   if ( xml_node != NULL )
   {
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Status of valve <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      valve_binvar = atoi(value);

      lb = (int) SCIPcomputeVarLbGlobal(scip, valve->valve_binvar);
      ub = (int) SCIPcomputeVarUbGlobal(scip, valve->valve_binvar);

      if ( lb <= valve_binvar && valve_binvar <= ub )
      {
         SCIP_CALL( SCIPfixVar(scip, valve->valve_binvar, (SCIP_Real) valve_binvar, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
   }
   else
   {
      SCIPerrorMessage("Status of valve <%s> not found.\n", id);
      return SCIP_READERROR;
   }

   /* set/find values for flow and flow binvars */
   if ( valve_binvar == 0 )
   {
      /* there can be no flow and flow binvars have to be 0 */
      q = 0.0;
      pbin = 0;
      nbin = 0;
      fix_pbin = TRUE;
      fix_nbin = TRUE;
   }
   else
   {
      /* find flow value */
      xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "flow", 0, 1);
      if ( xml_node == NULL )
      {
         SCIPerrorMessage("Could not find the flow on valve <%s>.\n", id);
         return SCIP_READERROR;
      }
      else
      {
         unit = SCIPxmlGetAttrval(xml_node, "unit");
         if ( unit == NULL )
         {
            SCIPerrorMessage("Flow unit of valve <%s> not found.\n", id);
            return SCIP_READERROR;
         }
         value = SCIPxmlGetAttrval(xml_node, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("Flow value of valve <%s> not found.\n", id);
            return SCIP_READERROR;
         }
         q = atof(value);

         SCIP_CALL( change_flow_unit(unit, &q, probdata->network->normDensity) );
      }

      /* find negative flow binvar if lsf contains it */
      xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "negFlowBinvar", 0, 1);
      if ( xml_node != NULL )
      {
         value = SCIPxmlGetAttrval(xml_node, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("Negative flow binvar of valve <%s> not found.\n", id);
            return SCIP_READERROR;
         }

         fix_nbin = TRUE;
         nbin = atoi(value);
      }

      /* find positive flow binvar if lsf contains it */
      xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "posFlowBinvar", 0, 1);
      if ( xml_node != NULL )
      {
         value = SCIPxmlGetAttrval(xml_node, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("Positive flow binvar of valve <%s> not found.\n", id);
            return SCIP_READERROR;
         }

         fix_pbin = TRUE;
         pbin = atoi(value);
      }
   }

   /* find pressure differential if lsf contains it */
   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "pressureDifferential", 0, 1);
   if ( xml_node != NULL )
   {
      unit = SCIPxmlGetAttrval(xml_node, "unit");
      if ( unit == NULL )
      {
         SCIPerrorMessage("Unit of pressureDifferential at valve <%s> not found.\n", id);
         return SCIP_READERROR;
      }
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Value of pressureDifferential at valve <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      pdiff = atof(value);
      SCIP_CALL( change_unit_bar(unit, &pdiff) );
   }
   else
   {
      pdiff = SCIPcomputeVarUbGlobal(scip, arc->sourcenode->pressurevar);
      pdiff -= SCIPcomputeVarUbGlobal(scip, arc->targetnode->pressurevar);
   }

   /* pressure differential consistency check */
   if ( valve_binvar )
   {
      if ( ! SCIPisZero(scip, pdiff) )
      {
         SCIPerrorMessage("Pressure differential on valve <%s> not 0, although valve is open.\n", id);
         return SCIP_INVALIDDATA;
      }
   }
   else
   {
      if ( ! SCIPisFeasLE(scip, ABS(pdiff), valve->pressureDifferentialMax) )
      {
         SCIPerrorMessage("Maximal pressure differential of valve <%s> not satisfied.\n", id);
         return SCIP_INVALIDDATA;
      }
   }

   /* check if values are possible and fix variables */
   if ( SCIPisFeasGT(scip, arc->flowMin,q) || SCIPisFeasLT(scip, arc->flowMax, q) )
   {
      SCIPerrorMessage("Flow on valve <%s> given in solution file does not satisfy the variable bounds!\n", id);
      return SCIP_INVALIDDATA;
   }
   else
   {
      SCIP_CALL( SCIPfixVar(scip, arc->flowvar, q, &infeasible, &fixed) );
      assert( !infeasible );
      assert( fixed );
   }

   if ( nbin + pbin == 2 )
   {
      SCIPerrorMessage("Positive and Negative flow binvars of valve <%s> cannot be both 1.\n", id);
      return SCIP_INVALIDDATA;
   }

   if ( fix_nbin )
   {
      if ( nbin == 0 && SCIPisFeasNegative(scip, q) )
      {
         SCIPerrorMessage("Cannot fix flowBinvar <%s> to 0, since flow given in solution file is negative.\n", SCIPvarGetName(arc->negativeFlowBinvar));
         return SCIP_INVALIDDATA;
      }

      lb = (int) SCIPcomputeVarLbGlobal(scip, arc->negativeFlowBinvar);
      ub = (int) SCIPcomputeVarUbGlobal(scip, arc->negativeFlowBinvar);
      if ( lb <= nbin && nbin <= ub )
      {
         SCIP_CALL( SCIPfixVar(scip, arc->negativeFlowBinvar, (SCIP_Real) nbin, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
      else
      {
         SCIPerrorMessage("Cannot fix flowBinvar <%s> to %d, since this contradicts the variable bounds.\n", SCIPvarGetName(arc->negativeFlowBinvar), nbin);
         return SCIP_INVALIDDATA;
      }
   }

   if ( fix_pbin )
   {
      if ( pbin == 0 && SCIPisFeasPositive(scip, q) )
      {
         SCIPerrorMessage("Cannot fix flowBinvar <%s> to 0, since flow given in solution file is positive.\n", SCIPvarGetName(arc->positiveFlowBinvar));
         return SCIP_INVALIDDATA;
      }

      lb = (int) SCIPcomputeVarLbGlobal(scip, arc->positiveFlowBinvar);
      ub = (int) SCIPcomputeVarUbGlobal(scip, arc->positiveFlowBinvar);

      if ( lb <= pbin && pbin <= ub )
      {
         SCIP_CALL( SCIPfixVar(scip, arc->positiveFlowBinvar, (SCIP_Real) pbin, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
      else
      {
         SCIPerrorMessage("Cannot fix flowBinvar <%s> to %d, since this contradicts the variable bounds.\n", SCIPvarGetName(arc->positiveFlowBinvar), pbin);
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** reads solution data of a controlvalve given by the xml_node */
static
SCIP_RETCODE readCvSolution(
   SCIP*                 scip,               /**< SCIP instance */
   const XML_NODE*       xml_arc             /**< XML block for arc */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_HASHTABLE* ArcsTable;
   GAS_Arc* arc;
   GAS_Arc* bypass_arc;
   GAS_Controlvalve* cv;
   GAS_Valve* bypass_valve;
   const XML_NODE* xml_node;
   const char* id;
   const char* value;
   const char* unit;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool fixed = FALSE;
   SCIP_Real q;
   SCIP_Real qlb;
   SCIP_Real qub;
   SCIP_Real flow_cv;
   SCIP_Real flow_bypass;
   SCIP_Real pdiff;
   int cv_open, cv_active;
   int lb;
   int ub;

   assert( scip != NULL );
   assert( xml_arc != NULL );

   probdata = SCIPgetProbData(scip);
   ArcsTable = probdata->network->Arcshashtable;

   /* find variable id */
   id = SCIPxmlGetAttrval(xml_arc, "id");
   if ( id == NULL )
   {
      SCIPerrorMessage("Attribute \"id\" of valve not found.\n");
      return SCIP_READERROR;
   }

   /* retrieve GAS_Node from hashtable */
   arc = (GAS_Arc*) SCIPhashtableRetrieve(ArcsTable, (void*) id);
   if ( arc == NULL )
   {
      SCIPerrorMessage("Could not find GAS_Arc <%s> given in the solution file. Maybe solution file does not match network?\n", id);
      return SCIP_READERROR;
   }
   cv = arc->detailed_info;

   if ( cv->internalBypassRequired )
   {
      bypass_arc = cv->bypass;
      bypass_valve = bypass_arc->detailed_info;
   }
   else
   {
      bypass_arc = NULL;
      bypass_valve = NULL;
   }

   /* find flow value */
   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "flow", 0, 1);
   if ( xml_node == NULL )
   {
      SCIPerrorMessage("Could not find the flow on controlvalve <%s>.\n", id);
      return SCIP_READERROR;
   }
   else
   {
      unit = SCIPxmlGetAttrval(xml_node, "unit");
      if ( unit == NULL )
      {
         SCIPerrorMessage("Flow unit of controlvalve <%s> not found.\n", id);
         return SCIP_READERROR;
      }
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Flow value of controlvalve <%s> not found.\n", id);
         return SCIP_READERROR;
      }
      q = atof(value);

      SCIP_CALL( change_flow_unit(unit, &q, probdata->network->normDensity) );
   }

   /* find status active of cv */
   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "active", 0, 1);
   if ( xml_node != NULL )
   {
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Status \"active\" of controlvalve <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      cv_active = atoi(value);
   }
   else
   {
      SCIPerrorMessage("Status \"active\" of controlvalve <%s> not found.\n", id);
      return SCIP_READERROR;
   }

   /* find status open of cv */
   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "open", 0, 1);
   if ( xml_node != NULL )
   {
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Status \"open\" of controlvalve <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      cv_open = atoi(value);
   }
   else
   {
      SCIPerrorMessage("Status \"open\" of controlvalve <%s> not found.\n", id);
      return SCIP_READERROR;
   }

   /* plausibility check and set binary variables */
   if ( (cv_active == 1) && (cv_open == 0) )
   {
      SCIPerrorMessage("Controlvalve <%s> can not be active if it is not open!\n", id);
      return SCIP_INVALIDDATA;
   }
   else if ( (cv_active == 0) && (cv_open == 1) && (bypass_arc == NULL) )
   {
      SCIPerrorMessage("Controlvalve <%s> can not be in bypass if there is none!\n", id);
      return SCIP_INVALIDDATA;
   }

   lb = (int) SCIPcomputeVarLbGlobal(scip, cv->binvar);
   ub = (int) SCIPcomputeVarUbGlobal(scip, cv->binvar);

   if ( (lb > cv_active) || (cv_active > ub) )
   {
      SCIPerrorMessage("Can not set binary variable of controlvalve <%s> to %d.\n", id, cv_active);
      return SCIP_INVALIDDATA;
   }
   else
   {
      SCIP_CALL( SCIPfixVar(scip, cv->binvar, (SCIP_Real) cv_active, &infeasible, &fixed) );
      assert( !infeasible );
      assert( fixed );
   }

   if ( bypass_valve != NULL )
   {
      int tmp = MIN(1 - cv_active, cv_open);

      lb = (int) SCIPcomputeVarLbGlobal(scip, bypass_valve->valve_binvar);
      ub = (int) SCIPcomputeVarUbGlobal(scip, bypass_valve->valve_binvar);

      if ( (lb > tmp) || (tmp > ub) )
      {
         SCIPerrorMessage("Can not set binary variable of controlvalve <%s> bypass to %d.\n", id, tmp);
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPfixVar(scip, bypass_valve->valve_binvar, (SCIP_Real) tmp, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
   }

   /* check and set flow variables */
   if ( cv_active )
   {
      flow_cv = q;
      flow_bypass = 0.0;

      if ( SCIPisFeasNegative(scip, q) )
      {
         SCIPerrorMessage("Flow over controlvalve <%s> has to be nonnegative.\n", id);
         return SCIP_INVALIDDATA;
      }
   }
   else if ( cv_open )
   {
      flow_cv = 0.0;
      flow_bypass = q;
   }
   else
   {
      flow_cv = 0.0;
      flow_bypass = 0.0;

      if ( ! SCIPisFeasZero(scip, q) )
      {
         SCIPerrorMessage("Flow over controlvalve <%s> has to be zero if it is closed.\n", id);
         return SCIP_INVALIDDATA;
      }
   }

   qlb = SCIPcomputeVarLbGlobal(scip, arc->flowvar);
   qub = SCIPcomputeVarUbGlobal(scip, arc->flowvar);

   if ( SCIPisFeasGT(scip, qlb, flow_cv) || SCIPisFeasLT(scip, qub, flow_cv) )
   {
      SCIPerrorMessage("Flow over controlvalve <%s> does not satisfy the variable bounds.\n", id);
      return SCIP_INVALIDDATA;
   }
   else
   {
      SCIP_CALL( SCIPfixVar(scip, arc->flowvar, flow_cv, &infeasible, &fixed) );
      assert( !infeasible );
      assert( fixed );
   }

   if ( bypass_arc != NULL )
   {
      qlb = SCIPcomputeVarLbGlobal(scip, bypass_arc->flowvar);
      qub = SCIPcomputeVarUbGlobal(scip, bypass_arc->flowvar);

      if ( SCIPisFeasGT(scip, qlb, flow_bypass) || SCIPisFeasLT(scip, qub, flow_bypass) )
      {
         SCIPerrorMessage("Flow over controlvalve <%s> bypass does not satisfy the variable bounds.\n", id);
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPfixVar(scip, bypass_arc->flowvar, flow_bypass, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
   }

   /* find pressure differential if lsf contains it */
   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "pressureDifferential", 0, 1);
   if ( xml_node != NULL )
   {
      unit = SCIPxmlGetAttrval(xml_node, "unit");
      if ( unit == NULL )
      {
         SCIPerrorMessage("Unit of pressureDifferential of controlvalve <%s> not found.\n", id);
         return SCIP_READERROR;
      }
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Value of pressureDifferential of  controlvalve <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      pdiff = atof(value);
      SCIP_CALL( change_unit_bar(unit, &pdiff) );
   }
   else
   {
      pdiff = SCIPcomputeVarUbGlobal(scip, arc->sourcenode->pressurevar);
      pdiff -= SCIPcomputeVarUbGlobal(scip, arc->targetnode->pressurevar);
   }

   /* check if pressure diff satisfies the bounds */
   if ( cv_active )
   {
      /* @note: this is only correct if there are no nonlinear resistors in the cv station */
      if ( SCIPisFeasGT(scip, pdiff, cv->pressureDifferentialMax + cv->pressureLossIn + cv->pressureLossOut) )
      {
         SCIPerrorMessage("Pressure differential over controlvalve <%s> is too big.\n", id);
         return SCIP_INVALIDDATA;
      }

      if ( SCIPisFeasLT(scip, pdiff, cv->pressureDifferentialMin + cv->pressureLossIn + cv->pressureLossOut) )
      {
         SCIPerrorMessage("Pressure differential over controlvalve <%s> is too small.\n", id);
         return SCIP_INVALIDDATA;
      }
   }
   else if ( cv_open )
   {
      if ( ! SCIPisFeasZero(scip, pdiff) )
      {
         SCIPerrorMessage("Pressure differential over controlvalve <%s> has to be zero if it is in bypass.\n", id);
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}

/** reads solution data of a resistor given by the xml_node */
static
SCIP_RETCODE readResistorSolution(
   SCIP*                 scip,               /**< SCIP instance */
   const XML_NODE*       xml_arc             /**< XML block for arc */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_HASHTABLE* ArcsTable;
   GAS_Arc* arc;
   GAS_Resistor* res;
   const XML_NODE* xml_node;
   const char* id;
   const char* value;
   const char* unit;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool fixed = FALSE;
   SCIP_Real q;
   SCIP_Real qlb;
   SCIP_Real qub;
   SCIP_Real pdiff;
   int pbin = 0;
   int nbin = 0;
   int lb;
   int ub;

   assert( scip != NULL );
   assert( xml_arc != NULL );

   probdata = SCIPgetProbData(scip);
   ArcsTable = probdata->network->Arcshashtable;

   /* find variable id */
   id = SCIPxmlGetAttrval(xml_arc, "id");
   if ( id == NULL )
   {
      SCIPerrorMessage("Attribute \"id\" of resistor not found.\n");
      return SCIP_READERROR;
   }

   /* retrieve GAS_Node from hashtable */
   arc = (GAS_Arc*) SCIPhashtableRetrieve(ArcsTable, (void*) id);
   if ( arc == NULL )
   {
      SCIPerrorMessage("Could not find GAS_Arc <%s> given in the solution file. Maybe solution file does not match network?\n", id);
      return SCIP_READERROR;
   }
   res = arc->detailed_info;

   /* find flow */
   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "flow", 0, 1);
   if ( xml_node == NULL )
   {
      SCIPerrorMessage("Could not find the flow on resistor <%s>.\n", id);
      return SCIP_READERROR;
   }
   else
   {
      unit = SCIPxmlGetAttrval(xml_node, "unit");
      if ( unit == NULL )
      {
         SCIPerrorMessage("Flow unit of resistor <%s> not found.\n", id);
         return SCIP_READERROR;
      }
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Flow value of resistor <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      q = atof(value);
      SCIP_CALL( change_flow_unit(unit, &q, probdata->network->normDensity) );

      qlb = SCIPcomputeVarLbGlobal(scip, arc->flowvar);
      qub = SCIPcomputeVarUbGlobal(scip, arc->flowvar);

      /* fix flow */
      if ( SCIPisFeasGT(scip, qlb, q) || SCIPisFeasLT(scip, qub, q) )
      {
         SCIPerrorMessage("Flow on resistor <%s> given in solution file does not satisfy the variable bounds!\n", id);
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPfixVar(scip, arc->flowvar, q, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
   }

   /* find binary flow values if existing */
   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "negFlowBinvar", 0, 1);
   if ( xml_node != NULL )
   {
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Negative flow binvar of resistor <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      nbin = atoi(value);

      lb = (int) SCIPcomputeVarLbGlobal(scip, arc->negativeFlowBinvar);
      ub = (int) SCIPcomputeVarUbGlobal(scip, arc->negativeFlowBinvar);

      if ( (nbin == 1 && SCIPisFeasPositive(scip, q)) || (nbin == 0 && SCIPisFeasNegative(scip, q)) )
      {
         SCIPerrorMessage("Sign of flow on resistor <%s> does not coincide with negative flow binvar.\n", id);
         return SCIP_INVALIDDATA;
      }

      if ( lb <= nbin && nbin <= ub )
      {
         SCIP_CALL( SCIPfixVar(scip, arc->negativeFlowBinvar, (SCIP_Real) nbin, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
      else
      {
         SCIPerrorMessage("Cannot fix flowBinvar <%s> to %d, since this contradicts the variable bounds.\n", SCIPvarGetName(arc->negativeFlowBinvar), nbin);
         return SCIP_INVALIDDATA;
      }
   }

   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "posFlowBinvar", 0, 1);
   if ( xml_node != NULL )
   {
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Positive flow binvar of resistor <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      pbin = atoi(value);

      lb = (int) SCIPcomputeVarLbGlobal(scip, arc->positiveFlowBinvar);
      ub = (int) SCIPcomputeVarUbGlobal(scip, arc->positiveFlowBinvar);

      if ( (pbin == 0 && SCIPisFeasPositive(scip, q)) || (pbin == 1 && SCIPisFeasNegative(scip, q)) )
      {
         SCIPerrorMessage("Sign of flow on resistor <%s> does not coincide with positive flow binvar.\n", id);
         return SCIP_INVALIDDATA;
      }

      if ( lb <= pbin && pbin <= ub )
      {
         SCIP_CALL( SCIPfixVar(scip, arc->positiveFlowBinvar, (SCIP_Real) pbin, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
      else
      {
         SCIPerrorMessage("Cannot fix flowBinvar <%s> to %d, since this contradicts the variable bounds.\n", SCIPvarGetName(arc->positiveFlowBinvar), pbin);
         return SCIP_INVALIDDATA;
      }
   }

   if ( nbin + pbin == 2 )
   {
      SCIPerrorMessage("Positive and Negative flow binvars of resistor <%s> cannot be both 1.\n", id);
      return SCIP_INVALIDDATA;
   }

   /* if resistor is nonlinear set the additional variables */
   if ( res->type == NONLINEAR )
   {
      /* find pressure differential if lsf contains it */
      xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "pressureDifferential", 0, 1);
      if ( xml_node != NULL )
      {
         unit = SCIPxmlGetAttrval(xml_node, "unit");
         if ( unit == NULL )
         {
            SCIPerrorMessage("Unit of pressureDifferential of resistor <%s> not found.\n", id);
            return SCIP_READERROR;
         }
         value = SCIPxmlGetAttrval(xml_node, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("Value of pressureDifferential of resistor <%s> not found.\n", id);
            return SCIP_READERROR;
         }

         pdiff = atof(value);
         SCIP_CALL( change_unit_bar(unit, &pdiff) );
      }
      else
      {
         SCIPerrorMessage("Resistor solution data of resistor <%s> does not contain the pressure differential.\n", id);
         pdiff = SCIPcomputeVarUbGlobal(scip, arc->sourcenode->pressurevar);
         pdiff -= SCIPcomputeVarUbGlobal(scip, arc->targetnode->pressurevar);
      }

      SCIP_CALL( SCIPfixVar(scip, probdata->NLR_DeltaVars[res->NLR_position], pdiff, &infeasible, &fixed) );
      assert( !infeasible );
      assert( fixed );

      pdiff *= ABS(pdiff);
      SCIP_CALL( SCIPfixVar(scip, probdata->NLR_AbsDeltaVars[res->NLR_position], pdiff, &infeasible, &fixed) );
      assert( !infeasible );
      assert( fixed );

      q *= ABS(q);
      SCIP_CALL( SCIPfixVar(scip, probdata->NLR_AbsFlowVars[res->NLR_position], q, &infeasible, &fixed) );
      assert( !infeasible );
      assert( fixed );
   }

   return SCIP_OKAY;
}

/** reads solution data for a pipe given by the xml_node */
static
SCIP_RETCODE readCSsolution(
   SCIP*                 scip,               /**< SCIP instance */
   const XML_NODE*       xml_arc             /**< XML block for arc */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_HASHTABLE* ArcsTable;
   GAS_Arc* arc;
   GAS_CS* cs;
   const XML_NODE* xml_node;
   const XML_NODE* xml_child;
   const char* id;
   const char* conf_id = NULL;
   const char* value;
   const char* unit;
   SCIP_Bool infeasible = FALSE;
   SCIP_Bool fixed = FALSE;
   SCIP_Real q;
   SCIP_Real p;
   SCIP_Real pdiff;
   int bin;

   assert( scip != NULL );
   assert( xml_arc != NULL );

   probdata = SCIPgetProbData(scip);
   ArcsTable = probdata->network->Arcshashtable;

   /* find variable id */
   id = SCIPxmlGetAttrval(xml_arc, "id");
   if ( id == NULL )
   {
      SCIPerrorMessage("Attribute \"id\" of compressorStation not found.\n");
      return SCIP_READERROR;
   }

   /* retrieve GAS_Node from hashtable */
   arc = (GAS_Arc*) SCIPhashtableRetrieve(ArcsTable, (void*) id);
   if ( arc == NULL )
   {
      SCIPerrorMessage("Could not find GAS_Arc <%s> given in the solution file. Maybe solution file does not match network?\n", id);
      return SCIP_READERROR;
   }
   cs = arc->detailed_info;

   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "flow", 0, 1);
   if ( xml_node == NULL )
   {
      SCIPerrorMessage("Could not find the flow on compressorStation <%s>.\n", id);
      return SCIP_READERROR;
   }
   else
   {
      unit = SCIPxmlGetAttrval(xml_node, "unit");
      if ( unit == NULL )
      {
         SCIPerrorMessage("Flow unit of compressorStation <%s> not found.\n", id);
         return SCIP_READERROR;
      }
      value = SCIPxmlGetAttrval(xml_node, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Flow value of compressorStation <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      q = atof(value);

      SCIP_CALL( change_flow_unit(unit, &q, probdata->network->normDensity) );
   }

   /* find configuration and set binvars accordingly */
   xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "configuration", 0, 1);
   if ( xml_node == NULL )
   {
      SCIPerrorMessage("Could not find the configuration of compressorStation <%s>.\n", id);
      return SCIP_READERROR;
   }
   else
   {
      /* find config */
      conf_id = SCIPxmlGetAttrval(xml_node, "id"); /*lint !e838*/
      if ( conf_id == NULL )
      {
         SCIPerrorMessage("Attribute \"id\" of configuration of compressorStation <%s> not found.\n", id);
         return SCIP_READERROR;
      }
   }
   assert( conf_id != NULL );

   if ( strcmp("bypass", conf_id) == 0 )
   {
      /* set bypass binvar to 1 if bypass exists */
      if ( ! cs->internalBypassRequired )
      {
         SCIPerrorMessage("Cannot set compressorStation <%s> to bypass mode, because there is none.\n", id);
         return SCIP_INVALIDDATA;
      }

      bin = (int) SCIPcomputeVarUbGlobal(scip, probdata->VALVE_binvars[cs->bypassPosition]);
      if ( bin == 0 )
      {
         SCIPerrorMessage("Cannot set bypass_binvar of compressorStation <%s> to 1.\n", id);
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPchgVarLb(scip, probdata->VALVE_binvars[cs->bypassPosition], 1.0) );
      }

      /* fix flow variable */
      SCIP_CALL( SCIPfixVar(scip, cs->bypass->flowvar, q, &infeasible, &fixed) );
      assert( !infeasible );
      assert( fixed );
   }
   else
   {
      /* set bypass binvar to 0 if bypass exists */
      if ( cs->internalBypassRequired )
      {
         bin = (int) SCIPcomputeVarLbGlobal(scip, probdata->VALVE_binvars[cs->bypassPosition]);
         if ( bin == 1 )
         {
            SCIPerrorMessage("Cannot set bypass_binvar of compressorStation <%s> to zero.\n", id);
            return SCIP_INVALIDDATA;
         }
         else
         {
            SCIP_CALL( SCIPchgVarUb(scip, probdata->VALVE_binvars[cs->bypassPosition], 0.0) );
         }

         SCIP_CALL( SCIPfixVar(scip, cs->bypass->flowvar, 0.0, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
   }

   if ( (strcmp("closed", conf_id) == 0) || (strcmp("bypass", conf_id) == 0) )
   {
      /* set binvar to 0 */
      bin = (int) SCIPcomputeVarLbGlobal(scip, cs->compressor_binvar);
      if ( bin == 1 )
      {
         SCIPerrorMessage("Cannot set compressor_binvar of compressorStation <%s> to zero.\n", id);
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(scip, cs->compressor_binvar, 0.0) );
      }

      /* fix flow variable */
      SCIP_CALL( SCIPfixVar(scip, arc->flowvar, 0.0, &infeasible, &fixed) );
      assert( !infeasible );
      assert( fixed );
   }
   else
   {
      /* set compressor binvar to 1 */
      bin = (int) SCIPcomputeVarUbGlobal(scip, cs->compressor_binvar);
      if ( bin == 0 )
      {
         SCIPerrorMessage("Cannot set compressor_binvar of compressorStation <%s> to 1.\n", id);
         return SCIP_INVALIDDATA;
      }
      else
      {
         SCIP_CALL( SCIPchgVarLb(scip, cs->compressor_binvar, 1.0) );
      }

      /* fix flow variable */
      SCIP_CALL( SCIPfixVar(scip, arc->flowvar, q, &infeasible, &fixed) );
      assert( !infeasible );
      assert( fixed );

      /* try to find configuration matching conf_id */
      if ( probdata->boxConstraintModel )
      {
         GAS_CSBoxCons* boxcons;
         int i;

         for (i = 0; i < cs->numconfigurations; ++i)
         {
            boxcons = &(cs->boxcons[i]);

            if ( strcmp(boxcons->id, conf_id) == 0 )
               break;
         }

         if ( i < cs->numconfigurations )
         {
            /* set configuration binvar to 1 */
            bin = (int) SCIPcomputeVarUbGlobal(scip, boxcons->config_binvar); /*lint !e771 */
            if ( bin == 0 )
            {
               SCIPerrorMessage("Cannot set config_binvar for configuration <%s> of compressorStation <%s> to 1.\n", conf_id, id);
               return SCIP_INVALIDDATA;
            }
            else
            {
               SCIP_CALL( SCIPchgVarLb(scip, boxcons->config_binvar, 1.0) );
            }

            /* fix flow variable */
            SCIP_CALL( SCIPfixVar(scip, boxcons->BCM_flow, q, &infeasible, &fixed) );
            assert( !infeasible );
            assert( fixed );
         }
         else if ( strcmp("active", conf_id) == 0 )
         {
            /* cannot set a config binvar, if solution was produced with idealized compressor model */
         }
         else
         {
            SCIPerrorMessage("Could not find a configuration <%s> of compressorStation <%s>.\n", conf_id, id);
            return SCIP_INVALIDDATA;
         }
      }
   }

   /* set resistor variables if existing */
   /* inlet resistor */
   if ( (cs->NLRin_pressure != NULL) && (cs->NLRin_pressureDiff != NULL) )
   {
      xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "resistorIn", 0, 1);
      xml_child = NULL;
      if ( xml_node != NULL )
         xml_child = SCIPxmlFindNodeMaxdepth(xml_node, "pressureOut", 0, 1);

      if ( xml_child != NULL )
      {
         unit = SCIPxmlGetAttrval(xml_child, "unit");
         if ( unit == NULL )
         {
            SCIPerrorMessage("Unit of pressure for inlet resistor of compressorStation <%s> not found.\n", id);
            return SCIP_READERROR;
         }
         value = SCIPxmlGetAttrval(xml_child, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("Pressure value for inlet resistor of compressorStation <%s> not found.\n", id);
            return SCIP_READERROR;
         }

         p = atof(value);
         SCIP_CALL( change_unit_bar(unit, &p) );

         SCIP_CALL( SCIPfixVar(scip, cs->NLRin_pressure, p, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }

      xml_child = NULL;
      if ( xml_node != NULL )
         xml_child = SCIPxmlFindNodeMaxdepth(xml_node, "pressureDifferential", 0, 1);
      if ( xml_child != NULL )
      {
         unit = SCIPxmlGetAttrval(xml_child, "unit");
         if ( unit == NULL )
         {
            SCIPerrorMessage("Unit of pressureDifferential for inlet resistor of compressorStation <%s> not found.\n", id);
            return SCIP_READERROR;
         }
         value = SCIPxmlGetAttrval(xml_child, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("Pressure differential for inlet resistor of compressorStation <%s> not found.\n", id);
            return SCIP_READERROR;
         }

         pdiff = atof(value);
         SCIP_CALL( change_unit_bar(unit, &pdiff) );

         SCIP_CALL( SCIPfixVar(scip, cs->NLRin_pressureDiff, pdiff, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
   }

   /* outlet resistor */
   if ( cs->NLRout_pressure != NULL && cs->NLRout_pressureDiff != NULL )
   {
      xml_node = SCIPxmlFindNodeMaxdepth(xml_arc, "resistorOut", 0, 1);
      xml_child = NULL;
      if ( xml_node != NULL )
         xml_child = SCIPxmlFindNodeMaxdepth(xml_node, "pressureIn", 0, 1);

      if ( xml_child != NULL )
      {
         unit = SCIPxmlGetAttrval(xml_child, "unit");
         if ( unit == NULL )
         {
            SCIPerrorMessage("Unit of pressure for outlet resistor of compressorStation <%s> not found.\n", id);
            return SCIP_READERROR;
         }
         value = SCIPxmlGetAttrval(xml_child, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("Pressure value for outlet resistor of compressorStation <%s> not found.\n", id);
            return SCIP_READERROR;
         }

         p = atof(value);
         SCIP_CALL( change_unit_bar(unit, &p) );

         SCIP_CALL( SCIPfixVar(scip, cs->NLRout_pressure, p, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }

      xml_child = NULL;
      if ( xml_node != NULL )
         xml_child = SCIPxmlFindNodeMaxdepth(xml_node, "pressureDifferential", 0, 1);
      if ( xml_child != NULL )
      {
         unit = SCIPxmlGetAttrval(xml_child, "unit");
         if ( unit == NULL )
         {
            SCIPerrorMessage("Unit of pressureDifferential for outlet resistor of compressorStation <%s> not found.\n", id);
            return SCIP_READERROR;
         }
         value = SCIPxmlGetAttrval(xml_child, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("Pressure differential for outlet resistor of compressorStation <%s> not found.\n", id);
            return SCIP_READERROR;
         }

         pdiff = atof(value);
         SCIP_CALL( change_unit_bar(unit, &pdiff) );

         SCIP_CALL( SCIPfixVar(scip, cs->NLRout_pressureDiff, pdiff, &infeasible, &fixed) );
         assert( !infeasible );
         assert( fixed );
      }
   }

   return SCIP_OKAY;
}

/** read the node data in a solution file and set variables */
static
SCIP_RETCODE readNodeSolution(
   SCIP*                 scip,               /**< SCIP instance */
   const XML_NODE*       xml_node            /**< XML block for node */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_HASHTABLE* NodesTable;
   const XML_NODE* xml_pressure;
   const XML_NODE* xml_flow;
   GAS_Node* node;
   const char* id;
   const char* value;
   const char* unit;
   SCIP_Real p;
   SCIP_Real psquare;
   SCIP_Real q;

   assert( scip != NULL );
   assert( xml_node != NULL );

   probdata = SCIPgetProbData(scip);
   NodesTable = probdata->network->Nodeshashtable;

   /* find variable id */
   id = SCIPxmlGetAttrval(xml_node, "id");
   if ( id == NULL )
   {
      SCIPerrorMessage("Attribute \"id\" not found.\n");
      return SCIP_READERROR;
   }

   /* retrieve GAS_Node from hashtable */
   node = (GAS_Node*) SCIPhashtableRetrieve(NodesTable, (void*) id);
   if ( node == NULL )
   {
      SCIPerrorMessage("Could not find GAS_Node <%s> given in the solution file. Maybe solution file does not match network?\n", id);
      return SCIP_READERROR;
   }

   /* read pressure value at node and set the corresponding variable bounds */
   xml_pressure = SCIPxmlFindNodeMaxdepth(xml_node, "pressure", 0, 1);
   if ( xml_pressure == NULL )
   {
      SCIPerrorMessage("Could not find the pressure at node <%s>.\n", id);
      return SCIP_READERROR;
   }
   else
   {
      unit = SCIPxmlGetAttrval(xml_pressure, "unit");
      if ( unit == NULL )
      {
         SCIPerrorMessage("Unit of pressure at node <%s> not found.\n", id);
         return SCIP_READERROR;
      }
      value = SCIPxmlGetAttrval(xml_pressure, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Pressure value at node <%s> not found.\n", id);
         return SCIP_READERROR;
      }

      p = atof(value);
      SCIP_CALL( change_unit_bar(unit, &p) );

      if ( SCIPisFeasLE(scip, node->pressureMin,p) && SCIPisFeasGE(scip, node->pressureMax, p))
      {
         SCIP_CALL( SCIPchgVarLb(scip, node->pressurevar, p) );
         SCIP_CALL( SCIPchgVarUb(scip, node->pressurevar, p) );
      }
      else
      {
         SCIPerrorMessage("Pressure at node <%s> given in solution file does not satisfy the variable bounds!\n", id);
         return SCIP_INVALIDDATA;
      }
   }

   /* find and set squared pressure value */
   if ( node->needPIvar )
   {
      xml_pressure = SCIPxmlFindNodeMaxdepth(xml_node, "pressureSquare", 0, 1);
      if ( xml_pressure == NULL )
      {
         unit = SCIPxmlGetAttrval(xml_pressure, "unit");
         if ( unit == NULL )
         {
            SCIPerrorMessage("Unit of pressureSquare at node <%s> not found.\n", id);
            return SCIP_READERROR;
         }
         value = SCIPxmlGetAttrval(xml_pressure, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("Squared pressure value at node <%s> not found.\n", id);
            return SCIP_READERROR;
         }

         psquare = atof(value);

         if ( strcmp(unit, "bar_square") != 0 )
         {
            SCIPerrorMessage("Unknown unit for squared pressure variable at node <%s>.\n", id);
            return SCIP_INVALIDDATA;
         }

         if ( ! SCIPisEQ(scip, p, psquare) )
         {
            SCIPerrorMessage("Pressure does not match squared pressure at node <%s>.\n", id);
            return SCIP_INVALIDDATA;
         }

         if ( SCIPisFeasLE(scip, node->pressureMin * node->pressureMin, psquare) && SCIPisFeasGE(scip, node->pressureMax * node->pressureMax, psquare) )
         {
            SCIP_CALL( SCIPchgVarLb(scip, node->PIvar, psquare) );
            SCIP_CALL( SCIPchgVarUb(scip, node->PIvar, psquare) );
         }
         else
         {
            SCIPerrorMessage("Pressure at node <%s> given in solution file does not satisfy the variable bounds!\n", id);
            return SCIP_INVALIDDATA;
         }
      }
   }

   if ( node->type != INNODE )
   {
      xml_flow = SCIPxmlFindNodeMaxdepth(xml_node, "flow", 0, 1);
      if ( xml_flow == NULL )
      {
         SCIPerrorMessage("Could not find the in/outflow at node <%s>.\n", id);
         return SCIP_READERROR;
      }
      else
      {
         unit = SCIPxmlGetAttrval(xml_flow, "unit");
         if ( unit == NULL )
         {
            SCIPerrorMessage("Unit of flow at node <%s> not found.\n", id);
            return SCIP_READERROR;
         }
         value = SCIPxmlGetAttrval(xml_flow, "value");
         if ( value == NULL )
         {
            SCIPerrorMessage("In/Outflow value at node <%s> not found.\n", id);
            return SCIP_READERROR;
         }

         q = atof(value);
         SCIP_CALL( change_flow_unit(unit, &q, probdata->network->normDensity) );

         if ( node->type == EXIT )
            q *= - 1.0;

         if ( !SCIPisEQ(scip, node->flow, q) )
         {
            SCIPerrorMessage("In/Outflow at node <%s> provided in the solution file <%.15f> does not match the scenario <%.15f>.\n", id,q,node->flow);
            return SCIP_INVALIDDATA;
         }
      }
   }

   return SCIP_OKAY;
}

/** reads a given solution file and fixes variables to the given values */
SCIP_RETCODE readLSF(
   SCIP*                 scip,               /**< SCIP instance */
   const char*           pathToLsf           /**< path to LSF file */
   )
{  /*lint --e{668}*/
   SCIP_PROBDATA*  probdata;
   XML_NODE*       xml_start;                /* points to first xml_node in Gas file */
   const XML_NODE* xml_node_section;
   const XML_NODE* xml_node_subsection;
   const XML_NODE* xml_node;
   const XML_NODE* xml_arc_section;
   const XML_NODE* xml_arc;

   assert( scip != NULL );
   assert( pathToLsf != NULL );

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   if ( probdata->noFlowBinvars )
   {
      SCIPerrorMessage("Currently reading LS-files is not compatible with the model option noFlowBinvars.");
      return SCIP_ERROR;
   }

   /* read xml file */
   xml_start = SCIPxmlProcess(pathToLsf);

   if ( xml_start == NULL )
   {
      SCIPerrorMessage("Some error occured during parsing <%s>.\n", pathToLsf);
      return SCIP_READERROR;
   }

   /* find nodes and set the pressure values */
   xml_node_section = SCIPxmlFindNodeMaxdepth(xml_start, "nodes", 0, 2);
   xml_node_subsection = SCIPxmlFirstChild(xml_node_section); /* sources, sinks and innodes */

   while ( xml_node_subsection != NULL )
   {
      xml_node = SCIPxmlFirstChild(xml_node_subsection);
      while ( xml_node != NULL )
      {
         SCIP_CALL( readNodeSolution(scip, xml_node) );

         /* go to next node */
         xml_node = SCIPxmlNextSibl(xml_node);
      }
      xml_node_subsection = SCIPxmlNextSibl(xml_node_subsection);
   }

   /* find pipe section */
   xml_arc_section = SCIPxmlFindNodeMaxdepth(xml_start, "pipes", 0, 3);

   if (xml_arc_section != NULL )
      xml_arc = SCIPxmlFirstChild(xml_arc_section);
   else
      xml_arc = NULL;

   while ( xml_arc != NULL )
   {
      SCIP_CALL( readPipeSolution(scip, xml_arc) );
      xml_arc = SCIPxmlNextSibl(xml_arc);
   }

   /* find shortPipe section */
   xml_arc_section = SCIPxmlFindNodeMaxdepth(xml_start, "shortPipes", 0, 3);

   if (xml_arc_section != NULL )
      xml_arc = SCIPxmlFirstChild(xml_arc_section);
   else
      xml_arc = NULL;

   while ( xml_arc != NULL )
   {
      SCIP_CALL( readPipeSolution(scip, xml_arc) );
      xml_arc = SCIPxmlNextSibl(xml_arc);
   }

   /* find compressor station section */
   xml_arc_section = SCIPxmlFindNodeMaxdepth(xml_start, "compressorStations", 0, 3);

   if (xml_arc_section != NULL )
      xml_arc = SCIPxmlFirstChild(xml_arc_section);
   else
      xml_arc = NULL;

   while ( xml_arc != NULL )
   {
      SCIP_CALL( readCSsolution(scip, xml_arc) );
      xml_arc = SCIPxmlNextSibl(xml_arc);
   }

   /* find valve section */
   xml_arc_section = SCIPxmlFindNodeMaxdepth(xml_start, "valves", 0, 3);

   if (xml_arc_section != NULL )
      xml_arc = SCIPxmlFirstChild(xml_arc_section);
   else
      xml_arc = NULL;

   while ( xml_arc != NULL )
   {
      SCIP_CALL( readValveSolution(scip, xml_arc) );
      xml_arc = SCIPxmlNextSibl(xml_arc);
   }

   /* find controlvalve section */
   xml_arc_section = SCIPxmlFindNodeMaxdepth(xml_start, "controlValves", 0, 3);

   if (xml_arc_section != NULL )
      xml_arc = SCIPxmlFirstChild(xml_arc_section);
   else
      xml_arc = NULL;

   while ( xml_arc != NULL )
   {
      SCIP_CALL( readCvSolution(scip, xml_arc) );
      xml_arc = SCIPxmlNextSibl(xml_arc);
   }

   /* find resistor section */
   xml_arc_section = SCIPxmlFindNodeMaxdepth(xml_start, "resistors", 0, 3);

   if (xml_arc_section != NULL )
      xml_arc = SCIPxmlFirstChild(xml_arc_section);
   else
      xml_arc = NULL;

   while ( xml_arc != NULL )
   {
      SCIP_CALL( readResistorSolution(scip, xml_arc) );
      xml_arc = SCIPxmlNextSibl(xml_arc);
   }

   /* free xml data */
   SCIPxmlFreeNode(xml_start);

   return SCIP_OKAY;
}
