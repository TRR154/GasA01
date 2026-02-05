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

/**@file   probdata_gas.c
 * @brief  handling of data needed for solving stationary gas transport problems
 * @author Marc Pfetsch
 * @author Alexandra Stern
 * @author Oliver Habeck
 * @author Dennis Gabriel
 */

#include "probdata_gas.h"
#include "characteristics_gas.h"
#include "unitConversion_gas.h"
#include "myxmldefs.h"
#include <scip/misc.h>
#include <xml/xml.h>
#include <string.h>
#include <math.h>
#include <time.h>



/* remove the following define, if default units of numbers should be assumed - otherwise an error will be issued */
#define ERROR_ON_UNSPECIFIED_UNIT


/*------------------local functions------------------------------------- */

/** get problem name */
static
SCIP_RETCODE getProblemName(
   const char*           filename,           /**< filename */
   char*                 probname,           /**< returned problem name */
   int                   maxsize             /**< maximal size of probname */
   )
{
   int i = 0;
   int l;
   int j;

   /* first find end of string */
   while ( filename[i] != 0) /* find end of string */
      ++i;
   l = i;                /* end of string */

   /* if we found ".gz" */
   if ( l > 3 && filename[l-3] == '.' && filename[l-2] == 'g' && filename[l-1] == 'z' )
   {
      l = l - 4;
      i = l;
   }

   /* go back until '.' or '/' appears */
   while ((i > 0) && (filename[i] != '.') && (filename[i] != '/'))
      --i;
   assert( i > 0 );

   /* if we found '.', search for '/' */
   if (filename[i] == '.')
   {
      l = i;
      while ((i > 0) && (filename[i] != '/'))
	 --i;
   }

   /* correct counter */
   if (filename[i] == '/')
      ++i;

   /* copy name */
   j = 0;
   while ( (i < l) && (filename[i] != 0) )
   {
      probname[j++] = filename[i++];
      if ( j >= maxsize-1 )
         return SCIP_READERROR;
   }
   probname[j] = 0;

   return SCIP_OKAY;
}



/*------------------hash function------------------------------------- */

/** hash function used to handle the id of a node, gets the element as the key */
static
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyNode)
{  /*lint --e{715} */
   GAS_Node* node = (GAS_Node*) elem;
   return (void*) node->id;
}

static
SCIP_DECL_HASHGETKEY(SCIPhashGetKeyArc)
{  /*lint --e{715} */
   GAS_Arc* arc = (GAS_Arc*) elem;
   return (void*) arc->id;
}

/* ----------------- SCIP interface functions ------------------------ */

/** delete problem data */
static
SCIP_DECL_PROBDELORIG(probdelorigGAS)
{  /*lint --e{831} */
   int i;
   int j;
   int k;

   assert( probdata != NULL );
   assert( *probdata != NULL );
   assert( (*probdata)->network != NULL );
   assert( (*probdata)->flowvars != NULL || ! SCIPisTransformed(scip) );
   assert( (*probdata)->pressurevars != NULL || ! SCIPisTransformed(scip) );

   /* release all variables and free memory in reverse order of probdata */

   /*
    * release variables for compressor stations and free variable memory
    */
   if ( (*probdata)->network->numcompressor > 0 )
   {
      /* release all variables associated with the box constraint model */
      if ( (*probdata)->BCM_flowCS != NULL )
      {
         for (k = (*probdata)->network->numcsconfigurations - 1; k >= 0 ; --k)
         {
            if ( (*probdata)->BCM_pressureOut[k] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->BCM_pressureOut[k]) );
            }
            if ( (*probdata)->BCM_pressureIn[k] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->BCM_pressureIn[k]) );
            }
            if ( (*probdata)->BCM_flowCS[k] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->BCM_flowCS[k]) );
            }
         }
         SCIPfreeBlockMemoryArray(scip, &((*probdata)->BCM_pressureOut), (*probdata)->network->numcsconfigurations);
         SCIPfreeBlockMemoryArray(scip, &((*probdata)->BCM_pressureIn), (*probdata)->network->numcsconfigurations);
         SCIPfreeBlockMemoryArray(scip, &((*probdata)->BCM_flowCS), (*probdata)->network->numcsconfigurations);
      }

      /* release variables for resistors in compressor stations and free memory */
      if ( (*probdata)->CS_resistorvars != NULL && (*probdata)->network->numcsvarsfornlrs != 0 )
      {
         for (k = (*probdata)->network->numcsvarsfornlrs -1; k >= 0 ; --k)
         {
            if ( (*probdata)->CS_resistorvars[k] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->CS_resistorvars[k]) );
            }
         }
         SCIPfreeBlockMemoryArray(scip, &((*probdata)->CS_resistorvars), (*probdata)->network->numcsvarsfornlrs);
      }

      /* release binary variables */
      if ( (*probdata)->CS_binvars != NULL )
      {
         int numcsbinvars = (*probdata)->network->numcompressor;

         if ((*probdata)->boxConstraintModel )
         {
            numcsbinvars += (*probdata)->network->numcsconfigurations;
         }

         for (k = numcsbinvars - 1; k >= 0; --k)
         {
            if ( (*probdata)->CS_binvars[k] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->CS_binvars[k]) );
            }
         }
         SCIPfreeBlockMemoryArray(scip, &((*probdata)->CS_binvars), numcsbinvars);
      }
   }

   /*
    * release variables for linear resistors and free variable memory
    */
   if ( (*probdata)->LR_negFlowDir != NULL && (*probdata)->network->numlinresistor > 0 )
   {
      assert( (*probdata)->LR_smoothingFlow != NULL );
      assert( (*probdata)->LR_posFlowDir != NULL );
      for (k = (*probdata)->network->numlinresistor - 1; k >= 0 ; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->LR_negFlowDir[k]) );
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->LR_smoothingFlow[k]) );
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->LR_posFlowDir[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->LR_negFlowDir), (*probdata)->network->numlinresistor);
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->LR_smoothingFlow), (*probdata)->network->numlinresistor);
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->LR_posFlowDir), (*probdata)->network->numlinresistor);
   }

   /*
    * release variables for nonlinear resistors and free variable memory
    */
   if ( (*probdata)->NLR_AbsFlowVars != NULL && (*probdata)->network->numnonlinresistor > 0 )
   {
      assert( (*probdata)->NLR_AbsDeltaVars != NULL );
      assert( (*probdata)->NLR_DeltaVars != NULL );
      for (k = (*probdata)->network->numnonlinresistor - 1; k >= 0 ; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->NLR_AbsFlowVars[k]) );
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->NLR_AbsDeltaVars[k]) );
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->NLR_DeltaVars[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->NLR_AbsFlowVars), (*probdata)->network->numnonlinresistor);
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->NLR_AbsDeltaVars), (*probdata)->network->numnonlinresistor);
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->NLR_DeltaVars), (*probdata)->network->numnonlinresistor);
   }

   /*
    * release the binary varibales for controlvalves and free memory
    */
   if ( (*probdata)->CV_binvars != NULL && (*probdata)->network->numcontrolvalves > 0 )
   {
      for (k = (*probdata)->network->numcontrolvalves - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->CV_binvars[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->CV_binvars), (*probdata)->network->numcontrolvalves);
   }

   /*
    * release binary varibales for valves and free memory
    */
   if ( (*probdata)->VALVE_binvars != NULL && (*probdata)->network->numvalves > 0 )
   {
      for (k = (*probdata)->network->numvalves - 1; k >= 0 ; --k)
      {
         if ( (*probdata)->VALVE_binvars[k] != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->VALVE_binvars[k]) );
         }
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->VALVE_binvars), (*probdata)->network->numvalves);
   }

   /*
    * release PIdiffvars for the algebraic model and free memory
    */
   if ( (*probdata)->PIdiffvars != NULL )
   {
      for (k = (*probdata)->network->numpipes - 1; k >= 0; --k)
      {
         if ( (*probdata)->PIdiffvars[k] != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->PIdiffvars[k]) );
         }
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->PIdiffvars), (*probdata)->network->numpipes);
   }

   /*
    * release pSlackVars and free memory
    */
   if ( (*probdata)->pSlackVars != NULL )
   {
      for (k = (*probdata)->network->numnodes - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->pSlackVars[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->pSlackVars), (*probdata)->network->numnodes);
   }

   /*
    * release the squared pressure variables and free memory
    */
   if ( (*probdata)->PIvars != NULL )
   {
      for (k = (*probdata)->network->numpivars - 1; k >= 0; --k)
      {
         if ( (*probdata)->PIvars[k] != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->PIvars[k]) );
         }
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->PIvars), (*probdata)->network->numpivars);
   }

   /*
    * release the pressure variables and free memory
    */
   if ( (*probdata)->pressurevars != NULL )
   {
      for (k = (*probdata)->network->numnodes - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->pressurevars[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->pressurevars), (*probdata)->network->numnodes);
   }

   /*
    * release flow direction variables and free memory
    */
   if ( (*probdata)->positiveFlowBinvars != NULL )
   {
      assert( (*probdata)->negativeFlowBinvars != NULL );
      for (k = (*probdata)->network->numflowdirvars - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->negativeFlowBinvars[k]) );
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->positiveFlowBinvars[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->negativeFlowBinvars), (*probdata)->network->numflowdirvars);
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->positiveFlowBinvars), (*probdata)->network->numflowdirvars);
   }

   /*
    * release flow variables and free memory
    */
   if ( (*probdata)->flowvars != NULL )
   {
      for (k = (*probdata)->network->numarcs - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->flowvars[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->flowvars), (*probdata)->network->numarcs);
   }

   /*
    * release mixing variables and free memory
    */
   if ( (*probdata)->mixingRatioMass != NULL )
   {
      for (k = (*probdata)->network->numarcs - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->mixingRatioMass[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->mixingRatioMass), (*probdata)->network->numarcs);
   }

   if ( (*probdata)->mixingRatio != NULL )
   {
      for (k = (*probdata)->network->numarcs - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->mixingRatio[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->mixingRatio), (*probdata)->network->numarcs);
   }

   if ( (*probdata)->Lambda != NULL )
   {
      for (k = (*probdata)->network->numarcs - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->Lambda[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->Lambda), (*probdata)->network->numarcs);
   }

   if ( (*probdata)->nodeMixingRatio != NULL )
   {
      for (k = (*probdata)->network->numnodes - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->nodeMixingRatio[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->nodeMixingRatio), (*probdata)->network->numnodes);
   }

   if ( (*probdata)->nodeMixingRatioMol != NULL )
   {
      for (k = (*probdata)->network->numnodes - 1; k >= 0; --k)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->nodeMixingRatioMol[k]) );
      }
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->nodeMixingRatioMol), (*probdata)->network->numnodes);
   }

   /*
    * release netflow variables and free memory
    */
   if ( (*probdata)->netflowvars != NULL )
   {
      for (k = (*probdata)->network->numnodes - 1; k >= 0; --k)
      {
         if ( (*probdata)->netflowvars[k] != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->netflowvars[k]) );
         }
      }
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->netflowvars, (*probdata)->network->numnodes);
   }

   /*
    * release objective variable
    */
   if ( (*probdata)->objective != NULL )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->objective) );
   }

   /*
    * all variables are released and memory is freed
    * nodes and arc data is only freed if this is the original scip
    */
   if ( SCIPgetSubscipDepth(scip) > 0 )
   {
      /* free probdata */
      SCIPfreeBlockMemory(scip, probdata);

      return SCIP_OKAY;
   }

   if ( (*probdata)->scnname != NULL)
   {
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->scnname), strlen((*probdata)->scnname) + 1);
   }

   if ( (*probdata)->netname != NULL)
   {
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->netname), strlen((*probdata)->netname) + 1);
   }

   /*
    * free data of cycles and spanning tree
    */
   while ( (*probdata)->network->lastCycle != NULL )
   {
      if ( (*probdata)->network->lastCycle->arcsInCycle != NULL )
      {
         SCIPfreeBlockMemoryArray(scip, &((*probdata)->network->lastCycle->arcsInCycle), (*probdata)->network->numarcs);
      }
      if ( (*probdata)->network->lastCycle->combinedBasisCycles != NULL )
      {
         SCIPfreeBlockMemoryArray(scip, &((*probdata)->network->lastCycle->combinedBasisCycles), (*probdata)->network->numBasisCycles);
      }

      if ( (*probdata)->network->lastCycle->previous_cycle != NULL )
      {
         (*probdata)->network->lastCycle =  (*probdata)->network->lastCycle->previous_cycle;
         SCIPfreeBlockMemory(scip, &((*probdata)->network->lastCycle->next_cycle) );
      }
      else
      {
         SCIPfreeBlockMemory(scip, &((*probdata)->network->lastCycle) );
         (*probdata)->network->lastCycle = NULL;
      }
   }

   if ( (*probdata)->network->spanningTree != NULL )
   {
      /* free spanningTree data */
      if ( (*probdata)->network->spanningTree->arcsInTree != NULL )
      {
         SCIPfreeBlockMemoryArray(scip, &((*probdata)->network->spanningTree->arcsInTree), (*probdata)->network->numarcs);
      }

      if ( (*probdata)->network->spanningTree->wayToNode != NULL )
      {
         for (k = (*probdata)->network->numarcs - 1; k >=0; --k)
         {
            SCIPfreeBlockMemoryArray(scip, &((*probdata)->network->spanningTree->wayToNode[k]), (*probdata)->network->numnodes);
         }
         SCIPfreeBlockMemoryArray(scip, &((*probdata)->network->spanningTree->wayToNode), (*probdata)->network->numarcs);
      }
      SCIPfreeBlockMemory(scip, &((*probdata)->network->spanningTree) );
   }

   /*
    * free data from csfile:
    * note that this is not really in reverse order of allocation,
    * since the order in csfile is not always in order of netfile.
    * Also pistonCS and TurboCS need not appear in this order.
    */
   if ( (*probdata)->boxConstraintModel || (*probdata)->charDiagram )
   {
      GAS_CS* cs_struct;

      for (i = (*probdata)->network->numarcs -1; i >= 0; --i)
      {
         if ( (*probdata)->network->arcs_ptr[i].type != CS )
            continue;

         cs_struct = (*probdata)->network->arcs_ptr[i].detailed_info;

         if ( cs_struct->configurations != NULL )
         {
            GAS_CSConfig* config;
            for (k = cs_struct->numconfigurations - 1; k >= 0; --k)
            {
               config = &(cs_struct->configurations[k]);
               if ( config->stageOneCsTwo != NULL )
               {
                  SCIPfreeBlockMemoryArray(scip, &(config->stageOneCsTwo), strlen(config->stageOneCsTwo) + 1);
               }
               if ( config->stageOneCsOne != NULL )
               {
                  SCIPfreeBlockMemoryArray(scip, &(config->stageOneCsOne), strlen(config->stageOneCsOne) + 1);
               }

               if ( config->id != NULL )
               {
                  SCIPfreeBlockMemoryArray(scip, &(config->id), strlen(config->id) + 1);
               }
            }
            SCIPfreeBlockMemoryArray(scip, &(cs_struct->configurations), cs_struct->numconfigurations);
         }

         if ( cs_struct->pistoncs != NULL )
         {
            GAS_PistonCS* pistoncs;
            for (k = cs_struct->numPistonCS - 1; k >= 0; --k)
            {
               pistoncs = &(cs_struct->pistoncs[k]);
               SCIPfreeBlockMemoryArray(scip, &(pistoncs->id), strlen(pistoncs->id) + 1);
            }
            SCIPfreeBlockMemoryArray(scip, &(cs_struct->pistoncs), cs_struct->numPistonCS);
         }

         if ( cs_struct->turbocs != NULL )
         {
            GAS_TurboCS* turbocs;
            for (k = cs_struct->numTurboCS - 1; k >= 0; --k)
            {
               turbocs = &(cs_struct->turbocs[k]);
               SCIPfreeBlockMemoryArray(scip, &(turbocs->id), strlen(turbocs->id) + 1);
            }
            SCIPfreeBlockMemoryArray(scip, &(cs_struct->turbocs), cs_struct->numTurboCS);
         }

         if ( cs_struct->boxcons != NULL )
         {
            GAS_CSBoxCons* boxcons;
            for (k = cs_struct->numconfigurations - 1; k >= 0; --k)
            {
               boxcons = &(cs_struct->boxcons[k]);

               if ( boxcons->c_unit != NULL )
               {
                  SCIPfreeBlockMemoryArray(scip, &(boxcons->c_unit), strlen(boxcons->c_unit) + 1);
                  SCIPfreeBlockMemoryArray(scip, &(boxcons->c), strlen(boxcons->c) + 1);
                  SCIPfreeBlockMemoryArray(scip, &(boxcons->b_unit), strlen(boxcons->b_unit) + 1);
                  SCIPfreeBlockMemoryArray(scip, &(boxcons->b), strlen(boxcons->b) + 1);
                  SCIPfreeBlockMemoryArray(scip, &(boxcons->a_unit), strlen(boxcons->a_unit) + 1);
                  SCIPfreeBlockMemoryArray(scip, &(boxcons->a), strlen(boxcons->a) + 1);
                  SCIPfreeBlockMemoryArray(scip, &(boxcons->rhs), boxcons->numfacets);
               }

               if ( boxcons->facetcoeff != NULL)
               {
                  for (j = boxcons->numfacets - 1; j >= 0 ; --j)
                  {
                     SCIPfreeBlockMemoryArray(scip, &(boxcons->facetcoeff[j]), 3);
                  }
                  SCIPfreeBlockMemoryArray(scip, &(boxcons->facetcoeff), boxcons->numfacets);
               }

               SCIPfreeBlockMemoryArray(scip, &(boxcons->id), strlen(boxcons->id) + 1);
            }
            SCIPfreeBlockMemoryArray(scip, &(cs_struct->boxcons), cs_struct->numconfigurations);
         }
      }
   }

   /*
    * free arc data
    */
   for (k = (*probdata)->network->numarcs - 1; k >= 0; --k)
   {
      if ( (*probdata)->network->arcs_ptr[k].detailed_info != NULL )
      {
         switch ( (*probdata)->network->arcs_ptr[k].type )
         {
         case PIPE:
         {
            GAS_Pipe* pipe = (*probdata)->network->arcs_ptr[k].detailed_info;
            SCIPfreeBlockMemory(scip, &pipe);
            break;
         }
         case VALVE:
         {
            GAS_Valve* valve = (*probdata)->network->arcs_ptr[k].detailed_info;
            SCIPfreeBlockMemory(scip, &valve);
            break;
         }
         case CONTROLVALVE:
         {
            GAS_Controlvalve* cv = (*probdata)->network->arcs_ptr[k].detailed_info;
            SCIPfreeBlockMemory(scip, &cv);
            break;
         }
         case CS:
         {
            GAS_CS* cs = (*probdata)->network->arcs_ptr[k].detailed_info;
            SCIPfreeBlockMemory(scip, &cs);
            break;
         }
         case RESISTOR:
         {
            GAS_Resistor* resistor = (*probdata)->network->arcs_ptr[k].detailed_info;
            SCIPfreeBlockMemory(scip, &resistor);
            break;
         }
         case SHORTPIPE:
            /* nothting to do here */
            break;
         case UNKNOWNARCTYPE:
         default:
            SCIPerrorMessage("Wrong arc type %d.\n", (*probdata)->network->arcs_ptr[k].type);
            return SCIP_ERROR;
         }
         (*probdata)->network->arcs_ptr[k].detailed_info = NULL; /* for debugging */
      }
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->network->arcs_ptr[k].id, strlen((*probdata)->network->arcs_ptr[k].id) + 1);
   }
   SCIPfreeBlockMemoryArray(scip, &((*probdata)->network->arcs_ptr), (*probdata)->network->numarcs);

   /*
    * free node data
    */
   for (k = (*probdata)->network->numnodes - 1; k >= 0; --k)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->network->nodes_ptr[k].id, strlen((*probdata)->network->nodes_ptr[k].id) + 1);
   }
   SCIPfreeBlockMemoryArray(scip, &((*probdata)->network->nodes_ptr), (*probdata)->network->numnodes);

   /*
    * free network
    */
   SCIPfreeBlockMemory(scip, &((*probdata)->network));

   /*
    * free probdata
    */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}


/** copies user data of source SCIP for the target SCIP */
static
SCIP_DECL_PROBCOPY(probcopyGAS)
{
   SCIP_Bool valid = TRUE;
   SCIP_VAR* var;
   int i;
#if SCIP_VERSION < 900
   SCIP_Bool original = TRUE;

   if ( SCIPisTransformed(sourcescip) )
      original = FALSE;
#endif

   assert( scip != NULL );
   assert( sourcescip != NULL );
   assert( sourcedata != NULL );
   assert( targetdata != NULL );

   /* set up data */
   SCIP_CALL( SCIPallocBlockMemory(scip, targetdata) );

   /* fill in data and allocate memory in the order given in probdata */
   (*targetdata)->network = sourcedata->network; /* Be aware that variable pointer in gas arcs and nodes still point to original variables! */
   (*targetdata)->netname = NULL;
   (*targetdata)->scnname = NULL;

   if ( sourcedata->objective != NULL )
   {
      if ( ! original )
      {
         SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->objective, &var) );
      }
      else
         var = sourcedata->objective;

      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->objective), varmap, consmap, global, &valid) );
      SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->objective) );
      assert( valid );
   }
   else
   {
      (*targetdata)->objective = NULL;
   }

   if ( sourcedata->network->numarcs > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->flowvars), sourcedata->network->numarcs) );

      for ( i = 0; i < sourcedata->network->numarcs; ++i )
      {
         /* get transformed flowvars */
         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar( sourcescip, sourcedata->flowvars[i], &var) );
         }
         else
            var = sourcedata->flowvars[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->flowvars[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->flowvars[i]) );
         assert( valid );
      }

      if ( sourcedata->noFlowBinvars )
      {
         (*targetdata)->positiveFlowBinvars = NULL;
         (*targetdata)->negativeFlowBinvars = NULL;
      }
      else
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray( scip, &( (*targetdata)->positiveFlowBinvars), sourcedata->network->numflowdirvars ) );
         SCIP_CALL( SCIPallocClearBlockMemoryArray( scip, &( (*targetdata)->negativeFlowBinvars), sourcedata->network->numflowdirvars ) );

         for ( i = 0; i < sourcedata->network->numflowdirvars; ++i )
         {
            /* get transformed flow binvars */
            if ( sourcedata->positiveFlowBinvars[i] != NULL )
            {
               if ( ! original )
               {
                  SCIP_CALL( SCIPgetTransformedVar( sourcescip, sourcedata->positiveFlowBinvars[i], &var ) );
               }
               else
                  var = sourcedata->positiveFlowBinvars[i];

               SCIP_CALL( SCIPgetVarCopy( sourcescip, scip, var, &( ( *targetdata )->positiveFlowBinvars[i] ), varmap, consmap, global, &valid ) );
               SCIP_CALL( SCIPcaptureVar( scip, ( *targetdata )->positiveFlowBinvars[i] ) );
               assert( valid );
            }
            else
            {
               ( *targetdata )->positiveFlowBinvars[i] = NULL;
            }

            if ( sourcedata->negativeFlowBinvars[i] != NULL )
            {
               if ( ! original )
               {
                  SCIP_CALL( SCIPgetTransformedVar( sourcescip, sourcedata->negativeFlowBinvars[i], &var ) );
               }
               else
                  var = sourcedata->negativeFlowBinvars[i];

               SCIP_CALL( SCIPgetVarCopy( sourcescip, scip, var, &( ( *targetdata )->negativeFlowBinvars[i] ), varmap, consmap, global, &valid ) );
               SCIP_CALL( SCIPcaptureVar( scip, ( *targetdata )->negativeFlowBinvars[i] ) );
               assert( valid );
            }
            else
            {
               ( *targetdata )->negativeFlowBinvars[i] = NULL;
            }
         }
      }
   }

   if ( sourcedata->maxFlow )
   {
      assert( sourcedata->netflowvars != NULL );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->netflowvars), sourcedata->network->numnodes) );
      for (i = 0; i < sourcedata->network->numnodes; i++)
      {
         if ( sourcedata->netflowvars[i] != NULL )
         {
            if ( ! original )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->netflowvars[i], &var) );
            }
            else
               var = sourcedata->netflowvars[i];

            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->netflowvars[i]), varmap, consmap, global, &valid) );
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->netflowvars[i]) );
            assert( valid );
         }
      }
   }
   else
   {
      (*targetdata)->netflowvars = NULL;
   }

   if ( sourcedata->network->numnodes > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->pressurevars), sourcedata->network->numnodes) );
      for (i = 0; i < sourcedata->network->numnodes; i++)
      {
         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->pressurevars[i], &var) );
         }
         else
            var = sourcedata->pressurevars[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->pressurevars[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->pressurevars[i]) );
         assert( valid );
      }
   }
   else
   {
      (*targetdata)->pressurevars = NULL;
   }

   if ( sourcedata->PIvars != NULL )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->PIvars), sourcedata->network->numpivars) );

      for (i = 0; i < sourcedata->network->numpivars; i++)
      {
         if ( sourcedata->PIvars[i] != NULL )
         {
            if ( ! original )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->PIvars[i], &var) );
            }
            else
               var = sourcedata->PIvars[i];

            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->PIvars[i]), varmap, consmap, global, &valid) );
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->PIvars[i]) );
            assert( valid );
         }
      }
   }
   else
   {
      (*targetdata)->PIvars = NULL;
   }

   if ( sourcedata->minSlackPerBound && sourcedata->pSlackVars != NULL )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->pSlackVars), sourcedata->network->numnodes) );

      for (i = 0; i < sourcedata->network->numnodes; ++i)
      {
         if ( sourcedata->pSlackVars[i] != NULL )
         {
            if ( ! original )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->pSlackVars[i], &var) );
            }
            else
               var = sourcedata->pSlackVars[i];

            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->pSlackVars[i]), varmap, consmap, global, &valid) );
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->pSlackVars[i]) );
            assert( valid );
         }
      }
   }
   else
   {
      (*targetdata)->pSlackVars = NULL;
   }

   if ( sourcedata->network->numpipes > 0 && sourcedata->PIdiffvars != NULL )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(*targetdata)->PIdiffvars, sourcedata->network->numpipes) );
      for (i = 0; i < sourcedata->network->numpipes; ++i)
      {
         if ( sourcedata->PIdiffvars[i] != NULL )
         {
            if ( ! original )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->PIdiffvars[i], &var) );
            }
            else
               var = sourcedata->PIdiffvars[i];

            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->PIdiffvars[i]), varmap, consmap, global, &valid) );
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->PIdiffvars[i]) );
            assert( valid );
         }
      }
   }
   else
      (*targetdata)->PIdiffvars = NULL;

   if ( sourcedata->network->numvalves > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->VALVE_binvars), sourcedata->network->numvalves) );
      for (i = 0; i < sourcedata->network->numvalves; ++i)
      {
         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->VALVE_binvars[i], &var) );
         }
         else
            var = sourcedata->VALVE_binvars[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->VALVE_binvars[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->VALVE_binvars[i]) );
         assert( valid );
      }
   }
   else
      (*targetdata)->VALVE_binvars = NULL;

   if ( sourcedata->network->numcontrolvalves > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->CV_binvars), sourcedata->network->numcontrolvalves) );
      for (i = 0; i < sourcedata->network->numcontrolvalves; ++i)
      {
         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->CV_binvars[i], &var) );
         }
         else
            var = sourcedata->CV_binvars[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->CV_binvars[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->CV_binvars[i]) );
         assert( valid );
      }
   }
   else
      (*targetdata)->CV_binvars = NULL;

   if ( sourcedata->network->numnonlinresistor > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->NLR_DeltaVars), sourcedata->network->numnonlinresistor) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->NLR_AbsDeltaVars), sourcedata->network->numnonlinresistor) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->NLR_AbsFlowVars), sourcedata->network->numnonlinresistor) );

      for (i = 0; i < sourcedata->network->numnonlinresistor; ++i)
      {
         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->NLR_DeltaVars[i], &var) );
         }
         else
            var = sourcedata->NLR_DeltaVars[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->NLR_DeltaVars[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->NLR_DeltaVars[i]) );
         assert( valid );

         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->NLR_AbsDeltaVars[i], &var) );
         }
         else
            var = sourcedata->NLR_AbsDeltaVars[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->NLR_AbsDeltaVars[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->NLR_AbsDeltaVars[i]) );
         assert( valid );

         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->NLR_AbsFlowVars[i], &var) );
         }
         else
            var = sourcedata->NLR_AbsFlowVars[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->NLR_AbsFlowVars[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->NLR_AbsFlowVars[i]) );
         assert( valid );
      }
   }
   else
   {
      (*targetdata)->NLR_DeltaVars = NULL;
      (*targetdata)->NLR_AbsDeltaVars = NULL;
      (*targetdata)->NLR_AbsFlowVars = NULL;
   }

   if ( sourcedata->network->numlinresistor > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->LR_posFlowDir), sourcedata->network->numlinresistor) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->LR_smoothingFlow), sourcedata->network->numlinresistor) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->LR_negFlowDir), sourcedata->network->numlinresistor) );

      for (i = 0; i < sourcedata->network->numlinresistor; ++i)
      {
         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->LR_posFlowDir[i], &var) );
         }
         else
            var = sourcedata->LR_posFlowDir[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->LR_posFlowDir[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->LR_posFlowDir[i]) );
         assert( valid );

         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->LR_smoothingFlow[i], &var) );
         }
         else
            var = sourcedata->LR_smoothingFlow[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->LR_smoothingFlow[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->LR_smoothingFlow[i]) );
         assert( valid );

         if ( ! original )
         {
            SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->LR_negFlowDir[i], &var) );
         }
         else
            var = sourcedata->LR_negFlowDir[i];

         SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->LR_negFlowDir[i]), varmap, consmap, global, &valid) );
         SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->LR_negFlowDir[i]) );
         assert( valid );
      }
   }
   else
   {
      (*targetdata)->LR_posFlowDir    = NULL;
      (*targetdata)->LR_negFlowDir    = NULL;
      (*targetdata)->LR_smoothingFlow = NULL;
   }

   if ( sourcedata->network->numcompressor > 0 && sourcedata->CS_binvars != NULL )
   {
      int numcsbinvars = sourcedata->network->numcompressor;
      if ( sourcedata->boxConstraintModel )
      {
         numcsbinvars += sourcedata->network->numcsconfigurations;
      }
      /* transform binary variables for compressors */
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->CS_binvars), numcsbinvars) );
      for (i = 0; i < numcsbinvars; ++i)
      {
         if ( sourcedata->CS_binvars[i] != NULL )
         {
            if ( ! original )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->CS_binvars[i], &var) );
            }
            else
               var = sourcedata->CS_binvars[i];

            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->CS_binvars[i]), varmap, consmap, global, &valid) );
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->CS_binvars[i]) );
            assert( valid );
         }
      }
   }
   else
   {
      (*targetdata)->CS_binvars = NULL;
   }

   /* transform variables for nonlinear resistors in compressor stations */
   if ( sourcedata->network->numcsvarsfornlrs > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->CS_resistorvars), sourcedata->network->numcsvarsfornlrs) );
      for (i = 0; i < sourcedata->network->numcsvarsfornlrs; ++i)
      {
         if ( sourcedata->CS_resistorvars[i] != NULL )
         {
            if ( ! original )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->CS_resistorvars[i], &var) );
            }
            else
               var = sourcedata->CS_resistorvars[i];

            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->CS_resistorvars[i]), varmap, consmap, global, &valid) );
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->CS_resistorvars[i]) );
            assert( valid );
         }
      }
   }
   else
   {
      (*targetdata)->CS_resistorvars = NULL;
   }

   /* transform variables for the boxconstraint model */
   if ( sourcedata->network->numcompressor > 0 && sourcedata->BCM_flowCS != NULL )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->BCM_flowCS), sourcedata->network->numcsconfigurations) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->BCM_pressureIn), sourcedata->network->numcsconfigurations) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &((*targetdata)->BCM_pressureOut), sourcedata->network->numcsconfigurations) );

      for (i = 0; i < sourcedata->network->numcsconfigurations; ++i)
      {
         if ( sourcedata->BCM_flowCS[i] != NULL )
         {
            if ( ! original )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->BCM_flowCS[i], &var) );
            }
            else
               var = sourcedata->BCM_flowCS[i];

            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->BCM_flowCS[i]), varmap, consmap, global, &valid) );
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->BCM_flowCS[i]) );
            assert( valid );
         }

         if ( sourcedata->BCM_pressureIn[i] != NULL )
         {
            if ( ! original )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->BCM_pressureIn[i], &var) );
            }
            else
               var = sourcedata->BCM_pressureIn[i];

            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->BCM_pressureIn[i]), varmap, consmap, global, &valid) );
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->BCM_pressureIn[i]) );
            assert( valid );
         }
         if ( sourcedata->BCM_pressureOut[i] != NULL )
         {
            if ( ! original )
            {
               SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcedata->BCM_pressureOut[i], &var) );
            }
            else
               var = sourcedata->BCM_pressureOut[i];

            SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, var, &((*targetdata)->BCM_pressureOut[i]), varmap, consmap, global, &valid) );
            SCIP_CALL( SCIPcaptureVar(scip, (*targetdata)->BCM_pressureOut[i]) );
            assert( valid );
         }
      }
   }
   else
   {
      (*targetdata)->BCM_flowCS = NULL;
      (*targetdata)->BCM_pressureIn = NULL;
      (*targetdata)->BCM_pressureOut = NULL;
   }

   /* set remaining data */
   (*targetdata)->noFlowBinvars      = sourcedata->noFlowBinvars;
   (*targetdata)->binaryFlowCons     = sourcedata->binaryFlowCons;
   (*targetdata)->noCycles           = sourcedata->noCycles;
   (*targetdata)->allCycles          = sourcedata->allCycles;
   (*targetdata)->algebraic          = sourcedata->algebraic;
   (*targetdata)->papay              = sourcedata->papay;
   (*targetdata)->boxConstraintModel = sourcedata->boxConstraintModel;
   (*targetdata)->additionalFacets   = sourcedata->additionalFacets;
   (*targetdata)->charDiagram        = sourcedata->charDiagram;
   (*targetdata)->relaxLowerBounds   = sourcedata->relaxLowerBounds;
   (*targetdata)->relaxUpperBounds   = sourcedata->relaxUpperBounds;
   (*targetdata)->relaxCVbounds      = sourcedata->relaxCVbounds;
   (*targetdata)->relaxLBvalue       = sourcedata->relaxLBvalue;
   (*targetdata)->relaxUBvalue       = sourcedata->relaxUBvalue;
   (*targetdata)->relaxCVvalue       = sourcedata->relaxCVvalue;
   (*targetdata)->resistorCVInternal = sourcedata->resistorCVInternal;
   (*targetdata)->minCompSum         = sourcedata->minCompSum;
   (*targetdata)->minPressSum        = sourcedata->minPressSum;
   (*targetdata)->noObjective        = sourcedata->noObjective;
   (*targetdata)->minLambda        = sourcedata->minLambda;
   (*targetdata)->powerLoss          = sourcedata->powerLoss;
   (*targetdata)->minSlack           = sourcedata->minSlack;
   (*targetdata)->minSlackPerBound   = sourcedata->minSlackPerBound;
   (*targetdata)->idealCSminIncrease = sourcedata->idealCSminIncrease;
   (*targetdata)->idealCSmaxIncrease = sourcedata->idealCSmaxIncrease;
   (*targetdata)->csusemaxincrease   = sourcedata->csusemaxincrease;
   (*targetdata)->maxFlow            = sourcedata->maxFlow;
   (*targetdata)->addComponentCuts   = sourcedata->addComponentCuts;
   (*targetdata)->addParallelCSCuts  = sourcedata->addParallelCSCuts;
   (*targetdata)->reduceFlowConservation = sourcedata->reduceFlowConservation;

   /* currently do not copy mixing variables */
   (*targetdata)->mixingRatioMass = NULL;
   (*targetdata)->Lambda = NULL;
   (*targetdata)->mixingRatio = NULL;
   (*targetdata)->mixingRatioMass = NULL;
   (*targetdata)->nodeMixingRatio = NULL;
   (*targetdata)->nodeMixingRatioMol = NULL;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/* ----------------- auxiliary functions ------------------------ */

/** function to bind an arc as an outgoing arc to a node*/
static
void insert_outarc(
   GAS_Node*             G,                  /**< pointer to GAS_Node */
   GAS_Arc*              P                   /**< pointer to GAS_Arc */
   )
{
   if ( G->outarcs == NULL )
   {
      G->outarcs = P;
      G->outarcs->next_outarc = NULL;
   }
   else
   {
      GAS_Arc* ptr;
      ptr = (G->outarcs);
      G->outarcs = P;
      G->outarcs->next_outarc = ptr;
   }
}

/** function to bind an arc as an ingoing arc to a node */
static
void insert_inarc(
   GAS_Node*             G,                  /**< pointer to GAS_Node */
   GAS_Arc*              P                   /**< pointer to GAS_Arc */
   )
{
   if ( G->inarcs == NULL )
   {
      G->inarcs = P;
      G->inarcs->next_inarc = NULL;
   }
   else
   {
      GAS_Arc* ptr;
      ptr = (G->inarcs);
      G->inarcs = P;
      G->inarcs->next_inarc = ptr;
   }
}

/** takes string of arc types and returns corresponding enum */
static
Arc_type convertToArcType(
   const char*           name                /**< string with arc type */
   )
{
   if ( strcmp(name, "pipe") == 0 )
      return PIPE;
   else if ( strcmp(name, "shortPipe") == 0 )
      return SHORTPIPE;
   else if ( strcmp(name, "valve") == 0 )
      return VALVE;
   else if ( strcmp(name, "controlValve") == 0 )
      return CONTROLVALVE;
   else if ( strcmp(name, "compressorStation") == 0 )
      return CS;
   else if ( strcmp(name, "resistor") == 0 )
      return RESISTOR;

   return UNKNOWNARCTYPE;
}

/** takes string of node types and returns corresponding enum */
static
Node_type convertToNodeType(
   const char*           name                /**< string with node type */
   )
{
   if ( strcmp(name, "exit") == 0 )
      return EXIT;
   else if ( strcmp(name, "entry") == 0 )
      return ENTRY;
   else if ( strcmp(name, "innode") == 0 )
      return INNODE;

   return UNKNOWNNODETYPE;
}

/** read node data (gaslib format) */
static
SCIP_RETCODE readNodeData(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network,            /**< network of the node */
   const XML_NODE*       source,             /**< XML node to read */
   int                   nodeposition,       /**< position in node array */
   int*                  foundGasConstants   /**< counter for number times we found the gas constants */
   )
{
   const XML_NODE* source_elem;    /* source element in xml file */
   size_t          attrnamelen;    /* length of attribute name */
   SCIP_Real*      valu;           /* pointer to a real number */
   SCIP_Real       convalue;       /* real number */
   const char*     attrname;       /* attribute name of xml node */
   const char*     attrvalue;      /* attribute value of xml node */
   GAS_Node*       N;              /* gas node */
   int height      = 0;            /* ints for checking if elements found */
   int pressureMax = 0;
   int pressureMin = 0;
   int normDensity = 0;
   int pcPressure  = 0;
   int pcTemp      = 0;
   int gasTemp     = 0;
   int molarMass   = 0;
   int specificHeatCapacity = 0;

   /* find variable id */
   attrname = SCIPxmlGetAttrval(source, "id");
   if ( attrname == NULL )
   {
      SCIPerrorMessage("Attribute \"id\" not found.\n");
      return SCIP_READERROR;
   }

   /* check whether node is already present in hashtable */
   if ( SCIPhashtableRetrieve(network->Nodeshashtable, (void*) attrname) )
   {
      SCIPerrorMessage("Node <%s> already exists!\n", attrname);
      return SCIP_READERROR;
   }

   /* determine current node */
   N = &(network->nodes_ptr[nodeposition]);
   N->outarcs = NULL;
   N->inarcs = NULL;

   /* set default values */
   N->nodeposition    = nodeposition;
   N->geoWGS84Long    = 0.0;
   N->geoWGS84Lat     = 0.0;
   N->x               = 0.0;
   N->y               = 0.0;
   N->flow            = 0.0;
   N->scnFlowMax      = 0.0;
   N->scnFlowMin      = 0.0;
   N->pressureMin     = 1.01325;
   N->pressureMax     = 200.0;
   N->netPressureMin  = 1.01325;
   N->netPressureMax  = 200.0;
   N->scnPressureMin  = 1.01325;
   N->scnPressureMax  = 200.0;
   N->arcPressureMax  = 200.0;
   N->numneighbors    = 0;
   N->numincidentarcs = 0;
   N->pressurevar     = NULL;
   N->needPIvar       = FALSE;
   N->PIvar           = NULL;
   N->pSlackVar       = NULL;
   N->type            = INNODE;  /* we initialize all nodes as innodes; when reading the SCN-file, this type gets changed. */
   N->molarMass       = -1;      /* initalized as -1 to check later if node is an entry*/
   N->specificHeatCap = -1;      /* initalized as -1 to check later if node is an entry*/

   /* copy id */
   attrnamelen = strlen(attrname);
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(N->id), attrname, attrnamelen + 1) );

   /* read geographical information */
   attrvalue = SCIPxmlGetAttrval(source, "geoWGS84Long");
   if ( attrvalue == NULL )
      SCIPdebugMessage("Attribute \"geoWGS84Long\" not found.\n");
   else
      N->geoWGS84Long = atof(attrvalue);

   /* read geographical information */
   attrvalue = SCIPxmlGetAttrval(source, "geoWGS84Lat");
   if ( attrvalue == NULL )
      SCIPdebugMessage("Attribute \"geoWGS84Lat\" not found.\n");
   else
      N->geoWGS84Lat = atof(attrvalue);

   /* read coordinate x */
   attrvalue = SCIPxmlGetAttrval(source, "x");
   if ( attrvalue == NULL )
      SCIPdebugMessage("Attribute \"x\" not found.\n");
   else
      N->x = atof(attrvalue);

   /* read coordinate y */
   attrvalue = SCIPxmlGetAttrval(source, "y");
   if ( attrvalue == NULL )
      SCIPdebugMessage("Attribute \"y\" not found.\n");
   else
      N->y = atof(attrvalue);

   /* find variable height, need loop through the elements of source because the information is saved in a child node of source */
   for (source_elem = SCIPxmlFirstChild(source); source_elem != NULL; source_elem = SCIPxmlNextSibl(source_elem)) /*lint !e2840*/
   {
      const char* elemvalue;
      const char* elemunit;

      /* find element height*/
      if ( strcmp(SCIPxmlGetName(source_elem), "height") == 0)
      {
         ++height;
         elemvalue = SCIPxmlGetAttrval(source_elem, "value");
         if ( elemvalue == NULL )
         {
            SCIPerrorMessage(" \"Height\" not found.\n");
            return SCIP_READERROR;
         }
         else
            N->height = atof(elemvalue);
      }
      /* find element pressureMin */
      else if ( strcmp(SCIPxmlGetName(source_elem), "pressureMin") == 0 )
      {
         ++pressureMin;
         elemvalue = SCIPxmlGetAttrval(source_elem, "value");
         if ( elemvalue == NULL )
         {
            SCIPerrorMessage("Element \"pressureMin\" of node %s not found.\n", N->id);
            return SCIP_READERROR;
         }
         else
         {
            N->netPressureMin = atof(elemvalue);

            /* detemine the unit of the pressure*/
            elemunit = SCIPxmlGetAttrval(source_elem, "unit");
            if ( elemunit == NULL )
            {
               SCIPerrorMessage("Unit of pressureMin not found.\n");
               return SCIP_READERROR;
            }
            /* change unit, if necessary */
            valu = &(N->netPressureMin);
            SCIP_CALL( change_unit_bar(elemunit, valu) );
         }
      }
      /* find element pressureMax */
      else if ( strcmp(SCIPxmlGetName(source_elem), "pressureMax") == 0 )
      {
         ++pressureMax;
         elemvalue = SCIPxmlGetAttrval(source_elem, "value");
         if ( elemvalue == NULL )
         {
            SCIPerrorMessage("Element \"pressureMax\" of node %s not found.\n", N->id);
            return SCIP_READERROR;
         }
         else
         {
            N->netPressureMax = atof(elemvalue);

            /* detemine the unit of the pressure */
            elemunit = SCIPxmlGetAttrval(source_elem, "unit");
            if (elemunit == NULL )
            {
               SCIPerrorMessage("Unit of pressureMax not found.\n");
               return SCIP_READERROR;
            }
            /* change unit, if necessary */
            valu = &(N->netPressureMax);
            SCIP_CALL( change_unit_bar(elemunit, valu) );
         }
      }
      /* if the xml file contains more than one value for normDensity we use the first found value */
      else if ( strcmp(SCIPxmlGetName(source_elem), "normDensity") == 0 )
      {
         ++normDensity;
         elemvalue = SCIPxmlGetAttrval(source_elem, "value");

         /* determine the unit */
         elemunit = SCIPxmlGetAttrval(source_elem, "unit");
         if ( elemunit == NULL )
         {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
            SCIPerrorMessage("Unit of normDensity not found.\n");
            return SCIP_READERROR;
#else
            SCIPwarningMessage(scip, "Unit of normDensity not found. Assuming 'kg_per_m_cube'.\n");
            convalue = atof(elemvalue);
#endif
         }
         else
         {
            /* change unit, if necessary */
            convalue = atof(elemvalue);
            SCIP_CALL( change_density_unit(elemunit, &convalue) );
         }

         if ( *foundGasConstants == 0 )
         {
            network->normDensity = convalue;
            SCIPdebugMessage("NormDensity: <%f> \n", network->normDensity);
         }
      }
      /* reads pseudocritical pressure */
      else if ( strcmp(SCIPxmlGetName(source_elem), "pseudocriticalPressure") == 0 )
      {
         ++pcPressure;
         elemvalue = SCIPxmlGetAttrval(source_elem, "value");
         elemunit = SCIPxmlGetAttrval(source_elem, "unit");
         if ( elemunit == NULL )
         {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
            SCIPerrorMessage("Unit of pseudocriticalPressure not found.\n");
            return SCIP_READERROR;
#else
            SCIPwarningMessage(scip, "Unit of pseudocriticalPressure not found. Assuming 'bar'.\n");
            convalue = atof(elemvalue);
#endif
         }
         else
         {
            /* change unit, if necessary */
            convalue = atof(elemvalue);
            SCIP_CALL( change_unit_bar(elemunit, &convalue) );
         }

         if ( *foundGasConstants == 0 )
         {
            network->pseudocriticalPressure = convalue;
            SCIPdebugMessage("pseudocriticalPressure: <%f> \n", network->pseudocriticalPressure);
         }
      }
      /* reads pseudocritical temperature */
      else if ( strcmp(SCIPxmlGetName(source_elem), "pseudocriticalTemperature") == 0 )
      {
         ++pcTemp;
         elemvalue = SCIPxmlGetAttrval(source_elem, "value");
         elemunit = SCIPxmlGetAttrval(source_elem, "unit");
         if ( elemunit == NULL )
         {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
            SCIPerrorMessage("Unit of pseudocriticalTemperature not found.\n");
            return SCIP_READERROR;
#else
            SCIPwarningMessage(scip, "Unit of pseudocriticalTemperature not found. Assuming 'K'.\n");
            convalue = atof(elemvalue);
#endif
         }
         else
         {
            /* change unit, if necessary */
            convalue = atof(elemvalue);
            SCIP_CALL( change_temperature_unit(elemunit, &convalue) );
         }

         if ( *foundGasConstants == 0 )
         {
            network->pseudocriticalTemperature = convalue;
            SCIPdebugMessage("pseudocriticalTemperature: <%f> \n", network->pseudocriticalTemperature);
         }
      }
      /* reads gas temperature */
      else if ( strcmp(SCIPxmlGetName(source_elem), "gasTemperature") == 0 )
      {
         ++gasTemp;
         elemvalue = SCIPxmlGetAttrval(source_elem, "value");
         elemunit = SCIPxmlGetAttrval(source_elem, "unit");
         if ( elemunit == NULL )
         {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
            SCIPerrorMessage("Unit of gasTemperature not found.\n");
            return SCIP_READERROR;
#else
            SCIPwarningMessage(scip, "Unit of gasTemperature not found. Assuming 'K'.\n");
            convalue = atof(elemvalue);
#endif
         }
         else
         {
            /* change unit, if necessary */
            convalue = atof(elemvalue);
            SCIP_CALL( change_temperature_unit(elemunit, &convalue) );
         }

         if ( *foundGasConstants == 0 )
         {
            network->gasTemperature = convalue;
            SCIPdebugMessage("gasTemperature: <%f> \n", network->gasTemperature);
         }
      }
      /* If molar mass is spefified in the SCN file, the molar mass set here will be overwritten and replaced by SCN molar mass */
      else if ( strcmp(SCIPxmlGetName(source_elem), "molarMass") == 0 )
      {
         ++molarMass;
         elemvalue = SCIPxmlGetAttrval(source_elem, "value");
         elemunit = SCIPxmlGetAttrval(source_elem, "unit");
         if ( elemunit == NULL )
         {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
            SCIPerrorMessage("Unit of molarMass not found.\n");
            return SCIP_READERROR;
#else
            SCIPwarningMessage(scip, "Unit of molarMass not found. Assuming 'kg_per_kmol'.\n");
            convalue = atof(elemvalue);
#endif
         }
         else
         {
            /* change unit, if necessary */
            convalue = atof(elemvalue);
            SCIP_CALL( change_molarmass_unit(elemunit, &convalue) );
            N->molarMass = convalue;
         }

         if ( *foundGasConstants == 0 )
         {
            network->molarMass1 = convalue;
            SCIPdebugMessage("molarMass: <%f> \n", convalue);
         }
      }
      /* specificHeatCapacity is set */
      else if ( strcmp(SCIPxmlGetName(source_elem), "specificHeatCapacity") == 0 )
      {
         ++specificHeatCapacity;
         elemvalue = SCIPxmlGetAttrval(source_elem, "value");
         convalue = atof(elemvalue);
         N->specificHeatCap = convalue;
      }
   }

   /* check if a property of the node not found or found more than once.
    * necessary properties should appear once and optional properties at most once */
   if ( height != 1 || pressureMin != 1 || pressureMax != 1 || normDensity > 1 || pcPressure > 1
      || pcTemp > 1 || gasTemp > 1 || molarMass > 1 || specificHeatCapacity > 1 )
   {
      SCIPerrorMessage("Invalid data of node %s.\n", N->id);
      return SCIP_READERROR;
   }

   /* if all optional properties were found, increase the associated counter.
    * else if some but not all were found, return an error */
   if ( normDensity == 1 && pcPressure == 1 && pcTemp == 1 && gasTemp == 1 && molarMass == 1 )
      ++(*foundGasConstants);
   else if ( normDensity != 0 || pcPressure != 0 || pcTemp != 0 || gasTemp != 0 || molarMass != 0 )
   {
      SCIPerrorMessage("Only parts of the gas constants were present in node %s.\n", N->id);
      return SCIP_READERROR;
   }

   /* now insert node into hashtable*/
   SCIP_CALL( SCIPhashtableInsert(network->Nodeshashtable, (void*) N) );

   /* xmlFreeNode(source_elem); */
   return SCIP_OKAY;
}

/** reads the data in attribute flowMax */
static
SCIP_RETCODE readFlowMax(
   const XML_NODE*       arc_elem,           /**< XML Node to read */
   GAS_Arc*              gas_arc,            /**< arc to store data */
   SCIP_Real             normDensity         /**< norm density */
   )
{
   const char* p_elemvalue;
   const char* p_elemunit;

   assert( strcmp(SCIPxmlGetName(arc_elem), "flowMax") == 0 );
   p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
   if ( p_elemvalue == NULL )
   {
      SCIPdebugMessage("Value of element \"flowMax\" not found.\n");
      gas_arc->flowMax = SCIP_DEFAULT_INFINITY;
   }
   else
   {
      gas_arc->flowMax = atof(p_elemvalue);

      /* reads unit */
      p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
      if ( p_elemunit == NULL )
      {
         SCIPerrorMessage(" \"Unit\" of \"flowMax\" not found.\n");
         return SCIP_READERROR;
      }

      /* changes flow unit, if necessary */
      SCIP_CALL( change_flow_unit(p_elemunit, &(gas_arc->flowMax), normDensity) );
   }
   return SCIP_OKAY;
}

/** reads the data in attribute flowMin */
static
SCIP_RETCODE readFlowMin(
   const XML_NODE*       arc_elem,           /**< XML Node to read */
   GAS_Arc*              gas_arc,            /**< arc to store data */
   SCIP_Real             normDensity         /**< norm density */
   )
{
   const char* p_elemvalue;
   const char* p_elemunit;

   assert( strcmp(SCIPxmlGetName(arc_elem), "flowMin") == 0 );
   p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
   if ( p_elemvalue == NULL )
   {
      SCIPdebugMessage("Value of element \"flowMin\" not found.\n");
      gas_arc->flowMin = -SCIP_DEFAULT_INFINITY;
   }
   else
   {
      gas_arc->flowMin = atof(p_elemvalue);

      /* reads unit */
      p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
      if ( p_elemunit == NULL )
      {
         SCIPerrorMessage(" \"Unit\" of \"flowMin\" not found.\n");
         return SCIP_READERROR;
      }

      /* changes flow unit, if necessary */
      SCIP_CALL( change_flow_unit(p_elemunit, &(gas_arc->flowMin), normDensity) );
   }
   return SCIP_OKAY;
}

/** reads resistor specific data */
static
SCIP_RETCODE readResistorData(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       xml_arc,            /**< xml node to read */
   GAS_Arc*              gas_arc,            /**< arc to complete */
   GAS_Network*          network             /**< network */
   )
{
   const XML_NODE* arc_elem;
   GAS_Resistor* resistor;

   assert( scip != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, &resistor) );
   gas_arc->detailed_info = resistor;
   gas_arc->type = RESISTOR;

   resistor->pressureLoss = 0.0;
   resistor->dragFactor = 0.0;
   resistor->diameter = 1.0;
   resistor->LR_position = -1;
   resistor->NLR_position = -1;

   for (arc_elem = SCIPxmlFirstChild(xml_arc); arc_elem != NULL; arc_elem = SCIPxmlNextSibl(arc_elem)) /*lint !e2840*/
   {
      const char* p_elemvalue;
      const char* p_elemunit;

      /* find element drag factor */
      if ( strcmp(SCIPxmlGetName(arc_elem), "dragFactor") == 0)
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"dragFactor\" not found.\n");
            return SCIP_READERROR;
         }
         else
            resistor->dragFactor = atof(p_elemvalue);
      }
      /* find element pressureLoss */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureLoss") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureLoss\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            resistor->pressureLoss = atof(p_elemvalue);
            /* determine the unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
               SCIPerrorMessage("Unit of \"pressureLoss\" not found.\n");
               return SCIP_READERROR;
#else
               SCIPwarningMessage(scip, "Unit of \"pressureLoss\" not found. Assuming 'bar'.\n");
#endif
            }
            else
            {
               /* change unit, if necessary */
               SCIP_CALL( change_unit_bar(p_elemunit, &resistor->pressureLoss) );
            }
         }
      }
      /* find element diameter */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "diameter") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"diameter\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            resistor->diameter = atof(p_elemvalue);
            /* determine the unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" Unit of \"diameter\" not found.\n");
               return SCIP_READERROR;
            }
            /* change unit, if necessary */
            SCIP_CALL( change_unit_meter(p_elemunit, &resistor->diameter) );
         }
      }
      /* find element flowMin by using function readFlowMin */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMin") == 0 )
      {
         SCIP_CALL( readFlowMin(arc_elem, gas_arc, network->normDensity) );
      }
      /* find element flowMax by using function readFlowMax */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMax") == 0 )
      {
         SCIP_CALL( readFlowMax(arc_elem, gas_arc, network->normDensity) );
      }
   }

   /* decide whether the resistor is linear or nonlinear;
      Linear: dragFactor = 0.0;
      Nonlinear: pressureLoss = 0.0 */
   if ( resistor->dragFactor == 0.0 )
   {
      resistor->type = LINEAR;
      resistor->LR_position = network->numlinresistor;
      ++(network->numlinresistor);
   }
   else if ( resistor->pressureLoss == 0.0 )
   {
      resistor->type = NONLINEAR;
      resistor->NLR_position = network->numnonlinresistor;
      if ( ! gas_arc->sourcenode->needPIvar )
      {
         gas_arc->sourcenode->needPIvar = TRUE;
         ++(network->numpivars);
      }
      if ( ! gas_arc->targetnode->needPIvar )
      {
         gas_arc->targetnode->needPIvar = TRUE;
         ++(network->numpivars);
      }
      ++(network->numnonlinresistor);
   }
   else
   {
      SCIPerrorMessage("Failed to decide if resistor is linear or nonlinear.\n");
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}


/** read controlvalve specific data */
static
SCIP_RETCODE readControlvalveData(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       xml_arc,            /**< xml node to read */
   GAS_Arc*              gas_arc,            /**< arc to complete */
   GAS_Network*          network,            /**< network data */
   SCIP_Bool*            create_bypass       /**< boolean wheter or not a bypass should be added */
   )
{
   const XML_NODE* arc_elem;
   GAS_Controlvalve* controlvalve;
   const char* bypass;

   assert( scip != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, &controlvalve) );
   gas_arc->detailed_info = controlvalve;
   gas_arc->type = CONTROLVALVE;

   /* initialize gas valve */
   controlvalve->pressureLossIn          = 0.0;
   controlvalve->pressureLossOut         = 0.0;
   controlvalve->dragFactorIn            = 0.0;
   controlvalve->dragFactorOut           = 0.0;
   controlvalve->pressureInMin           = 1.0;
   controlvalve->pressureOutMax          = 120.0;
   controlvalve->pressureDifferentialMin = 0.0;
   controlvalve->pressureDifferentialMax = 120.0;
   controlvalve->binvar                  = NULL;
   controlvalve->internalBypassRequired  = 1;
   controlvalve->bypass                  = NULL;
   *create_bypass                        = TRUE;

   /* read attribute internalBypassRequired */
   bypass = SCIPxmlGetAttrval(xml_arc, "internalBypassRequired");
   if ( bypass == NULL )
   {
      SCIPdebugMessage("Use default value \"internalBypassRequired = 1\" for controlvalve.\n");
      controlvalve->internalBypassRequired = 1;
   }
   else
   {
      controlvalve->internalBypassRequired = atoi(bypass);
      if ( ! controlvalve->internalBypassRequired )
         *create_bypass = FALSE;
   }

   for (arc_elem = SCIPxmlFirstChild(xml_arc); arc_elem != NULL; arc_elem = SCIPxmlNextSibl(arc_elem)) /*lint !e2840*/
   {
      const char* p_elemvalue;
      const char* p_elemunit;

      /* reads element pressureLossIn */
      if ( strcmp(SCIPxmlGetName(arc_elem), "pressureLossIn") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureLossIn\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            controlvalve->pressureLossIn = atof(p_elemvalue);

            /* determine unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" Unit of \"pressureLossIn\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(p_elemunit, &controlvalve->pressureLossIn) );
         }
      }
      /* reads element pressureLossOut */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureLossOut") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureLossOut\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            controlvalve->pressureLossOut = atof(p_elemvalue);

            /* determine unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" Unit of \"pressureLossOut\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(p_elemunit, &controlvalve->pressureLossOut) );
         }
      }
      /* reads element pressureInMin */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureInMin") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureInMin\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            controlvalve->pressureInMin = atof(p_elemvalue);

            /* determine unit */
            p_elemunit=SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" Unit of \"pressureInMin\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(p_elemunit, &controlvalve->pressureInMin) );
         }
      }
      /* reads element pressureOutMax */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureOutMax") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureOutMax\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            controlvalve->pressureOutMax = atof(p_elemvalue);

            /* determine unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" \"Unit\" of \"pressureOutMax\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(p_elemunit, &controlvalve->pressureOutMax) );
         }
      }
      /* reads element pressureDifferentialMin */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureDifferentialMin") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureDiffferentialMin\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            controlvalve->pressureDifferentialMin = atof(p_elemvalue);

            /* determine unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
               SCIPerrorMessage(" \"Unit\" of \"pressureDifferentialMin\" not found.\n");
               return SCIP_READERROR;
#else
               SCIPwarningMessage(scip, " \"Unit\" of \"pressureDifferentialMin\" not found. Assuming 'bar'.\n");
#endif
            }
            else
            {
               SCIP_CALL( change_unit_bar(p_elemunit, &controlvalve->pressureDifferentialMin) );
            }
         }
      }
      /* find element pressureDifferentialMax */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureDifferentialMax") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureDifferentialMax\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            controlvalve->pressureDifferentialMax = atof(p_elemvalue);

            /* determine unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
               SCIPerrorMessage(" \"Unit\" of \"pressureDifferentialMax\" not found.\n");
               return SCIP_READERROR;
#else
               SCIPwarningMessage(scip, " \"Unit\" of \"pressureDifferentialMax\" not found. Assuming 'bar'.\n");
#endif
            }
            else
            {
               SCIP_CALL( change_unit_bar(p_elemunit, &controlvalve->pressureDifferentialMax) );
            }
         }
      }
      /* find element pressureSet */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureSet") == 0 )
      {
         SCIPerrorMessage("Model for control valves with preset pressure not implemented.\n");
         return SCIP_READERROR;
      }
      /* find element flowMin using function readFlowMin */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMin") == 0 )
      {
         SCIP_CALL( readFlowMin(arc_elem, gas_arc, network->normDensity) );
      }
      /* find element flowMax using function readFlowMax */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMax") == 0 )
      {
         SCIP_CALL( readFlowMax(arc_elem, gas_arc, network->normDensity) );
         if ( gas_arc->flowMax <= 0.0 )
         {
            SCIPerrorMessage("Maximal flow on controlvalve %s has to be positive!\n", gas_arc->id);
            return SCIP_READERROR;
         }
      }
      /* find element dragFactorIn */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "dragFactorIn") == 0)
      {
         p_elemvalue= SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"dragFactorIn\" not found.\n");
            return SCIP_READERROR;
         }
         else
            controlvalve->dragFactorIn = atof(p_elemvalue);

         if ( REALABS(controlvalve->dragFactorIn) < 1e-9 )
            SCIPinfoMessage(scip, NULL, "Found zero DragFactorIn for control valve <%s>, leading to 0 resistor.\n", gas_arc->id);
         else
            SCIPwarningMessage(scip, "DragFactorIn for resistor in control valve <%s> is currently ignored.\n", gas_arc->id);
      }
      /* find element dragFactorOut */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "dragFactorOut") == 0)
      {
         p_elemvalue= SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"dragFactorOut\" not found.\n");
            return SCIP_READERROR;
         }
         else
            controlvalve->dragFactorOut = atof(p_elemvalue);

         if ( REALABS(controlvalve->dragFactorOut) < 1e-9 )
            SCIPinfoMessage(scip, NULL, "Found zero DragFactorOut for control valve <%s>, leading to 0 resistor.\n", gas_arc->id);
         else
            SCIPwarningMessage(scip, "DragFactorOut for resistor in control valve <%s> is currently ignored.\n", gas_arc->id);
      }
   }

   /* flow on controlvalve station can only be negative if cv is in bypass mode */
   if ( controlvalve->internalBypassRequired == 0 )
   {
      gas_arc->flowMin = MAX(gas_arc->flowMin, 0.0);
   }

   return SCIP_OKAY;
}


/** read compressorstation specific data in the net file */
static
SCIP_RETCODE readCSData(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       xml_arc,            /**< xml node to read */
   GAS_Arc*              gas_arc,            /**< arc to complete */
   GAS_Network*          network,            /**< network */
   SCIP_Bool*            create_bypass,      /**< boolean wheter or not a bypass should be added  */
   SCIP_Bool             boxConstraintModel  /**< wheter or not the box constraint model is used */
   )
{
   const XML_NODE* arc_elem;
   const char* elemvalue;
   GAS_CS* cs;

   SCIP_CALL( SCIPallocBlockMemory(scip, &cs) );
   gas_arc->detailed_info = cs;
   gas_arc->type = CS;

   /* initialize gas_cs */
   cs->pressureLossIn         = 0.0;
   cs->pressureLossOut        = 0.0;
   cs->dragFactorIn           = 0.0;
   cs->dragFactorOut          = 0.0;
   cs->diameterIn             = 0.0;
   cs->diameterOut            = 0.0;
   cs->pressureInMin          = 1.0;
   cs->pressureOutMax         = 120.0;
   cs->numPistonCS            = 0;
   cs->pistoncs               = NULL;
   cs->numTurboCS             = 0;
   cs->turbocs                = NULL;
   cs->numconfigurations      = 0;
   cs->configurations         = NULL;
   cs->boxcons                = NULL;
   cs->compressor_binvar      = NULL;
   cs->NLRin_pressure         = NULL;
   cs->NLRin_pressureDiff     = NULL;
   cs->NLRout_pressure        = NULL;
   cs->NLRout_pressureDiff    = NULL;

   /* set default values */
   cs->internalBypassRequired = 0;
   *create_bypass = FALSE;
   cs->bypassPosition = -1;
   cs->bypass = NULL;

   /* The value internalBypassRequired tells us if we have to model a bypass in the compressor model.
    * The default value is 1. */
   elemvalue = SCIPxmlGetAttrval(xml_arc, "internalBypassRequired");
   if ( elemvalue == NULL  || atoi(elemvalue) == 1 )
   {
      cs->internalBypassRequired = 1;
      *create_bypass = TRUE;
      cs->bypassPosition = network->numvalves;
   }
   else if ( atoi(elemvalue) != 0 )
   {
      SCIPerrorMessage("Unknown value of \"internalBypassRequired\" of compressor %s.\n", gas_arc->id);
      return SCIP_READERROR;
   }

   for (arc_elem = SCIPxmlFirstChild(xml_arc); arc_elem != NULL; arc_elem = SCIPxmlNextSibl(arc_elem)) /*lint !e2840*/
   {
      const char* p_elemvalue;
      const char* p_elemunit;

      /* find element pressureLossIn */
      if ( strcmp(SCIPxmlGetName(arc_elem), "pressureLossIn") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureLossIn\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            cs->pressureLossIn = atof(p_elemvalue);

            /* determine unit*/
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage("Unit of \"pressureLossIn\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(p_elemunit, &cs->pressureLossIn) );
         }
      }
      /* find element pressureLossOut */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureLossOut") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureLossOut\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            cs->pressureLossOut = atof(p_elemvalue);

            /* determines unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage("Unit of \"pressureLossOut\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(p_elemunit, &cs->pressureLossOut) );
         }
      }
      /* find element pressureInMin */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureInMin") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureInMin\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            cs->pressureInMin = atof(p_elemvalue);

            /* determines unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage("Unit of \"pressureInMin\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(p_elemunit, &cs->pressureInMin) );
         }
      }
      /* find element pressureOutMax */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureOutMax") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureOutMax\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            cs->pressureOutMax = atof(p_elemvalue);

            /* determines unit */
            p_elemunit=SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage("Unit of \"pressureOutMax\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(p_elemunit, &cs->pressureOutMax) );
         }
      }
      /* find element flowMin */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMin") == 0 )
      {
         SCIP_CALL( readFlowMin(arc_elem, gas_arc, network->normDensity) );
      }
      /* find element flowMax */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMax") == 0 )
      {
         SCIP_CALL( readFlowMax(arc_elem, gas_arc, network->normDensity) );
         if ( gas_arc->flowMax < 0.0 )
         {
            SCIPerrorMessage("Negative maximal flow on compressor %s is not possible.\n", gas_arc->id);
            return SCIP_READERROR;
         }
      }
      /* find element dragFactorIn */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "dragFactorIn") == 0)
      {
         p_elemvalue= SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"dragFactorIn\" not found.\n");
            return SCIP_READERROR;
         }
         else
            cs->dragFactorIn = atof(p_elemvalue);
      }
      /* find element dragFactorOut */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "dragFactorOut") == 0)
      {
         p_elemvalue= SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"dragFactorOut\" not found.\n");
            return SCIP_READERROR;
         }
         else
            cs->dragFactorOut = atof(p_elemvalue);
      }
      /* find element diameterIn */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "diameterIn") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"diameterIn\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            cs->diameterIn = atof(p_elemvalue);

            /* determine the unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" Unit of \"diameterIn\" not found.\n");
               return SCIP_READERROR;
            }
            /* change unit, if necessary */
            SCIP_CALL( change_unit_meter(p_elemunit, &cs->diameterIn) );
         }
      }
      /* find element diameterOut */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "diameterOut") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"diameterOut\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            cs->diameterOut = atof(p_elemvalue);

            /* determine the unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" Unit of \"diameterOut\" not found.\n");
               return SCIP_READERROR;
            }
            /* change unit, if necessary */
            SCIP_CALL( change_unit_meter(p_elemunit, &cs->diameterOut) );
         }
      }
   }

   /* flow on compressor station can only be negative if cs is in bypass mode */
   if ( cs->internalBypassRequired == 0 )
   {
      gas_arc->flowMin = MAX(gas_arc->flowMin, 0.0);
   }

   /* cannot have both linear and nonlinear resistors at entry and exit */
   if ( ! SCIPisZero(scip, cs->pressureLossIn) && !SCIPisZero(scip, cs->dragFactorIn) )
   {
      SCIPerrorMessage("Cannot have a nonlinear and a linear resistor at the entry of a compressor!\n");
      return SCIP_READERROR;
   }

   if ( ! SCIPisZero(scip, cs->pressureLossOut) && !SCIPisZero(scip, cs->dragFactorOut) )
   {
      SCIPerrorMessage("Cannot have a nonlinear and a linear resistor at the exit of a compressor!\n");
      return SCIP_READERROR;
   }

   /* set the correct number of additional variables needed for the nonlinear resistors */
   if ( boxConstraintModel )
   {
      if ( ! SCIPisZero(scip, cs->dragFactorIn) )
         network->numcsvarsfornlrs += 2;

      if ( ! SCIPisZero(scip, cs->dragFactorOut) )
         network->numcsvarsfornlrs += 2;
   }

   return SCIP_OKAY;
}


/** read bypass data from a compressor or a controlvalve arc */
static
SCIP_RETCODE readBypassData(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network,            /**< network of the node */
   GAS_Arc*              arc,                /**< pointer to the compressor or controlvalve arc to read from */
   int                   arcposition         /**< position in the arc array */
   )
{
   char             bypass_name[SCIP_MAXSTRLEN];
   GAS_Arc*         gas_arc;
   GAS_Valve*       valve;

   /* determine current arc*/
   gas_arc = &(network->arcs_ptr[arcposition]);
   gas_arc->next_inarc = NULL;
   gas_arc->next_outarc = NULL;
   gas_arc->arcposition = arcposition;
   gas_arc->detailed_info = NULL;
   gas_arc->flowMax = 0.0;
   gas_arc->flowMin = 0.0;
   gas_arc->type = VALVE;

   /* create link of cv to bypass */
   if ( arc->type == CONTROLVALVE )
   {
      GAS_Controlvalve* cv;
      cv = arc->detailed_info;
      cv->bypass = gas_arc;
   }
   else if ( arc->type == CS )
   {
      GAS_CS* cs;
      cs = arc->detailed_info;
      cs->bypass = gas_arc;
   }

   /* create name "id_bypass" with id from the arc */
   (void) SCIPsnprintf(bypass_name, SCIP_MAXSTRLEN, "%s_bypass", arc->id);
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(gas_arc->id), bypass_name, strlen(bypass_name)+1) );

   /* set the right source and target node */
   gas_arc->sourcenode = arc->sourcenode;
   gas_arc->targetnode = arc->targetnode;

   /* add the bypass to the in- and outarcs */
   insert_outarc(arc->sourcenode, gas_arc);
   insert_inarc(arc->targetnode, gas_arc);

   /* update the numincidentarcs counter */
   ++(gas_arc->sourcenode->numincidentarcs);
   ++(gas_arc->targetnode->numincidentarcs);

   /* determine valve specific data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &valve) );
   gas_arc->detailed_info = valve;

   /* initialize gas valve */
   valve->pressureDifferentialMax = 120.0;                  /* default value, because a compressor has no such data */
   valve->valve_binvar = NULL;

   if ( arc->type == CONTROLVALVE )
   {
      valve->bypassed_cv = arc;
      valve->bypassed_cs = NULL;
   }
   else if ( arc->type == CS )
   {
      valve->bypassed_cv = NULL;
      valve->bypassed_cs = arc;
   }

   /* set flowMin and flowMax */
   gas_arc->flowMin = MIN(arc->flowMin, 0.0);  /* flow 0.0 has to be possible such that valve can be closed */
   gas_arc->flowMax = arc->flowMax;

   /* negative flow is only possible over the bypass! */
   arc->flowMin = MAX(0.0, arc->flowMin);

   SCIP_CALL( SCIPhashtableInsert(network->Arcshashtable, (void*) gas_arc) );

   return SCIP_OKAY;
}

/** read valve specific data */
static
SCIP_RETCODE readValveData(
   SCIP*                 scip,               /**< SCIP instance */
   const XML_NODE*       xml_arc,            /**< xml node to read */
   GAS_Arc*              gas_arc,            /**< arc to complete */
   GAS_Network*          network             /**< network */
   )
{
   const XML_NODE* arc_elem;
   GAS_Valve* valve;

   SCIP_CALL( SCIPallocBlockMemory(scip, &valve) );
   gas_arc->detailed_info = valve;
   gas_arc->type = VALVE;

   /* initialize gas valve */
   valve->pressureDifferentialMax = 120.0;
   valve->valve_binvar = NULL;
   valve->bypassed_cv = NULL;
   valve->bypassed_cs = NULL;

   for (arc_elem = SCIPxmlFirstChild(xml_arc); arc_elem != NULL; arc_elem = SCIPxmlNextSibl(arc_elem)) /*lint !e2840*/
   {
      const char* p_elemvalue;
      const char* p_elemunit;

      /* reads element pressureDifferentialMax */
      if ( strcmp(SCIPxmlGetName(arc_elem), "pressureDifferentialMax") == 0 )
      {
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value for element \"pressureDifferentialMax\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            valve->pressureDifferentialMax = atof(p_elemvalue);

            /* determines unit */
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" \"Unit\" of \"pressureDifferentialMax\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(p_elemunit, &valve->pressureDifferentialMax) );
         }
      }
      /* reads element flowMin */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMin") == 0 )
      {
         SCIP_CALL( readFlowMin(arc_elem, gas_arc, network->normDensity) );
      }
      /* reads element flowMax */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMax") == 0 )
      {
         SCIP_CALL( readFlowMax(arc_elem, gas_arc, network->normDensity) );
      }
   }

   return SCIP_OKAY;
}


/** read pipe specific data */
static
SCIP_RETCODE readPipeData(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       xml_arc,            /**< xml node to read */
   GAS_Arc*              gas_arc,            /**< arc to complete */
   GAS_Network*          network             /**< network */
   )
{
   const XML_NODE* arc_elem;
   GAS_Pipe*       pipe;
   int             length      = 0;
   int             diameter    = 0;
   int             roughness   = 0;
   int             flowMin     = 0;
   int             flowMax     = 0;
   int             pressureMax = 0;

   SCIP_CALL( SCIPallocBlockMemory(scip, &pipe) );
   gas_arc->detailed_info = pipe;
   gas_arc->type = PIPE;

   /* initialize values */
   pipe->length      = 0.0;
   pipe->diameter    = 0.0;
   pipe->roughness   = 0.0;
   pipe->pressureMax = 200.0;
   pipe->slope       = 0.0;
   pipe->PIdiffVar   = NULL;

   for (arc_elem = SCIPxmlFirstChild(xml_arc); arc_elem != NULL; arc_elem = SCIPxmlNextSibl(arc_elem)) /*lint !e2840*/
   {
      const char* p_elemvalue;
      const char* p_elemunit;
      SCIP_Real heightdiff;

      /*finds element length*/
      if ( strcmp(SCIPxmlGetName(arc_elem), "length") == 0 )
      {
         ++length;
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"length\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            pipe->length = atof(p_elemvalue);
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage("Unit of \"length\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_meter(p_elemunit, &pipe->length) );
            /* calculate the slope */
            heightdiff = gas_arc->targetnode->height - gas_arc->sourcenode->height;
            pipe->slope = calculateSlope( heightdiff, pipe->length);
         }
      }
      /* reads element diameter */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "diameter") == 0 )
      {
         ++diameter;
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"diameter\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            pipe->diameter = atof(p_elemvalue);
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" \"Unit\" of \"diameter\" not found.\n");
               return SCIP_READERROR;
            }

            SCIP_CALL( change_unit_meter(p_elemunit, &pipe->diameter) );
         }
      }
      /* reads element roughness */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "roughness") == 0 )
      {
         ++roughness;
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"roughness\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            pipe->roughness = atof(p_elemvalue);
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
               SCIPerrorMessage(" \"Unit\" of \"roughness\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_meter(p_elemunit, &pipe->roughness) );
         }
      }
      /* reads element pressureMax */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "pressureMax") == 0 )
      {
         ++pressureMax;
         p_elemvalue = SCIPxmlGetAttrval(arc_elem, "value");
         if ( p_elemvalue == NULL )
         {
            SCIPerrorMessage("Element \"pressureMax\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            pipe->pressureMax = atof(p_elemvalue);
            p_elemunit = SCIPxmlGetAttrval(arc_elem, "unit");
            if ( p_elemunit == NULL )
            {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
               SCIPerrorMessage("\"Unit\" of \"pressureMax\" not found.\n");
               return SCIP_READERROR;
#else
               SCIPwarningMessage(scip, "\"Unit\" of \"pressureMax\" not found. Assuming 'barg'.\n");
               SCIP_CALL( change_unit_bar("barg", &pipe->pressureMax) );
#endif
            }
            else
            {
               SCIP_CALL( change_unit_bar(p_elemunit, &pipe->pressureMax) );
            }
         }
      }
      /* reads element flowMin using function readFlowMin */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMin") == 0 )
      {
         ++flowMin;
         SCIP_CALL( readFlowMin(arc_elem, gas_arc, network->normDensity) );
      }
      /* reads element flowMax using function readFlowMax */
      else if ( strcmp(SCIPxmlGetName(arc_elem), "flowMax") == 0 )
      {
         ++flowMax;
         SCIP_CALL( readFlowMax(arc_elem, gas_arc, network->normDensity) );
      }
   }

   if ( length != 1 || diameter != 1 || roughness != 1 || flowMin != 1 || flowMax != 1 || pressureMax > 1 )
   {
      SCIPerrorMessage("Invalid data of pipe %s.\n", gas_arc->id);
      return SCIP_READERROR;
   }

   /* ajust pressure bounds of source- and targetnode if necessary */
   if ( pressureMax == 1 )
   {
      assert( SCIPisPositive(scip, pipe->pressureMax) );

      if ( SCIPisGT(scip, gas_arc->sourcenode->arcPressureMax, pipe->pressureMax) )
         gas_arc->sourcenode->arcPressureMax = pipe->pressureMax;

      if ( SCIPisGT(scip, gas_arc->targetnode->arcPressureMax, pipe->pressureMax) )
         gas_arc->targetnode->arcPressureMax = pipe->pressureMax;
   }

   return SCIP_OKAY;
}


/** read pipe specific data */
static
SCIP_RETCODE readShortPipeData(
   const XML_NODE*       xml_arc,            /**< xml node to read */
   GAS_Arc*              gas_arc,            /**< arc to complete */
   GAS_Network*          network             /**< network */
   )
{
   const XML_NODE* arc_elem;
   int             flowMin = 0;
   int             flowMax = 0;

   gas_arc->detailed_info = NULL;
   gas_arc->type = SHORTPIPE;

   for (arc_elem = SCIPxmlFirstChild(xml_arc); arc_elem != NULL; arc_elem= SCIPxmlNextSibl(arc_elem)) /*lint !e2840*/
   {
      /* reads element flowMin using function readFlowMin*/
      if ( strcmp(SCIPxmlGetName(arc_elem), "flowMin") == 0 )
      {
         ++flowMin;
         SCIP_CALL( readFlowMin(arc_elem, gas_arc, network->normDensity) );
      }

      /* reads element flowMax using function readFlowMax */
      if ( strcmp(SCIPxmlGetName(arc_elem), "flowMax") == 0 )
      {
         ++flowMax;
         SCIP_CALL( readFlowMax(arc_elem, gas_arc, network->normDensity) );
      }
   }

   if ( flowMin != 1 || flowMax != 1 )
   {
      SCIPerrorMessage("Invalid net-data of shortpipe %s.\n", gas_arc->id);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** read arc data (gaslib format) */
static
SCIP_RETCODE readArcData(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network,            /**< network of the node */
   const XML_NODE*       arc,                /**< XML node to read */
   int*                  arcposition,        /**< position in the arc array */
   SCIP_Bool             boxConstraintModel, /**< wheter or not the box constraint model is used */
   SCIP_Bool             algebraic           /**< wheter or not the algebraic model is used */
   )
{
   const char*      arc_id;
   const char*      sourceid;
   const char*      targetid;
   GAS_Node*        source_of_arc;
   GAS_Node*        target_of_arc;
   GAS_Arc*         gas_arc;
   Arc_type         arc_type;
   SCIP_Bool        create_bypass = FALSE;

   /* determine current arc */
   gas_arc = &(network->arcs_ptr[(*arcposition)]);
   gas_arc->next_inarc  = NULL;
   gas_arc->next_outarc = NULL;
   gas_arc->arcposition= (*arcposition);
   gas_arc->detailed_info = NULL;
   gas_arc->flowMax = 0.0;
   gas_arc->flowMin = 0.0;
   gas_arc->flowvar = NULL;
   gas_arc->positiveFlowBinvar = NULL;
   gas_arc->negativeFlowBinvar = NULL;
   gas_arc->mixingRatio = NULL;

   /* find variable name */
   arc_id = SCIPxmlGetAttrval(arc, "id");
   if ( arc_id == NULL )
   {
      SCIPerrorMessage("Attribute \"id\" of %d. arc not found.\n", network->numarcs);
      return SCIP_READERROR;
   }

   /* copy arc_id */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(gas_arc->id), arc_id, strlen(arc_id)+1) );

   /* find source node */
   sourceid = SCIPxmlGetAttrval(arc, "from");
   if ( sourceid == NULL )
   {
      SCIPerrorMessage("Attribute \"from\" of arc %s not found.\n", arc_id);
      return SCIP_READERROR;
   }

   /* find target node */
   targetid = SCIPxmlGetAttrval(arc, "to");
   if ( targetid == NULL )
   {
      SCIPerrorMessage("Attribute \"to\" of arc %s not found.\n", arc_id);
      return SCIP_READERROR;
   }

   /* uses hashtable to associate nodes with incidents arcs */
   if ( ! SCIPhashtableRetrieve(network->Nodeshashtable, (void*) sourceid) )
   {
      SCIPerrorMessage("Source <%s> does not exist!\n", sourceid);
      return SCIP_READERROR;
   }
   else
   {
      source_of_arc = (GAS_Node*) SCIPhashtableRetrieve(network->Nodeshashtable, (void*) sourceid);
      gas_arc->sourcenode = source_of_arc;
      insert_outarc(source_of_arc, gas_arc);
   }

   if ( ! SCIPhashtableRetrieve(network->Nodeshashtable, (void*) targetid) )
   {
      SCIPerrorMessage("Target <%s> does not exist!\n", targetid);
      return SCIP_READERROR;
   }
   else
   {
      target_of_arc = (GAS_Node*) SCIPhashtableRetrieve(network->Nodeshashtable, (void*) targetid);
      gas_arc->targetnode = target_of_arc;
      insert_inarc(target_of_arc, gas_arc);
   }

   /* check that we do not have loops */
   if ( strcmp(sourceid, targetid) == 0 )
   {
      SCIPerrorMessage("Arc <%s> is a loop with node <%s>.\n", arc_id, sourceid);
      return SCIP_READERROR;
   }

   SCIP_CALL( SCIPhashtableInsert(network->Arcshashtable, (void*) gas_arc) );

   /* determines the arc type and calls functions to read type specific information */
   arc_type = convertToArcType(SCIPxmlGetName(arc));
   if ( arc_type == UNKNOWNARCTYPE )
   {
      SCIPerrorMessage("Unknown network element: %s.\n", SCIPxmlGetName(arc));
      return SCIP_READERROR;
   }
   else if ( arc_type == PIPE )
   {
      SCIP_CALL( readPipeData(scip, arc, gas_arc, network) );
      ++(network->numpipes);
      if ( algebraic && ! gas_arc->sourcenode->needPIvar )
      {
         gas_arc->sourcenode->needPIvar = TRUE;
         ++(network->numpivars);
      }
      if ( algebraic && ! gas_arc->targetnode->needPIvar )
      {
         gas_arc->targetnode->needPIvar = TRUE;
         ++(network->numpivars);
      }
   }
   else if ( arc_type == SHORTPIPE )
   {
      SCIP_CALL( readShortPipeData(arc, gas_arc, network) );
      ++(network->numshortpipes);
   }
   else if ( arc_type == VALVE )
   {
      SCIP_CALL( readValveData(scip, arc, gas_arc, network) );
      ++(network->numvalves);
   }
   else if ( arc_type == CONTROLVALVE )
   {
      SCIP_CALL( readControlvalveData(scip, arc, gas_arc, network, &create_bypass) );
      ++(network->numcontrolvalves);
      if ( create_bypass )
      {
         SCIP_CALL( readBypassData(scip, network, gas_arc, (*arcposition)+1) );
         ++(network->numvalves);
         ++(*arcposition);
      }
   }
   else if ( arc_type == CS )
   {
      SCIP_CALL( readCSData(scip, arc, gas_arc, network, &create_bypass, boxConstraintModel) );
      if ( create_bypass )
      {
         SCIP_CALL( readBypassData(scip, network, gas_arc, (*arcposition)+1) );
         ++(network->numvalves);
         ++(*arcposition);
      }
      ++(network->numcompressor);
   }
   else
   {
      assert( arc_type == RESISTOR );
      SCIP_CALL( readResistorData(scip, arc, gas_arc, network) );
   }

   ++(target_of_arc->numincidentarcs);
   ++(source_of_arc->numincidentarcs);

   return SCIP_OKAY;
}

/** determine number of neighbors and maximal number */
static
SCIP_RETCODE determineNeighbors(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network             /**< network */
   )
{
   int* markneighbors;
   int i;

   assert( scip != NULL );
   assert( network != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &markneighbors, network->numnodes) );

   for (i = 0; i < network->numnodes; ++i)
      markneighbors[i] = -1;

   /* loop through nodes */
   for (i = 0; i < network->numnodes; ++i)
   {
      GAS_Arc* arc;
      GAS_Node* node;
      GAS_Node* neighbor;
      int numneighbors = 0;

      node = &network->nodes_ptr[i];

      arc = node->inarcs;
      while ( arc != NULL )
      {
         neighbor = arc->sourcenode;

         if ( markneighbors[neighbor->nodeposition] != i )
         {
            markneighbors[neighbor->nodeposition] = i;
            ++numneighbors;
         }

         arc = arc->next_inarc;
      }

      arc = node->outarcs;
      while ( arc != NULL )
      {
         neighbor = arc->targetnode;

         if ( markneighbors[neighbor->nodeposition] != i )
         {
            markneighbors[neighbor->nodeposition] = i;
            ++numneighbors;
         }

         arc = arc->next_outarc;
      }

      node->numneighbors = numneighbors;
   }

   SCIPfreeBufferArray(scip, &markneighbors);

   return SCIP_OKAY;
}

/** read gas data (gaslib format) */
static
SCIP_RETCODE GASLIBreadFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of file to read */
   SCIP_PROBDATA*        probdata,           /**< problem data to be filled */
   SCIP_Bool             boxConstraintModel, /**< wheter or not the box constraint model is used */
   SCIP_Bool             algebraic           /**< wheter or not the algebraic model is used */
   )
{
   XML_NODE*       start;                    /* points to first xml_node in Gas file */
   const char*     tag;
   const XML_NODE* frame_arcs;
   const XML_NODE* arc;
   const XML_NODE* frame_nodes;
   const XML_NODE* source;
   GAS_Network*    network;
   int             nodeposition      = 0;
   int             arcposition       = 0;
   int             foundGasConstants = 0;
   int             i;
   Arc_type        arc_type;
   const char*     elemvalue;

   assert( scip != NULL );
   assert( filename != NULL );
   assert( probdata != NULL );

   /* get memory for GAS_Network and set physical values as default */
   SCIP_CALL( SCIPallocBlockMemory(scip, &network) );
   network->numnodes                  = 0;
   network->numarcs                   = 0;
   network->maxnumincidentarcs        = 0;
   network->numflowdirvars            = 0;
   network->numpipes                  = 0;
   network->numshortpipes             = 0;
   network->numvalves                 = 0;
   network->numcontrolvalves          = 0;
   network->numcompressor             = 0;
   network->numcsconfigurations       = 0;
   network->numcsvarsfornlrs          = 0;
   network->numnonlinresistor         = 0;
   network->numlinresistor            = 0;
   network->numpivars                 = 0;
   network->normDensity               = 0.783;
   network->pseudocriticalPressure    = 46.4512;
   network->pseudocriticalTemperature = 192.033;
   network->gasTemperature            = 283.15;
   network->molarMass1                = -1.0;
   network->molarMass2                = -1.0;
   network->specificHeatCap1          = -1.0;
   network->specificHeatCap2          = -1.0;
   network->nodes_ptr                 = NULL;
   network->arcs_ptr                  = NULL;
   network->spanningTree              = NULL;
   network->numBasisCycles            = 0;
   network->firstCycle                = NULL;
   network->lastCycle                 = NULL;
   network->firstCombinedCycle        = NULL;

   /* read xml file */
   start = SCIPxmlProcess(filename);

   if ( start == NULL )
   {
      SCIPerrorMessage("Some error occured during parsing <%s>.\n", filename);
      return SCIP_READERROR;
   }

#if SCIP_VERSION >= 400
   /* creates hastable for finding node ids */
   SCIP_CALL( SCIPhashtableCreate(&(network->Nodeshashtable), SCIPblkmem(scip), 1000, SCIPhashGetKeyNode, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );
    /* creates hastable for finding arc ids */
   SCIP_CALL( SCIPhashtableCreate(&(network->Arcshashtable), SCIPblkmem(scip), 1000, SCIPhashGetKeyArc, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );
#else
   /* creates hastable for finding node ids */
   SCIP_CALL( SCIPhashtableCreate(&(network->Nodeshashtable), SCIPblkmem(scip), SCIPcalcHashtableSize(1000), SCIPhashGetKeyNode, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );
   /* creates hashtable for finding arc ids */
   SCIP_CALL( SCIPhashtableCreate(&(network->Arcshashtable), SCIPblkmem(scip), SCIPcalcHashtableSize(1000), SCIPhashGetKeyArc, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );
#endif

   /* finds sources */
   tag = "framework:nodes";
   frame_nodes = SCIPxmlFindNodeMaxdepth(start, tag, 0, 3);
   if ( frame_nodes == NULL )
   {
      /* free xml data */
      SCIPxmlFreeNode(start);

      SCIPerrorMessage("Nodes section not found.\n");
      return SCIP_READERROR;
   }

   /* counts number of nodes */
   for (source = SCIPxmlFirstChild(frame_nodes); source != NULL; source = SCIPxmlNextSibl(source)) /*lint !e2840*/
      ++(network->numnodes);

   /* allocates memory for nodes */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &network->nodes_ptr, network->numnodes) );

   for (source = SCIPxmlFirstChild(frame_nodes); source != NULL; source = SCIPxmlNextSibl(source)) /*lint !e2840*/
   {
      SCIP_CALL( readNodeData(scip, network, source, nodeposition, &foundGasConstants) );
      ++nodeposition;
   }
   assert( nodeposition <= network->numnodes );

   /* check if norm density and so on were found */
   if ( foundGasConstants == 0 )
   {
      SCIPerrorMessage("No gas constants, e.g., norm density, were found!\n");
      return SCIP_READERROR;
   }
   else if ( foundGasConstants > 1 )
   {
      SCIPdebugMessage("Multiple definitions of the gas constants in the net-file.\n");
   }

   /* finds connections */
   tag = "framework:connections";
   frame_arcs = SCIPxmlFindNodeMaxdepth(start, tag, 0, 3);
   if ( frame_arcs == NULL )
   {
      /* free xml data */
      SCIPxmlFreeNode(start);

      SCIPerrorMessage("Arc section not found.\n");
      return SCIP_READERROR;
   }

   /* count number of all arcs */
   for (arc = SCIPxmlFirstChild(frame_arcs); arc != NULL; arc = SCIPxmlNextSibl(arc)) /*lint !e2840*/
   {
      arc_type = convertToArcType(SCIPxmlGetName(arc));
      if ( arc_type == CS || arc_type == CONTROLVALVE )
      {
         elemvalue = SCIPxmlGetAttrval(arc, "internalBypassRequired");
         if ( elemvalue == NULL  || atoi(elemvalue) == 1 )
            ++(network->numarcs);
      }
      ++(network->numarcs);
   }

   /* allocate memory for arcs */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &network->arcs_ptr, network->numarcs) );

   for (arc = SCIPxmlFirstChild(frame_arcs); arc != NULL; arc = SCIPxmlNextSibl(arc)) /*lint !e2840*/
   {
      SCIP_CALL( readArcData(scip, network, arc, &arcposition, boxConstraintModel, algebraic) );
      ++arcposition;
   }

   /* now set the number of positive and negative flow direction binvars which is needed */
   if ( probdata->noFlowBinvars )
   {
      network->numflowdirvars = 0;
   }
   else
   {  /* changed the number in order to generate flowbinvars for every arc */
      network->numflowdirvars = network->numarcs; // - (network->numcontrolvalves + network->numcompressor);
   }

   /* set the maximal number of incident arcs to a node */
   for (i = 0; i < network->numnodes; ++i)
   {
      if ( network->nodes_ptr[i].numincidentarcs > network->maxnumincidentarcs )
         network->maxnumincidentarcs = network->nodes_ptr[i].numincidentarcs;
   }

   /* determine neighbors of all nodes (used for preprocessing later) */
   SCIP_CALL( determineNeighbors(scip, network) );

   /* free xml data */
   SCIPxmlFreeNode(start);

   probdata->network = network;

   return SCIP_OKAY;
}


/** read scenario (gaslib format) */
static
SCIP_RETCODE SCNreadFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of file to read */
   SCIP_PROBDATA*        probdata            /**< problem data to be filled */
   )
{
   XML_NODE*       start;                    /* points to first xml_node in Gas file */
   const char*     tag;
   const XML_NODE* scnnode;
   const XML_NODE* scen;
   const XML_NODE* nodeprop;
   SCIP_Real       flowChecksum = 0.0;
   SCIP_Real       convalue;

   assert( scip != NULL );
   assert( filename != NULL );
   assert( probdata != NULL );
   assert( probdata->network != NULL );

   /* read xml file */
   start = SCIPxmlProcess(filename);

   if ( start == NULL )
   {
      SCIPerrorMessage("Some error occured during parsing <%s>.\n", filename);
      return SCIP_READERROR;
   }

   /* finds nodes */
   tag = "scenario";
   scen = SCIPxmlFindNodeMaxdepth(start, tag, 0, 3);
   if ( scen == NULL )
   {
      /* free xml data */
      SCIPxmlFreeNode(start);
      SCIPerrorMessage("Nodes section not found.\n");
      return SCIP_READERROR;
   }

   /* loop through all nodes */
   for (scnnode = SCIPxmlFirstChild(scen); scnnode != NULL; scnnode = SCIPxmlNextSibl(scnnode)) /*lint !e2840*/
   {
      const char* attrvalue;
      Node_type   attrtype;
      const char* propname;
      const char* propvalue;
      const char* unit;
      GAS_Node*   N;
      SCIP_Real   lowerflowbound = -1.0;
      SCIP_Real   upperflowbound = -1.0;
      int molarMassCounter = 0;

      /* determines current node */
      if ( strcmp(SCIPxmlGetName(scnnode), "node") == 0 )
      {
         if ( strcmp(SCIPxmlGetAttrval(scnnode, "type"), "gasProperties") == 0 )
         {
            /* loop through node attributes */
            for (nodeprop = SCIPxmlFirstChild(scnnode); nodeprop != NULL && molarMassCounter < 3; nodeprop = SCIPxmlNextSibl(nodeprop))  /*lint !e2840*/
            {
               if ( strcmp(SCIPxmlGetName(nodeprop), "molarMass1") == 0 )
               {
                  propvalue = SCIPxmlGetAttrval(nodeprop, "value");
                  unit = SCIPxmlGetAttrval(nodeprop, "unit");

                  /* change unit, if necessary */
                  convalue = atof(propvalue);
                  SCIP_CALL( change_molarmass_unit(unit, &convalue) );
                  probdata->molarMass1 = convalue;
                  probdata->network->molarMass1 = convalue;

                  ++molarMassCounter;
               }
               else if ( strcmp(SCIPxmlGetName(nodeprop), "molarMass2") == 0 )
               {
                  propvalue = SCIPxmlGetAttrval(nodeprop, "value");
                  unit = SCIPxmlGetAttrval(nodeprop, "unit");

                  /* change unit, if necessary */
                  convalue = atof(propvalue);
                  SCIP_CALL( change_molarmass_unit(unit, &convalue) );
                  probdata->molarMass2 = convalue;
                  probdata->network->molarMass2 = convalue;

                  ++molarMassCounter;
               }
            }
            assert( molarMassCounter == 2 );
            assert( probdata->molarMass1 >= probdata->molarMass2 );
         }
         else
         {
            attrvalue = SCIPxmlGetAttrval( scnnode, "id" );
            if ( attrvalue == NULL )
            {
               SCIPerrorMessage( "Attribute \"id\" of some node in the scn-file not found.\n" );
               return SCIP_READERROR;
            }
            else if ( ! SCIPhashtableRetrieve(probdata->network->Nodeshashtable, (void*)attrvalue) )
            {
               SCIPerrorMessage( "Node <%s> in the scn-file does not exist in the net-file!\n", attrvalue );
               return SCIP_READERROR;
            }

            N = (GAS_Node*)SCIPhashtableRetrieve(probdata->network->Nodeshashtable, (void*)attrvalue);
            assert( N != NULL );

            /* determines node type */
            attrvalue = SCIPxmlGetAttrval(scnnode, "type");
            if ( attrvalue == NULL )
            {
               SCIPerrorMessage("Attribute \"type\" not found.\n");
               return SCIP_READERROR;
            }

            attrtype = convertToNodeType(attrvalue);
            if ( attrtype == UNKNOWNNODETYPE )
            {
               SCIPerrorMessage( "Unkown node type.\n" );
               return SCIP_READERROR;
            }

            if ( attrtype != ENTRY && attrtype != EXIT )
            {
               SCIPerrorMessage("Node type %d not allowed in scn-file.\n", attrtype);
               return SCIP_READERROR;
            }

            N->type = attrtype;   /* store type */

            /* loop through node attributes */
            for (nodeprop = SCIPxmlFirstChild(scnnode); nodeprop != NULL; nodeprop = SCIPxmlNextSibl(nodeprop)) /*lint !e2840*/
            {
               /* looks for xml node "pressure" */
               if ( strcmp(SCIPxmlGetName(nodeprop), "pressure") == 0 )
               {
                  /* find bound type */
                  propname = SCIPxmlGetAttrval(nodeprop, "bound");
                  if ( propname == NULL )
                  {
                     SCIPerrorMessage("No bound on the pressure of node %s found.\n", N->id);
                     return SCIP_READERROR;
                  }

                  /* find bound */
                  propvalue = SCIPxmlGetAttrval(nodeprop, "value");
                  if ( propvalue == NULL )
                  {
                     SCIPerrorMessage("Lower or upper pressure bound of node %s not found.\n", N->id);
                     return SCIP_READERROR;
                  }

                  /* find pressure unit */
                  unit = SCIPxmlGetAttrval(nodeprop, "unit");
                  if ( unit == NULL )
                  {
                     SCIPerrorMessage("Unit of %s pressure bound of node %s not found.\n", propname, N->id);
                     return SCIP_READERROR;
                  }

                  if ( strcmp(propname, "lower") == 0 )
                  {
                     N->scnPressureMin = atof(propvalue);
                     SCIP_CALL( change_unit_bar(unit, &N->scnPressureMin) );

                     if ( SCIPisLT(scip, N->scnPressureMin, N->netPressureMin) )
                     {
                        SCIPdebugMessage( "Lower pressure from scn file (%f) is less than minimal pressure (%f).\n", N->scnPressureMin, N->netPressureMin );
                     }
                  }
                  else if ( strcmp(propname, "upper") == 0 )
                  {
                     N->scnPressureMax = atof(propvalue);
                     SCIP_CALL( change_unit_bar(unit, &N->scnPressureMax) );

                     if ( SCIPisGT(scip, N->scnPressureMax, N->netPressureMax) )
                     {
                        SCIPdebugMessage("Upper pressure on node %s from scn file (%f) is larger than minimal pressure (%f).\n",
                           N->id, N->scnPressureMax, N->netPressureMax);
                     }
                  }
                  else
                  {
                     SCIPerrorMessage("\"%s\" is not a valid bound type on the pressure of node %s.\n", propname, N->id);
                     return SCIP_READERROR;
                  }
               }
               /* looks for xml node "flow" */
               else if ( strcmp(SCIPxmlGetName(nodeprop), "flow") == 0 )
               {
                  /* determine unit of the flow */
                  unit = SCIPxmlGetAttrval(nodeprop, "unit");
                  if ( unit == NULL )
                  {
                     SCIPerrorMessage("No unit of the flow at node %s found.\n", N->id);
                     return SCIP_READERROR;
                  }

                  /* determine flow bound type, i.e., lower, upper or bound */
                  propname = SCIPxmlGetAttrval(nodeprop, "bound");
                  if ( propname == NULL )
                  {
                     SCIPerrorMessage("No bound on the flow of node %s found.\n", N->id);
                     return SCIP_READERROR;
                  }

                  /* check for flow value */
                  propvalue = SCIPxmlGetAttrval(nodeprop, "value");
                  if ( propvalue == NULL )
                  {
                     SCIPerrorMessage("No bound on the flow of node %s found.\n", N->id);
                     return SCIP_READERROR;
                  }

                  /* check the bound type */
                  if ( strcmp( propname, "both" ) == 0 )
                  {
                     N->flow = atof(propvalue);
                     SCIP_CALL( change_flow_unit(unit, &N->flow, probdata->network->normDensity) );
                     if ( probdata->maxFlow )
                     {
                        lowerflowbound = atof(propvalue);
                        upperflowbound = atof(propvalue);
                        SCIP_CALL( change_flow_unit(unit, &upperflowbound, probdata->network->normDensity) );
                        SCIP_CALL( change_flow_unit(unit, &lowerflowbound, probdata->network->normDensity) );
                     }
                  }
                  else if ( ! strcmp(propname, "lower") )
                  {
                     lowerflowbound = atof(propvalue);
                     SCIP_CALL( change_flow_unit(unit, &lowerflowbound, probdata->network->normDensity) );
                  }
                  else if ( ! strcmp(propname, "upper") )
                  {
                     upperflowbound = atof(propvalue);
                     SCIP_CALL( change_flow_unit(unit, &upperflowbound, probdata->network->normDensity) );
                  }
                  else
                  {
                     SCIPerrorMessage("\"%s\" is not a valid bound for the in/outflow at node %s.\n", propname, N->id);
                     return SCIP_READERROR;
                  }
               }
               /* overrides the molar mass set in readNodeData */
               /* molar mass */
               else if ( strcmp(SCIPxmlGetName(nodeprop), "molarMass") == 0 )
               {
                  propvalue = SCIPxmlGetAttrval(nodeprop, "value");
                  unit = SCIPxmlGetAttrval(nodeprop, "unit");
                  if ( unit == NULL )
                  {
#ifdef ERROR_ON_UNSPECIFIED_UNIT
                     SCIPerrorMessage("Unit of molarMass not found.\n");
                     return SCIP_READERROR;
#else
                     SCIPwarningMessage(scip, "Unit of molarMass not found. Assuming 'kg_per_kmol'.\n");
                     convalue = atof( elemvalue );
#endif
                  }
                  else
                  {
                     /* change unit, if necessary */
                     convalue = atof(propvalue);
                     SCIP_CALL( change_molarmass_unit(unit, &convalue) );

                     assert(convalue <= probdata->molarMass1);
                     assert(convalue >= probdata->molarMass2);
                     N->molarMass = convalue;
                  }
               }
            }

            /* set the correct flow bound */
            if ( probdata->maxFlow )
            {
               if ( ! SCIPisLE(scip, lowerflowbound, upperflowbound) )
               {
                  SCIPerrorMessage("Lower larger than  upper bound on flow of node %s.\n", N->id);
                  return SCIP_READERROR;
               }
               else if ( lowerflowbound >= 0.0 )
               {
                  N->flow = upperflowbound * probdata->scalescn;
                  N->scnFlowMax = upperflowbound * probdata->scalescn;
                  N->scnFlowMin = lowerflowbound * probdata->scalescn;
               }
            }
            else
            {
               if ( ! SCIPisEQ(scip, lowerflowbound, upperflowbound) )
               {
                  SCIPerrorMessage("Different lower and upper bound on flow of node %s.\n", N->id);
                  return SCIP_READERROR;
               }
               else if ( lowerflowbound >= 0.0 )
               {
                  N->flow = lowerflowbound * probdata->scalescn;
                  N->scnFlowMax = upperflowbound * probdata->scalescn;
                  N->scnFlowMin = lowerflowbound * probdata->scalescn;
               }
            }

            /* if the node is a sink, we need to multiply the flow value by -1 */
            if ( N->type == EXIT )
            {
               N->flow *= -1.0;
               N->scnFlowMax *= -1.0;
               N->scnFlowMin *= -1.0;
            }

            flowChecksum += N->flow;

            SCIPdebugMessage("Flow at boundary node %s is %f.\n", N->id, N->flow);
         }
      }
   }

   if ( ! probdata->maxFlow )
   {
      /* the inflow has to be equal to the outflow */
      if ( ! SCIPisFeasZero(scip, flowChecksum) )
      {
         SCIPerrorMessage("The in- and outflow doesn't add up to zero (%10.9f).\n", flowChecksum);
         return SCIP_READERROR;
      }
   }

   /* end xml process */
   SCIPxmlFreeNode(start);

   return SCIP_OKAY;
}


/** read the box constraints for a single configuration */
static
SCIP_RETCODE CSreadConfigBoxData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data to be filled*/
   const XML_NODE*       prop,               /**< XML node to read */
   GAS_Arc*              arc,                /**< gas_arc corresponding to the compressor */
   GAS_CSBoxCons*        boxcons             /**< configuration struct */
   )
{
   const XML_NODE* bound;                    /* element in xml file */
   const char*     attrvalue;                /* attribute value of xml node */
   const char*     boundunit;                /* unit of xml node */
   int             i;

   for (bound = SCIPxmlFirstChild(prop); bound != NULL; bound = SCIPxmlNextSibl(bound)) /*lint !e2840*/
   {
      /* reads massFlowMin */
      if ( strcmp(SCIPxmlGetName(bound), "massFlowMin") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"massFlowMin\" not found.\n");
            return SCIP_READERROR;
         }
         boxcons->massFlowMin = atof(attrvalue);
         /* determine unit */
         boundunit = SCIPxmlGetAttrval(bound, "unit");
         if ( boundunit == NULL )
         {
            SCIPerrorMessage("Unit of \"massFlowMin\" not found.\n");
            return SCIP_READERROR;
         }
         SCIP_CALL( change_flow_unit(boundunit, &boxcons->massFlowMin, probdata->network->normDensity) );
         assert( SCIPisGE(scip, boxcons->massFlowMin, 0.0) );
      }
      /* reads massFlowMax */
      else if ( strcmp(SCIPxmlGetName(bound), "massFlowMax") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"massFlowMax\" not found.\n");
            return SCIP_READERROR;
         }
         boxcons->massFlowMax = atof(attrvalue);

         /* determine unit */
         boundunit = SCIPxmlGetAttrval(bound, "unit");
         if ( boundunit == NULL )
         {
            SCIPerrorMessage("Unit of \"massFlowMax\" not found.\n");
            return SCIP_READERROR;
         }
         SCIP_CALL( change_flow_unit(boundunit, &boxcons->massFlowMax, probdata->network->normDensity) );
      }
      /* reads minimal incoming pressure */
      else if ( strcmp(SCIPxmlGetName(bound), "pressureInMin") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"pressureInMin\" not found.\n");
            return SCIP_READERROR;
         }
         else
         {
            boxcons->pressureInMin = atof(attrvalue);

            /* determines unit */
            boundunit = SCIPxmlGetAttrval(bound, "unit");
            if ( boundunit == NULL )
            {
               SCIPerrorMessage("Unit of \"pressureInMin\" not found.\n");
               return SCIP_READERROR;
            }
            SCIP_CALL( change_unit_bar(boundunit, &boxcons->pressureInMin) );
         }
      }
      /* reads maximal incoming pressure */
      else if ( strcmp(SCIPxmlGetName(bound), "pressureInMax") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"pressureInMax\" not found.\n");
            return SCIP_READERROR;
         }
         boxcons->pressureInMax = atof(attrvalue);

         /* determines unit */
         boundunit = SCIPxmlGetAttrval(bound, "unit");
         if ( boundunit == NULL )
         {
            SCIPerrorMessage("Unit of \"pressureInMax\" not found.\n");
            return SCIP_READERROR;
         }
         SCIP_CALL( change_unit_bar(boundunit, &boxcons->pressureInMax) );
         SCIPdebugMessage("Maximal incoming pressure of <%s> is <%f>\n", arc->id, boxcons->pressureInMax);
      }
      /* reads minimal outgoing pressure */
      else if ( strcmp(SCIPxmlGetName(bound), "pressureOutMin") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"pressureOutMin\" not found.\n");
            return SCIP_READERROR;
         }
         boxcons->pressureOutMin = atof(attrvalue);
         boundunit = SCIPxmlGetAttrval(bound, "unit");
         if ( boundunit == NULL )
         {
            SCIPerrorMessage("Unit of \"pressureOutMin\" not found.\n");
            return SCIP_READERROR;
         }
         SCIP_CALL( change_unit_bar(boundunit, &boxcons->pressureOutMin) );
         SCIPdebugMessage("Minimal outcoming pressure of <%s> is <%f>\n", arc->id, boxcons->pressureOutMin);
      }
      /* reads maximal outgoing pressure */
      else if ( strcmp(SCIPxmlGetName(bound), "pressureOutMax") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"pressureOutMax\" not found.\n");
            return SCIP_READERROR;
         }
         boxcons->pressureOutMax = atof(attrvalue);
         boundunit = SCIPxmlGetAttrval(bound, "unit");
         if ( boundunit == NULL )
         {
            SCIPerrorMessage("Unit of \"pressureOutMax\" not found.\n");
            return SCIP_READERROR;
         }
         SCIP_CALL( change_unit_bar(boundunit, &boxcons->pressureOutMax) );
         SCIPdebugMessage("Maximal outcoming pressure of <%s> is <%f>\n", arc->id, boxcons->pressureOutMax);
      }
      /* reads minimal absolute pressure increase */
      else if ( strcmp(SCIPxmlGetName(bound), "pressureIncAbsMin") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"pressureIncAbsMin\" not found.\n");
            return SCIP_READERROR;
         }
         boxcons->pressureIncAbsMin = atof(attrvalue);
         boundunit = SCIPxmlGetAttrval(bound, "unit");
         if ( boundunit == NULL )
         {
            SCIPerrorMessage("Unit of \"pressureIncAbsMin\" not found.\n");
            return SCIP_READERROR;
         }
         SCIP_CALL( change_unit_bar(boundunit, &boxcons->pressureIncAbsMin) );
         /* assert ( boxcons->pressureIncAbsMin > 0.0); */
         SCIPdebugMessage("Minimal bound on absolute pressure increase of <%s> is <%f>\n", arc->id, boxcons->pressureIncAbsMin);
      }
      /* reads maximal absolute pressure increase */
      else if ( strcmp(SCIPxmlGetName(bound), "pressureIncAbsMax") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"pressureIncAbsMax\" not found.\n");
            return SCIP_READERROR;
         }
         boxcons->pressureIncAbsMax = atof(attrvalue);
         boundunit = SCIPxmlGetAttrval(bound, "unit");
         if ( boundunit == NULL )
         {
            SCIPerrorMessage("Unit of \"pressureIncAbsMax\" not found.\n");
            return SCIP_READERROR;
         }
         SCIP_CALL( change_unit_bar(boundunit, &boxcons->pressureIncAbsMax) );
         assert ( boxcons->pressureIncAbsMax > 0.0);
         SCIPdebugMessage("Maximal bound on absolute pressure increase of <%s> is <%f>\n", arc->id, boxcons->pressureIncAbsMax);
      }
      /* reads minimal relative pressure increase */
      else if ( strcmp(SCIPxmlGetName(bound), "pressureIncRelMin") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"pressureIncRelMin\" not found.\n");
            return SCIP_READERROR;
         }
         boxcons->pressureIncRelMin = atof(attrvalue);
         boundunit = SCIPxmlGetAttrval(bound, "unit");
         if ( boundunit == NULL )
         {
            SCIPerrorMessage("Unit of \"pressureIncRelMin\" not found.\n");
            return SCIP_READERROR;
         }
         SCIP_CALL( change_unit_bar(boundunit, &boxcons->pressureIncRelMin) );
         /* assert ( boxcons->pressureIncRelMin > 1.0); */
         SCIPdebugMessage("Minimal bound on relative pressure increase of <%s> is <%f>\n", arc->id, boxcons->pressureIncRelMin);
      }
      /* reads maximal relative pressure increase */
      else if ( strcmp(SCIPxmlGetName(bound), "pressureIncRelMax") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(bound, "value");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Value of element \"pressureIncRelMax\" not found.\n");
            return SCIP_READERROR;
         }
         boxcons->pressureIncRelMax = atof(attrvalue);
         boundunit = SCIPxmlGetAttrval(bound, "unit");
         if ( boundunit == NULL )
         {
            SCIPerrorMessage("Unit of \"pressureIncRelMax\" not found.\n");
            return SCIP_READERROR;
         }
         SCIP_CALL( change_unit_bar(boundunit, &boxcons->pressureIncRelMax) );
         assert ( boxcons->pressureIncRelMax > 1.0);
         SCIPdebugMessage("Maximal bound on relative pressure increase of <%s> is <%f>\n", arc->id, boxcons->pressureIncRelMax);
      }
      /* reads information about additional facets */
      else if ( strcmp(SCIPxmlGetName(bound), "additionalFacets") == 0 )
      {
         const XML_NODE*       facet;              /* element in xml file */
         const char*           facetvalue;

         /* only read facets for box constraint model in variable ppq */
         if ( strcmp(SCIPxmlGetAttrval(bound, "space"), "ppq") != 0 )
            continue;

         boxcons->numfacets = 0;

         /* counts number of facets to add */
         for (facet = SCIPxmlFirstChild(bound); facet != NULL; facet = SCIPxmlNextSibl(facet)) /*lint !e2840*/
         {
            if ( strcmp(SCIPxmlGetName(facet), "facet") == 0 )
               ++boxcons->numfacets;
         }

         /* gets memory for facet coefficients */
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(boxcons->facetcoeff), boxcons->numfacets) );
         for (i = 0; i < boxcons->numfacets; i++)
         {
            SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(boxcons->facetcoeff[i]), 3));
         }
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(boxcons->rhs), boxcons->numfacets) );

         /* reads facet information */
         i = 0;
         for (facet = SCIPxmlFirstChild(bound); facet != NULL; facet = SCIPxmlNextSibl(facet)) /*lint !e2840*/
         {
            if ( strcmp(SCIPxmlGetName(facet), "variables") == 0 )
            {
               const XML_NODE*       variable;
               const char*           coeff;
               const char*           name;
               const char*           unit;

               for (variable = SCIPxmlFirstChild(facet); variable != NULL; variable = SCIPxmlNextSibl(variable)) /*lint !e2840*/
               {
                  if ( strcmp(SCIPxmlGetName(variable), "variable") == 0 )
                  {
                     coeff = SCIPxmlGetAttrval(variable, "coeff");
                     if ( coeff == NULL )
                     {
                        SCIPerrorMessage("variable coefficient not found.\n");
                        return SCIP_READERROR;
                     }
                     name = SCIPxmlGetAttrval(variable, "name");
                     if ( name == NULL )
                     {
                        SCIPerrorMessage("name not found.\n");
                        return SCIP_READERROR;
                     }
                     unit = SCIPxmlGetAttrval(variable, "unit");
                     if ( unit == NULL )
                     {
                        SCIPerrorMessage("unit not found.\n");
                        return SCIP_READERROR;
                     }
                     if ( strcmp(coeff, "a") == 0 )
                     {
                        SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(boxcons->a), name, strlen(name) + 1) );
                        SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(boxcons->a_unit), unit, strlen(unit) + 1) );
                     }
                     if ( strcmp(coeff, "b") == 0 )
                     {
                        SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(boxcons->b), name, strlen(name) + 1) );
                        SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(boxcons->b_unit), unit, strlen(unit) + 1) );
                     }
                     if ( strcmp(coeff, "c") == 0 )
                     {
                        SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(boxcons->c), name, strlen(name) + 1) );
                        SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(boxcons->c_unit), unit, strlen(unit) + 1) );
                     }
                  }
               }
            }
            else if ( strcmp(SCIPxmlGetName(facet), "facet") == 0 )
            {
               /* looking for coefficient a */
               facetvalue = SCIPxmlGetAttrval(facet, "a");
               if ( facetvalue == NULL )
               {
                  SCIPerrorMessage("coefficient \"a\" not found.\n");
                  return SCIP_READERROR;
               }
               boxcons->facetcoeff[i][0] = atof(facetvalue);

               /* looking for coefficient b */
               facetvalue = SCIPxmlGetAttrval(facet, "b");
               if ( facetvalue == NULL )
               {
                  SCIPerrorMessage("coefficient \"b\" not found.\n");
                  return SCIP_READERROR;
               }
               boxcons->facetcoeff[i][1] = atof(facetvalue);

               /* looking for coefficient c */
               facetvalue = SCIPxmlGetAttrval(facet, "c");
               if ( facetvalue == NULL )
               {
                  SCIPerrorMessage("coefficient \"c\" not found.\n");
                  return SCIP_READERROR;
               }
               boxcons->facetcoeff[i][2] = atof(facetvalue);

               /* checking if relation is "le" */
               facetvalue = SCIPxmlGetAttrval(facet, "rel");
               if ( facetvalue == NULL )
               {
                  SCIPerrorMessage("relation not found.\n");
                  return SCIP_READERROR;
               }
               if ( strcmp(facetvalue, "le") != 0)
               {
                  SCIPerrorMessage("expected 'le' not found.\n");
                  return SCIP_READERROR;
               }

               /* looking for coefficient rhs */
               facetvalue = SCIPxmlGetAttrval(facet, "rhs");
               if ( facetvalue == NULL )
               {
                  SCIPerrorMessage("righthand side not found.\n");
                  return SCIP_READERROR;
               }
               boxcons->rhs[i] = atof(facetvalue);
               i++;
            }
         }
      }
   }
   return SCIP_OKAY;
}

/** read data for boxConstraintModel for all compressor stations
 *  and all configurations from compressor station (CS) file
 */
static
SCIP_RETCODE CSreadBoxModelData(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of file to read */
   SCIP_PROBDATA*        probdata            /**< problem data to be filled */
   )
{
   XML_NODE* start;                          /* points to first xml_node in Gas file */
   const char* tag;
   const char* tag2;
   const XML_NODE* csnode;
   const XML_NODE* child;
   const XML_NODE* nodeprop;
   const XML_NODE* prop;
   const XML_NODE* configurations;
   const XML_NODE* config;
   SCIP_Real gasTemperature = 0.0;
   SCIP_Bool equalTemperature = FALSE;
   SCIP_Real tempdiff;
   SCIP_Real min;
   SCIP_Real closestTemperature = 0;
   int k;

   assert( scip != NULL );
   assert( filename != NULL );
   assert( probdata != NULL );
   assert( probdata->network != NULL );

   /* read xml file */
   start = SCIPxmlProcess(filename);
   if ( start == NULL )
   {
      SCIPerrorMessage("Some error occured during parsing <%s>.\n", filename);
      return SCIP_READERROR;
   }

   /* find nodes */
   tag = "compressorStations";
   child = SCIPxmlFindNodeMaxdepth(start, tag, 0, 4);
   if ( child == NULL )
   {
      /* free xml data */
      SCIPxmlFreeNode(start);

      SCIPerrorMessage("Nodes section not found.\n");
      return SCIP_READERROR;
   }

   /* loop through all sources */
   for ( csnode = SCIPxmlFirstChild(child); csnode != NULL; csnode = SCIPxmlNextSibl(csnode) ) /*lint !e2840*/
   {
      const char*         attrvalue;
      const char*         compressorname;
      const char*         propvalue;
      const char*         unit;
      GAS_CS*             compressor;
      GAS_Arc*            arc;
      GAS_CSBoxCons*      boxcons;

      /* Need to check, if the child of "child" represents a node etc. */
      if ( strcmp(SCIPxmlGetName(csnode), "compressorStation") == 0 )
      {
         attrvalue = SCIPxmlGetAttrval(csnode, "id");
         if ( attrvalue == NULL )
         {
            SCIPerrorMessage("Attribute \"id\" in cs file not found.\n");
            return SCIP_READERROR;
         }
         /* copy id */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(compressorname), attrvalue, strlen(attrvalue) + 1) );
         SCIPdebugMessage ("csname: <%s>\n", compressorname );

         arc = (GAS_Arc*) SCIPhashtableRetrieve(probdata->network->Arcshashtable, (void*) compressorname);
         if ( arc == NULL )
         {
            SCIPerrorMessage("Compressor <%s> does not exist!\n", compressorname);
            return SCIP_READERROR;
         }

         SCIPdebugMessage("arc id <%s>\n", arc->id);
         compressor = arc->detailed_info;

         configurations = SCIPxmlFindNodeMaxdepth(csnode, "configurations", 0, 1);
         for (config = SCIPxmlFirstChild(configurations); config != NULL; config = SCIPxmlNextSibl(config)) /*lint !e2840*/
         {
            ++(compressor->numconfigurations);
            ++(probdata->network->numcsconfigurations);
         }

         k = 0;
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &compressor->boxcons, compressor->numconfigurations) );
         for (config = SCIPxmlFirstChild(configurations); config != NULL; config = SCIPxmlNextSibl(config)) /*lint !e2840*/
         {
            boxcons = &(compressor->boxcons[k]);

            attrvalue = SCIPxmlGetAttrval(config, "confId");
            if ( attrvalue == NULL )
            {
               SCIPerrorMessage("Attribute \"confId\" in cs file not found.\n");
               return SCIP_READERROR;
            }

            boxcons->massFlowMin         = 0.0;
            boxcons->massFlowMax         = 0.0;
            boxcons->pressureInMin       = 1.0;
            boxcons->pressureInMax       = 120.0;
            boxcons->pressureOutMin      = 1.0;
            boxcons->pressureOutMax      = 120.0;
            boxcons->pressureIncAbsMin   = 0.0;
            boxcons->pressureIncAbsMax   = 120.0;
            boxcons->pressureIncRelMin   = 0.0;
            boxcons->pressureIncRelMax   = 0.0;
            boxcons->config_binvar       = NULL;
            boxcons->facetcoeff          = NULL;
            boxcons->rhs                 = NULL;
            boxcons->a                   = NULL;
            boxcons->a_unit              = NULL;
            boxcons->b                   = NULL;
            boxcons->b_unit              = NULL;
            boxcons->c                   = NULL;
            boxcons->c_unit              = NULL;
            boxcons->numfacets           = 0;

            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(boxcons->id), attrvalue, strlen(attrvalue) + 1) );
            tag2 = "boxModelBounds";
            nodeprop = SCIPxmlFindNodeMaxdepth(config, tag2, 0, 1);
            for (prop = SCIPxmlFirstChild(nodeprop); prop != NULL; prop = SCIPxmlNextSibl(prop)) /*lint !e2840*/
            {
               /* looks for xml node "gasTemperature" */
               if ( strcmp(SCIPxmlGetName(prop), "gasTemperature") == 0 )
               {
                  propvalue = SCIPxmlGetAttrval(prop, "value");
                  if ( propvalue == NULL )
                  {
                     SCIPerrorMessage("Value of element \"gasTemperature\" not found.\n");
                     return SCIP_READERROR;
                  }
                  gasTemperature = atof(propvalue);
                  unit = SCIPxmlGetAttrval(prop, "unit");
                  SCIP_CALL( change_temperature_unit(unit, &gasTemperature) );
                  if ( SCIPisEQ(scip, gasTemperature, probdata->network->gasTemperature) )
                  {
                     SCIP_CALL( CSreadConfigBoxData(scip, probdata, prop, arc, boxcons) );
                     equalTemperature = TRUE;
                  }
               }
            }

            if ( ! equalTemperature )
            {
               min = SCIPinfinity(scip);
               for (prop = SCIPxmlFirstChild(nodeprop); prop != NULL; prop = SCIPxmlNextSibl(prop)) /*lint !e2840*/
               {
                  /* looks for xml node "gasTemperature" */
                  if ( strcmp(SCIPxmlGetName(prop), "gasTemperature") == 0 )
                  {
                     propvalue = SCIPxmlGetAttrval(prop, "value");
                     gasTemperature = atof(propvalue);
                     unit = SCIPxmlGetAttrval(prop, "unit");
                     SCIP_CALL( change_temperature_unit(unit, &gasTemperature) );
                     tempdiff = fabs(gasTemperature - probdata->network->gasTemperature);
                     if ( tempdiff < min )
                     {
                        closestTemperature = gasTemperature;
                        min = tempdiff;
                     }
                     SCIPdebugMessage("tempdiff : <%f> \n", tempdiff);
                  }
               }

               for (prop = SCIPxmlFirstChild(nodeprop); prop != NULL; prop = SCIPxmlNextSibl(prop)) /*lint !e2840*/
               {
                  /* looks for xml node "gasTemperature" */
                  if ( strcmp(SCIPxmlGetName(prop), "gasTemperature") == 0 )
                  {
                     propvalue = SCIPxmlGetAttrval(prop, "value");
                     gasTemperature = atof(propvalue);
                     unit = SCIPxmlGetAttrval(prop, "unit");
                     SCIP_CALL( change_temperature_unit(unit, &gasTemperature) );
                     if ( SCIPisEQ(scip, gasTemperature, closestTemperature) )
                     {
                        SCIP_CALL( CSreadConfigBoxData(scip, probdata, prop, arc, boxcons) );
                     }
                  }
               }
            }

            ++k;
         }

         SCIPfreeBlockMemoryArray(scip, &(compressorname), strlen(compressorname) + 1);
      }
   }
   SCIPxmlFreeNode(start);

   return SCIP_OKAY;
}

/** read data of a turbo compressor machine from cs-file */
static
SCIP_RETCODE CSreadTurboCSData(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       xml_comp,           /**< XML_node of compressor machine */
   GAS_TurboCS*          turbocs             /**< turbo compressor struct */
   )
{
   const XML_NODE* xml_data;
   const char* value;
   const char* name;
   const char* unit;
   char tag_list[17][20] = {"speedMin", "speedMax", "n_isoline_coeff_1", "n_isoline_coeff_2", "n_isoline_coeff_3", "n_isoline_coeff_4", "n_isoline_coeff_5", "n_isoline_coeff_6", "n_isoline_coeff_7", "n_isoline_coeff_8", "n_isoline_coeff_9", "surgeline_coeff_1", "surgeline_coeff_2", "surgeline_coeff_3", "chokeline_coeff_1", "chokeline_coeff_2", "chokeline_coeff_3"};
   int i;

   assert( scip != NULL );
   assert( xml_comp != NULL );
   assert( turbocs != NULL );
   assert( strcmp(SCIPxmlGetName(xml_comp), "turboCompressor") == 0 );

   name = SCIPxmlGetAttrval(xml_comp, "id");
   if ( name == NULL )
   {
      SCIPerrorMessage("Name of turboCompressor not found.\n");
      return SCIP_READERROR;
   }
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(turbocs->id), name, strlen(name) + 1) );

   xml_data = SCIPxmlFirstChild(xml_comp);
   if ( xml_data == NULL )
   {
      SCIPerrorMessage("turboCompressor <%s> has no xml_children.\n", name);
      return SCIP_READERROR;
   }

   /* find properties of turbo compressor */
   for ( i = 0; i < 17; ++i )
   {
      xml_data = SCIPxmlFindNodeMaxdepth(xml_comp, tag_list[i], 0, 1);
      if ( xml_data == NULL )
      {
         SCIPerrorMessage("Could not find <%s> of compressor <%s>.\n", tag_list[i], turbocs->id);
         return SCIP_READERROR;
      }

      value = SCIPxmlGetAttrval(xml_data, "value");
      if ( value == NULL )
      {
         SCIPerrorMessage("Could not find value of <%s> of compressor <%s>.\n", tag_list[i], turbocs->id);
         return SCIP_READERROR;
      }

      if ( (strcmp(tag_list[i], "speedMin") == 0) || (strcmp(tag_list[i], "speedMax") == 0) )
      {
         unit = SCIPxmlGetAttrval(xml_data, "unit");
         if ( unit == NULL )
         {
            SCIPerrorMessage("Could not find unit of <%s> of compressor <%s>.\n", tag_list[i], turbocs->id);
            return SCIP_READERROR;
         }
         else if ( strcmp(unit, "per_min") != 0 )
         {
            SCIPerrorMessage("Speed of compressor machine <%s> not given as <per_min>, have to implement conversion.", turbocs->id);
            return SCIP_READERROR;
         }

         if ( strcmp(tag_list[i], "speedMin") == 0 )
            turbocs->speedMin = atof(value);
         else
            turbocs->speedMax = atof(value);
      }
      else if ( strcmp(tag_list[i], "n_isoline_coeff_1") == 0 )
         turbocs->speedMatrix[0] = atof(value);
      else if ( strcmp(tag_list[i], "n_isoline_coeff_2") == 0 )
         turbocs->speedMatrix[1] = atof(value);
      else if ( strcmp(tag_list[i], "n_isoline_coeff_3") == 0 )
         turbocs->speedMatrix[2] = atof(value);
      else if ( strcmp(tag_list[i], "n_isoline_coeff_4") == 0 )
         turbocs->speedMatrix[3] = atof(value);
      else if ( strcmp(tag_list[i], "n_isoline_coeff_5") == 0 )
         turbocs->speedMatrix[4] = atof(value);
      else if ( strcmp(tag_list[i], "n_isoline_coeff_6") == 0 )
         turbocs->speedMatrix[5] = atof(value);
      else if ( strcmp(tag_list[i], "n_isoline_coeff_7") == 0 )
         turbocs->speedMatrix[6] = atof(value);
      else if ( strcmp(tag_list[i], "n_isoline_coeff_8") == 0 )
         turbocs->speedMatrix[7] = atof(value);
      else if ( strcmp(tag_list[i], "n_isoline_coeff_9") == 0 )
         turbocs->speedMatrix[8] = atof(value);
      else if ( strcmp(tag_list[i], "surgeline_coeff_1") == 0 )
         turbocs->surgeline[0] = atof(value);
      else if ( strcmp(tag_list[i], "surgeline_coeff_2") == 0 )
         turbocs->surgeline[1] = atof(value);
      else if ( strcmp(tag_list[i], "surgeline_coeff_3") == 0 )
         turbocs->surgeline[2] = atof(value);
      else if ( strcmp(tag_list[i], "chokeline_coeff_1") == 0 )
         turbocs->chokeline[0] = atof(value);
      else if ( strcmp(tag_list[i], "chokeline_coeff_2") == 0 )
         turbocs->chokeline[1] = atof(value);
      else if ( strcmp(tag_list[i], "chokeline_coeff_3") == 0 )
         turbocs->chokeline[2] = atof(value);
   }

   return SCIP_OKAY;
}

/** read a configuration from cs-file */
static
SCIP_RETCODE CSreadConfiguration(
   SCIP*                 scip,               /**< SCIP data structure */
   const XML_NODE*       xml_config,         /**< XML_node of a configuration */
   GAS_CSConfig*         config              /**< configuration struct */
   )
{
   const XML_NODE* xml_stage;
   const XML_NODE* xml_cs;
   const char* value;
   const char* id_cf;
   const char* id_cs;

   assert( scip != NULL );
   assert( xml_config != NULL );
   assert( config != NULL );
   assert( strcmp(SCIPxmlGetName(xml_config), "configuration") == 0 );

   /* initialize config struct */
   config->id = NULL;
   config->numSerialStages = 0;
   config->numCsStageOne = 0;
   config->stageOneCsOne = NULL;
   config->stageOneCsTwo = NULL;

   id_cf = SCIPxmlGetAttrval(xml_config, "confId");
   if ( id_cf == NULL )
   {
      SCIPerrorMessage("Name of configuration not found.\n");
      return SCIP_READERROR;
   }
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(config->id), id_cf, strlen(id_cf) + 1) );

   value = SCIPxmlGetAttrval(xml_config, "nrOfSerialStages");
   if ( value == NULL )
   {
      SCIPerrorMessage("Number of serial stages of configuration not found.\n");
      return SCIP_READERROR;
   }
   config->numSerialStages = atoi(value);

   if ( config->numSerialStages > 1 )
   {
      SCIPerrorMessage("Currently configurations with serial stages are not supported.\n");
      return SCIP_READERROR;
   }
   else if ( config->numSerialStages == 0 )
   {
      SCIPerrorMessage("Configuration <%s> must have at least one stage!\n", id_cf);
      return SCIP_READERROR;
   }

   /* find first stage of configuration */
   xml_stage = SCIPxmlFindNodeMaxdepth(xml_config, "stage", 0, 1);
   if ( xml_stage == NULL )
   {
      SCIPerrorMessage("No stage of configuration <%s> found.\n", id_cf);
      return SCIP_READERROR;
   }

   value = SCIPxmlGetAttrval(xml_stage, "nrOfParallelUnits");
   if ( value == NULL )
   {
      SCIPerrorMessage("Number of parallel units of configuration <%s> not found.\n", id_cf);
      return SCIP_READERROR;
   }
   config->numCsStageOne = atoi(value);

   if ( config->numCsStageOne > 2 )
   {
      SCIPerrorMessage("Configuration <%s> should have at most two parallel units.\n", id_cf);
      return SCIP_READERROR;
   }
   else if ( config->numCsStageOne == 0 )
   {
      SCIPerrorMessage("Configuration <%s> must have at least one compressor in its stage.\n", id_cf);
      return SCIP_READERROR;
   }

   /* find first compressor in stage */
   xml_cs = SCIPxmlFirstChild(xml_stage);
   if ( xml_cs == NULL )
   {
      SCIPerrorMessage("First xml_stage of configuration <%s> does not have children.\n", id_cf);
      return SCIP_READERROR;
   }
   else if ( strcmp(SCIPxmlGetName(xml_cs), "compressor") != 0 )
   {
      SCIPerrorMessage("Stage number <%s> of config <%s> has unknown child node <%s>.\n", SCIPxmlGetAttrval(xml_stage, "stageNr"), id_cf, SCIPxmlGetName(xml_cs));
      return SCIP_READERROR;
   }

   id_cs = SCIPxmlGetAttrval(xml_cs, "id");
   if ( id_cs == NULL )
   {
      SCIPerrorMessage("Name of compressor of config <%s> not found.\n", id_cf);
      return SCIP_READERROR;
   }
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(config->stageOneCsOne), id_cs, strlen(id_cs) + 1) );

   /* find second compressor in stage */
   if ( config->numCsStageOne == 2 )
   {
      xml_cs = SCIPxmlNextSibl(xml_cs);
      if ( xml_cs == NULL )
      {
         SCIPerrorMessage("First xml_stage of configuration <%s> should contain two compressors.\n", id_cf);
         return SCIP_READERROR;
      }
      else if ( strcmp(SCIPxmlGetName(xml_cs), "compressor") != 0 )
      {
         SCIPerrorMessage("Stage number <%s> of config <%s> has unknown child node <%s>.\n", SCIPxmlGetAttrval(xml_stage, "stageNr"), id_cf, SCIPxmlGetName(xml_cs));
         return SCIP_READERROR;
      }

      id_cs = SCIPxmlGetAttrval(xml_cs, "id");
      if ( id_cs == NULL )
      {
         SCIPerrorMessage("Name of second compressor of config <%s> not found.\n", id_cf);
         return SCIP_READERROR;
      }
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(config->stageOneCsTwo), id_cs, strlen(id_cs) + 1) );
   }

   return SCIP_OKAY;
}

/** read data for detailed compressor station model
 *  from compressor station (CS) file
 */
static
SCIP_RETCODE CSreadDetailedCSData(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of file to read */
   SCIP_PROBDATA*        probdata            /**< problem data to be filled */
   )
{
   XML_NODE* start;                          /* points to first xml_node in Gas file */
   const XML_NODE* xml_css;
   const XML_NODE* xml_cs;                   /* compressor station node */
   const XML_NODE* xml_compressors;
   const XML_NODE* xml_comp;                 /* compressor node */
   const XML_NODE* xml_configurations;
   const XML_NODE* xml_config;               /* config node */
   GAS_Arc* arc;
   GAS_CS* cs_struct;
   const char* attrname;
   int counter_tc;
   int counter_cf;

   assert( scip != NULL );
   assert( filename != NULL );
   assert( probdata != NULL );
   assert( probdata->network != NULL );

   /* read xml file */
   start = SCIPxmlProcess(filename);
   if ( start == NULL )
   {
      SCIPerrorMessage("Some error occured during parsing <%s>.\n", filename);
      return SCIP_READERROR;
   }

   /* find nodes */
   xml_css = SCIPxmlFindNodeMaxdepth(start, "compressorStations", 0, 1);
   if ( xml_css == NULL )
   {
      /* free xml data */
      SCIPxmlFreeNode(start);

      SCIPerrorMessage("CompressorStations not found.\n");
      return SCIP_READERROR;
   }

   /* loop through all compressor stations */
   for (xml_cs = SCIPxmlFirstChild(xml_css); xml_cs != NULL; xml_cs = SCIPxmlNextSibl(xml_cs)) /*lint !e2840*/
   {
      /* check if the child of "compressorStations" is a compressor station. */
      if ( strcmp(SCIPxmlGetName(xml_cs), "compressorStation") != 0 )
         continue;

      attrname = SCIPxmlGetAttrval(xml_cs, "id");
      if ( attrname == NULL )
      {
         /* free xml data */
         SCIPxmlFreeNode(start);

         SCIPerrorMessage("Attribute \"id\" in cs file not found.\n");
         return SCIP_READERROR;
      }

      /* get cs arc */
      arc = (GAS_Arc*) SCIPhashtableRetrieve(probdata->network->Arcshashtable, (void*) attrname);
      if ( arc == NULL )
      {
         /* free xml data */
         SCIPxmlFreeNode(start);

         SCIPerrorMessage("CompressorStation <%s> does not exist!\n", attrname);
         return SCIP_READERROR;
      }

      cs_struct = arc->detailed_info;

      /* count compressors */
      xml_compressors = SCIPxmlFindNodeMaxdepth(xml_cs, "compressors", 0, 1);

      if ( xml_compressors == NULL )
      {
         SCIPerrorMessage("No compressor in compressor station <%s> found.\n", attrname);
         return SCIP_READERROR;
      }

      /* count number of compressors (and types)  */
      for (xml_comp = SCIPxmlFirstChild(xml_compressors); xml_comp != NULL; xml_comp = SCIPxmlNextSibl(xml_comp)) /*lint !e2840*/
      {
         if ( strcmp(SCIPxmlGetName(xml_comp), "turboCompressor") == 0 )
            ++(cs_struct->numTurboCS);
         else if ( strcmp(SCIPxmlGetName(xml_comp), "pistonCompressor") == 0 )
            ++(cs_struct->numPistonCS);
         else
         {
            SCIPerrorMessage("Found unkown compressor type <%s>.\n", SCIPxmlGetName(xml_comp));
            /* free xml data */
            SCIPxmlFreeNode(start);

            return SCIP_READERROR;
         }
      }

      /* allocate memory for compressor structs */
      if ( cs_struct->numTurboCS > 0 )
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &cs_struct->turbocs, cs_struct->numTurboCS) );
      }
      else if ( cs_struct->numPistonCS > 0 )
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &cs_struct->pistoncs, cs_struct->numPistonCS) );
      }

      /* read data of compressor machines */
      counter_tc = 0;
      for (xml_comp = SCIPxmlFirstChild(xml_compressors); xml_comp != NULL; xml_comp = SCIPxmlNextSibl(xml_comp)) /*lint !e2840*/
      {
         if ( strcmp(SCIPxmlGetName(xml_comp), "turboCompressor") == 0 )
         {
            SCIP_CALL( CSreadTurboCSData(scip, xml_comp, &(cs_struct->turbocs[counter_tc])) );
            ++counter_tc;
         }
         else if ( strcmp(SCIPxmlGetName(xml_comp), "pistonCompressor") == 0 )
         {
            /* SCIP_CALL( CSreadPistonCSData(scip, xml_comp, &(cs_struct->pistoncs[counter_pc])) ); */
            SCIPerrorMessage("Piston compressors are currently not supported.\n");
            return SCIP_READERROR;
         }
      }

      /* read configurations */
      xml_configurations = SCIPxmlFindNodeMaxdepth(xml_cs, "configurations", 0, 1);

      if ( xml_configurations == NULL )
      {
         SCIPerrorMessage("Compressor station <%s> does not have a configuration section.\n", SCIPxmlGetName(xml_cs));
         return SCIP_READERROR;
      }

      for (xml_config = SCIPxmlFirstChild(xml_configurations); xml_config != NULL; xml_config = SCIPxmlNextSibl(xml_config)) /*lint !e2840*/
      {
         if ( strcmp(SCIPxmlGetName(xml_config), "configuration") == 0 )
         {
            ++(cs_struct->numconfigurations);
         }
         else
         {
            SCIPerrorMessage("Found unkown child node <%s> of configurations of compressor station <%s>.\n", SCIPxmlGetName(xml_config), SCIPxmlGetName(xml_cs));
            return SCIP_READERROR;
         }
      }

      if ( cs_struct->numconfigurations == 0 )
      {
         SCIPerrorMessage("Compressor station <%s> does not have a configuration.\n", SCIPxmlGetName(xml_cs));
         return SCIP_READERROR;
      }

      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &cs_struct->configurations, cs_struct->numconfigurations) );

      counter_cf = 0;
      for (xml_config = SCIPxmlFirstChild(xml_configurations); xml_config != NULL; xml_config = SCIPxmlNextSibl(xml_config)) /*lint !e2840*/
      {
         SCIP_CALL( CSreadConfiguration(scip, xml_config, &(cs_struct->configurations[counter_cf])) );
         ++counter_cf;
      }
   }

   /* free xml data */
   SCIPxmlFreeNode(start);

   return SCIP_OKAY;
}

/** read compressor station (CS) file */
static
SCIP_RETCODE CSreadFile(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< name of file to read */
   SCIP_PROBDATA*        probdata            /**< problem data to be filled */
   )
{
   assert( scip != NULL );
   assert( filename != NULL );
   assert( probdata != NULL );

   if ( probdata->boxConstraintModel )
   {
      SCIP_CALL( CSreadBoxModelData(scip, filename, probdata) );
   }
   else if ( probdata->charDiagram )
   {
      SCIP_CALL( CSreadDetailedCSData(scip, filename, probdata) );
   }
   else
   {
      SCIPerrorMessage("There is a cs-file given, but no model should be used:\n <%s>\n", filename);
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** function that updates the spanning tree structe, i.e., the wayToNode matrix and the arcsInTree array */
static
SCIP_RETCODE updateSpanningTreeStructure(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network,            /**< network data */
   int                   parentNode,         /**< position of the parent node */
   int                   childNode,          /**< position of the child node */
   int                   arcposition,        /**< position of the arc connecting parent and child node */
   int                   orientation         /**< orientation of the arc such that it points from parent to child */
   )
{
   int i;

   assert( scip != NULL );
   assert( network != NULL );
   assert( (orientation == 1) || (orientation == -1) );

   /* update the way to node matrix */
   for (i = 0; i < network->numarcs; ++i)
   {
      (network->spanningTree->wayToNode)[i][childNode] = (network->spanningTree->wayToNode)[i][parentNode];
   }
   (network->spanningTree->wayToNode)[arcposition][childNode] = orientation;

   /* set the entry in arcsInTree */
   (network->spanningTree->arcsInTree)[arcposition] = 1;

   return SCIP_OKAY;
}


/** checks for the in- and outarcs of a node, whether or not the target and
 *  source node are already in a given hashtable, adds the target or source
 *  node to the the table and queue if not and calls the routine to update
 *  the tree
 */
static
SCIP_RETCODE addNeighborsToSpanningTree(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Node*             currentNode,        /**< node to check for neighbors not in the tree */
   GAS_Network*          network,            /**< network data */
   SCIP_HASHTABLE**      nodesInTree,        /**< nodes already in the tree */
   SCIP_QUEUE**          nodesToCheck,       /**< nodes not checked for neighbors not in the tree */
   int*                  counter             /**< counter for the nodes in the tree */
   )
{
   GAS_Arc*  arc;
   GAS_Node* neighbor;

   assert( scip != NULL );
   assert( currentNode != NULL );
   assert( network != NULL );
   assert( nodesInTree != NULL );
   assert( nodesToCheck != NULL );
   assert( counter != NULL );

   /* start with checking the inarcs */
   arc = currentNode->inarcs;
   while ( arc != NULL )
   {
      neighbor = arc->sourcenode;

      if ( SCIPhashtableRetrieve(*nodesInTree, (void*) (neighbor->id)) == NULL )
      {
         ++(*counter);
         SCIP_CALL( SCIPqueueInsert(*nodesToCheck, (void*) neighbor) );
         SCIP_CALL( SCIPhashtableInsert(*nodesInTree, (void*) neighbor) );
         SCIP_CALL( updateSpanningTreeStructure(scip, network, currentNode->nodeposition, neighbor->nodeposition, arc->arcposition, -1) );
      }

      arc = arc->next_inarc;
   }

   /* next we check the outarcs */
   arc = currentNode->outarcs;
   while ( arc != NULL )
   {
      neighbor = arc->targetnode;

      if ( SCIPhashtableRetrieve(*nodesInTree, (void*) (neighbor->id)) == NULL )
      {
         ++(*counter);
         SCIP_CALL( SCIPqueueInsert(*nodesToCheck, (void*) neighbor) );
         SCIP_CALL( SCIPhashtableInsert(*nodesInTree, (void*) neighbor) );
         SCIP_CALL( updateSpanningTreeStructure(scip, network, currentNode->nodeposition, neighbor->nodeposition, arc->arcposition, 1) );
      }

      arc = arc->next_outarc;
   }

   return SCIP_OKAY;
}

/** compute a spanning tree for finding a cycle basis
 *
 *  Returns a matrix wayToNode with rows corresponding to the arcs and columns corresponding to the nodes and
 *  furthermore an array which indicates the arcs in the tree.
 *
 *  @note: We assume that the network is still connected when removing the compressors
 */
static
SCIP_RETCODE computeSpanningTree(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network             /**< network data */
   )
{
   GAS_SpanningTree* spanningTree;
   GAS_Node*         root_node;
   GAS_Node*         nodeToCheck;
   SCIP_HASHTABLE*   nist;                   /* hashtable with nodes in the tree */
   SCIP_QUEUE*       nVisitedNodes;          /* nodes in the tree where out/ingoing arcs are not checked */
   int               counter;                /* counts the number of nodes in the tree */
   int**             wayToNode;              /* matrix with entries -1,0,1 which indicate whether or not an arc is on the way from the root node to the specific node; -1/1 for the orientation */
   int*              arcsInTree;             /* array of 0/1 for the arcs in the spanning tree */
   int               i;

   assert( scip != NULL );
   assert( network != NULL );

   /* get memory for the tree and initialize it */
   SCIP_CALL( SCIPallocBlockMemory(scip, &spanningTree) );
   spanningTree->root_node    = &(network->nodes_ptr[0]);
   root_node                  = &(network->nodes_ptr[0]);
   spanningTree->foundCycles = FALSE;
   spanningTree->wayToNode = NULL;
   spanningTree->arcsInTree = NULL;
   network->spanningTree = spanningTree;

   if ( network->numnodes > (network->numarcs) )
   {
      SCIPdebugMessage("Network has to few arcs to have cycles or is not connected without the compressors.\n");
      return SCIP_OKAY;
   }

   /* get memory for the wayToNode Matrix */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &wayToNode, network->numarcs) );
   for ( i = 0; i< network->numarcs; ++i )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(wayToNode[i]), network->numnodes) );
   }

   /* get memory for the arcsInTree array */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &arcsInTree, network->numarcs) );

   /* set the pointer in the spanning tree struct */
   spanningTree->wayToNode = wayToNode;
   spanningTree->arcsInTree = arcsInTree;

   /* create queue for the nodes to look at */
   SCIP_CALL( SCIPqueueCreate( &nVisitedNodes, 10, 2.0) );

   /* create hashtable for the nodes in the tree */
#if SCIP_VERSION >= 400
   SCIP_CALL( SCIPhashtableCreate(&nist, SCIPblkmem(scip), 1000, SCIPhashGetKeyNode, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );
#else
   SCIP_CALL( SCIPhashtableCreate(&nist, SCIPblkmem(scip), SCIPcalcHashtableSize(1000), SCIPhashGetKeyNode, SCIPhashKeyEqString, SCIPhashKeyValString, NULL) );
#endif

   /* add the first node to the queue and the hashtable */
   SCIP_CALL( SCIPqueueInsert(nVisitedNodes, (void*) root_node) );
   SCIP_CALL( SCIPhashtableInsert(nist, (void*) root_node) );

   counter = 1;
   while ( counter < network->numnodes )
   {
      if ( SCIPqueueIsEmpty(nVisitedNodes) )
      {
         /* Error: Something went wrong or the network is not connected */
         SCIPdebugMessage("Empty queue of unchecked nodes while trying to find a spanning tree.\n");
         break;
      }
      else
      {
         /* pick the first node in the queue and remove it */
         nodeToCheck = (GAS_Node*) SCIPqueueRemove(nVisitedNodes);
      }
      /* check the neighbors and add them to the tree */
      SCIP_CALL( addNeighborsToSpanningTree(scip, nodeToCheck, network, &nist, &nVisitedNodes, &counter) );
   }

   if ( counter == network->numnodes )
   {
      network->spanningTree->foundCycles = TRUE;
      /* printf("\nwayToNode Matrix of the network:\n"); */
      /* for ( i = 0; i < network->numarcs; ++i ) */
      /* { */
      /*    int j; */
      /*    for ( j = 0; j < network->numnodes; ++j ) */
      /*    { */
      /*       printf("%d ", spanningTree->wayToNode[i][j]); */
      /*    } */
      /*    printf("\n"); */
      /* } */
   }
   else if ( counter > network->numnodes )
   {
      SCIPerrorMessage("Error in computeSpanningTree: added too many nodes to the tree!\n");
   }

   /* we don't need the queue and hashtable anymore */
   SCIPqueueFree(&(nVisitedNodes));
   SCIPhashtableFree(&(nist));

   return SCIP_OKAY;
}

/** function, which finds the fundamental cycles in the network */
static
SCIP_RETCODE findFundamentalCycles(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network             /**< network data */
   )
{
   GAS_SpanningTree* tree;
   GAS_Controlvalve* cv;
   GAS_CS*           cs;
   GAS_Cycle*        newCycle;
   int               sourceposition;
   int               targetposition;
   int               i;
   int               j;

   assert( scip != NULL );
   assert( network != NULL);
   assert( network->spanningTree != NULL );
   assert( network->spanningTree->wayToNode != NULL );
   assert( network->spanningTree->arcsInTree != NULL );

   tree = network->spanningTree;

   for (i = 0; i < network->numarcs; ++i)
   {
      /* if arc i already is in the tree it cannot generate a cycles */
      if ( tree->arcsInTree[i] == 1 )
         continue;

      /* if a compressor or controlvalve have a bypass, use that for the cycle */
      if ( network->arcs_ptr[i].type == CS )
      {
         cs =  network->arcs_ptr[i].detailed_info;
         if ( cs->internalBypassRequired )
         {
            continue;
         }
      }
      else if ( network->arcs_ptr[i].type == CONTROLVALVE )
      {
         cv =  network->arcs_ptr[i].detailed_info;
         if ( cv->internalBypassRequired )
            continue;
      }

      ++(network->numBasisCycles);
      SCIP_CALL( SCIPallocBlockMemory(scip, &newCycle) );
      newCycle->combinedBasisCycles = NULL;
      newCycle->arcsInCycle = NULL;
      newCycle->next_cycle = NULL;
      newCycle->nr_cycle = network->numBasisCycles;
      newCycle->previous_cycle = network->lastCycle;
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &(newCycle->arcsInCycle), network->numarcs) );

      sourceposition = network->arcs_ptr[i].sourcenode->nodeposition;
      targetposition = network->arcs_ptr[i].targetnode->nodeposition;

      for (j = 0; j < network->numarcs; ++j)
      {
         if ( j == i )
            newCycle->arcsInCycle[j] = 1;
         else
            newCycle->arcsInCycle[j] = tree->wayToNode[j][sourceposition] - tree->wayToNode[j][targetposition];
      }

      if ( network->firstCycle == NULL )
         network->firstCycle = newCycle;

      if ( network->lastCycle != NULL )
         network->lastCycle->next_cycle =  newCycle;

      network->lastCycle = newCycle;
   }

   return SCIP_OKAY;
}

/** function to check whether two cycles can and should be combined.
 *
 *  If they should, then the function stores the combination of the basis cycle and the
 *  second cycle in the array arcsincycle.
 *  The first cycle given to this function has to be a basis cycle!
 */
static
SCIP_RETCODE checkPossibleCombination(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network,            /**< network data */
   GAS_Cycle*            basisCycle,         /**< cycle to possibly combine with the second cycle */
   GAS_Cycle*            secondCycle,        /**< cycle to possibly combine with the first cycle */
   int**                 arcsincycle,        /**< pointer to array to store the arcs of a new cycle */
   int**                 combinedcycles,     /**< pointer to array to store the combined basis cycles */
   SCIP_Bool*            combinable          /**< Boolean to store if cycles can be combined  */
   )
{  /*lint --e{715}*/
   GAS_Cycle* cmb_cycle;
   GAS_Node* source;
   GAS_Node* current_node;
   GAS_Node* next_node;
   int arc_num;
   GAS_Arc* arc;
   int  sign = 0;
   int  i;

   assert( scip != NULL );
   assert( arcsincycle != NULL );
   assert( combinedcycles != NULL );
   assert( combinable != NULL );

   *combinable = TRUE;

   /* make sure arcsincycle does not contain old information */
   for (i = 0; i < network->numarcs; ++i)
      (*arcsincycle)[i] = 0;

   if ( basisCycle->combinedBasisCycles != NULL )
   {
      *combinable = FALSE;
      SCIPerrorMessage("Function checkPossibleCombination called with a non basisCycle!\n");
      return SCIP_READERROR;
   }

   if ( secondCycle->combinedBasisCycles != NULL )
   {
      /* check if first cycle already is part of second cycle */
      if ( secondCycle->combinedBasisCycles[basisCycle->nr_cycle - 1] == 1 )
      {
         *combinable = FALSE;
         return SCIP_OKAY;
      }
   }
   else if ( secondCycle->combinedBasisCycles == NULL )
   {
      if ( basisCycle->nr_cycle == secondCycle->nr_cycle )
      {
         *combinable = FALSE;
         return SCIP_OKAY;
      }
   }

   if ( secondCycle->combinedBasisCycles == NULL )
   {
      for (i = 0; i < network->numBasisCycles; ++i)
      {
         (*combinedcycles)[i] = 0;
      }
      (*combinedcycles)[secondCycle->nr_cycle -1] = 1;
   }
   else
   {
      for (i = 0; i < network->numBasisCycles; ++i)
      {
         (*combinedcycles)[i] = secondCycle->combinedBasisCycles[i];
      }
   }
   (*combinedcycles)[basisCycle->nr_cycle -1] = 1;

   /* check if combination already exists */
   cmb_cycle = network->firstCombinedCycle;
   while ( cmb_cycle != NULL )
   {
      for (i = 0; i < network->numBasisCycles; ++i)
      {
         if ( cmb_cycle->combinedBasisCycles[i] != (*combinedcycles)[i] )
            break;
      }

      if ( i == network->numBasisCycles )
      {
         *combinable = FALSE;
         break;
      }

      cmb_cycle = cmb_cycle->next_cycle;
   }

   /* try to write arcsInTree of the new cycle as
    *    basisCycle->arcsInCycle + sign * secondCycle->arcsInCycle
    * therefore set sign to 1 if they use an arc in oppsite direction,
    * otherwise if they use it in the same direction set sign to -1
    * if they use multiple arcs in different directions then the cycles are not *combinable
    */
   if ( *combinable )
   {
      for (i = 0; i < network->numarcs; ++i)
      {
         if ( basisCycle->arcsInCycle[i] == 1 )
         {
            if ( secondCycle->arcsInCycle[i] == 1 )
            {
               if ( sign == 1 )
               {
                  *combinable = FALSE;
                  break;
               }
               sign = -1;
            }
            else if ( secondCycle->arcsInCycle[i] == -1 )
            {
               if ( sign == -1 )
               {
                  *combinable = FALSE;
                  break;
               }
               sign = 1;
            }
         }
         else if ( basisCycle->arcsInCycle[i] == - 1 )
         {
            if ( secondCycle->arcsInCycle[i] == -1 )
            {
               if ( sign == 1 )
               {
                  *combinable = FALSE;
                  break;
               }
               sign = -1;
            }
            else if ( secondCycle->arcsInCycle[i] == 1 )
            {
               if ( sign == -1 )
               {
                  *combinable = FALSE;
                  break;
               }
               sign = 1;
            }
         }
      }
   }

   if ( sign == 0 )
      *combinable = FALSE;

   if ( *combinable == FALSE )
      return SCIP_OKAY;

   /* combine cycles */
   for (i = 0; i < network->numarcs; ++i)
   {
      (*arcsincycle)[i] = basisCycle->arcsInCycle[i] + sign * (secondCycle->arcsInCycle[i]);
   }

   /* check if combination is simple and connected */
   /* Therefore find first arc in the cycle */
   arc_num = 0;
   while ( (*arcsincycle)[arc_num] == 0 && arc_num < network->numarcs )
      ++arc_num;

   if ( arc_num == network->numarcs )
   {
      SCIPerrorMessage("Something went completely wrong, combining cycles produced an empty array of arcs.\n");
      return SCIP_READERROR;
   }

   if ( (*arcsincycle)[arc_num] == -1 )
   {
      source = network->arcs_ptr[arc_num].targetnode;
      next_node = network->arcs_ptr[arc_num].sourcenode;
   }
   else if ( (*arcsincycle)[arc_num] == 1 )
   {
      source = network->arcs_ptr[arc_num].sourcenode;
      next_node = network->arcs_ptr[arc_num].targetnode;
   }
   else
   {
      SCIPerrorMessage("Something went completely wrong, when combining cycles. Current arc is %d-times in the cycle!\n", (*arcsincycle)[arc_num]);
      return SCIP_READERROR;
   }

   /* temporarily remove arc from cycle */
   (*arcsincycle)[arc_num] = 0;

   /* check if the other arcs (with the directions given in the array arcsincycle)
    * form a path from next_node to the source */
   while ( next_node != source )
   {
      current_node = next_node;
      next_node = NULL;

      /* check if inarcs of current node are in the cycle */
      arc = current_node->inarcs;
      while ( arc != NULL )
      {
         arc_num = arc->arcposition;

         if ( (*arcsincycle)[arc_num] == -1 )
         {
            if ( next_node == NULL )
            {
               next_node = arc->sourcenode;
               (*arcsincycle)[arc_num] = 0;
            }
            else
            {
               /* Combination of cycles is not simple! current_node has more than 2 incident arcs */
               *combinable = FALSE;
               break;
            }
         }
         else if ( (*arcsincycle)[arc_num] == 1 )
         {
            /* Combination of cycles is not simple or not oriented correctly!
             * current_node has at least 2 inarcs */
            *combinable = FALSE;
            break;
         }
         else if ( (*arcsincycle)[arc_num] != 0 )
         {
            SCIPerrorMessage("Something went completely wrong, when combining cycles. Current arc is %d-times in the cycle!\n", (*arcsincycle)[arc_num]);
            return SCIP_READERROR;
         }

         arc = arc->next_inarc;
      }

      if ( *combinable == FALSE )
         break;

      /* check if outarcs of current node are in the combination */
      arc = current_node->outarcs;
      while ( arc != NULL )
      {
         arc_num = arc->arcposition;

         if ( (*arcsincycle)[arc_num] == 1 )
         {
            if ( next_node == NULL )
            {
               next_node = arc->targetnode;
               (*arcsincycle)[arc_num] = 0;
            }
            else
            {
               /* Combination of cycles is not simple! current_node has more than 2 incident arcs */
               *combinable = FALSE;
               break;
            }
         }
         else if ( (*arcsincycle)[arc_num] == -1 )
         {
            /* Combination of cycles is not simple or not oriented correctly!
             * current_node has at least 2 inarcs */
            *combinable = FALSE;
            break;
         }
         else if ( (*arcsincycle)[arc_num] != 0 )
         {
            SCIPerrorMessage("Something went completely wrong, when combining cycles. Current arc is %d-times in the cycle!\n", (*arcsincycle)[arc_num]);
            return SCIP_READERROR;
         }

         arc = arc->next_outarc;
      }

      if ( *combinable == FALSE )
         break;

      if ( next_node == NULL )
      {
         SCIPerrorMessage("Combination of cycles contains a node of degree 1!\n");
         return SCIP_READERROR;
      }
   }

   if ( *combinable )
   {
      for (i = 0; i < network->numarcs; ++i)
      {
         if( (*arcsincycle)[i] != 0 )
         {
            /* Combination of cycles is not connected any more */
            *combinable = FALSE;
            break;
         }
      }
   }

   /* reset the arcsincycle array */
   if ( *combinable )
   {
      for (i = 0; i < network->numarcs; ++i)
         (*arcsincycle)[i] = basisCycle->arcsInCycle[i] + sign * (secondCycle->arcsInCycle[i]);
   }

   return SCIP_OKAY;
}

/** function to combine fundamental cycles if they share at least one edge */
static
SCIP_RETCODE combineCycles(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network             /**< network data */
   )
{
   GAS_Cycle* basisCycle;
   GAS_Cycle* cycleToCombine = NULL;
   GAS_Cycle* firstNewCycle  = NULL;
   GAS_Cycle* lastNewCycle   = NULL;
   GAS_Cycle* newCycle       = NULL;
   int i;
   int nr_cycles;
   SCIP_Bool combinable;
   int* arcsincycle;
   int* combinedcycles;

   SCIP_CALL( SCIPallocClearBufferArray(scip, &arcsincycle, network->numarcs) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &combinedcycles, network->numBasisCycles) );

   nr_cycles = network->numBasisCycles;
   basisCycle = network->firstCycle;

   for (i = 1; i <= network->numBasisCycles; ++i)
   {
      assert( basisCycle != NULL );

      cycleToCombine = basisCycle->next_cycle;
      while ( cycleToCombine != NULL )
      {
         combinable = TRUE;

         SCIP_CALL( checkPossibleCombination(scip, network, basisCycle, cycleToCombine, &arcsincycle, &combinedcycles, &combinable) );

         if ( combinable )
         {
            ++nr_cycles;
            SCIP_CALL( SCIPallocBlockMemory(scip, &newCycle) );
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(newCycle->combinedBasisCycles), combinedcycles, network->numBasisCycles) );
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(newCycle->arcsInCycle), arcsincycle, network->numarcs) );
            newCycle->next_cycle = NULL;
            newCycle->nr_cycle = nr_cycles;

            if ( firstNewCycle == NULL )
            {
               firstNewCycle = newCycle;
               newCycle->previous_cycle = network->lastCycle;
               lastNewCycle = newCycle;
            }
            else
            {
               assert( lastNewCycle != NULL );
               newCycle->previous_cycle = lastNewCycle;
               lastNewCycle->next_cycle = newCycle;
               lastNewCycle = newCycle;
            }
         }
         cycleToCombine = cycleToCombine->next_cycle;
      } /* end while*/

      if ( firstNewCycle != NULL )
      {
         network->lastCycle->next_cycle = firstNewCycle;
         network->lastCycle = lastNewCycle;
         if ( network->firstCombinedCycle == NULL )
            network->firstCombinedCycle = firstNewCycle;

         firstNewCycle = NULL;
         lastNewCycle = NULL;
      }

      basisCycle = basisCycle->next_cycle;
   }

   SCIPfreeBufferArray(scip, &combinedcycles);
   SCIPfreeBufferArray(scip, &arcsincycle);

   return SCIP_OKAY;
}

/* ----------------- public interface functions ------------------------ */

/** initialize probdata structure */
SCIP_RETCODE GASinitProb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< problem data structure */
   )
{
   assert( scip != NULL );

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );

   /* initialize nomination specific data */
   (*probdata)->network             = NULL;
   (*probdata)->netname             = NULL;
   (*probdata)->scnname             = NULL;

   /* variables */
   (*probdata)->objective           = NULL;
   (*probdata)->flowvars            = NULL;
   (*probdata)->netflowvars         = NULL;
   (*probdata)->positiveFlowBinvars = NULL;
   (*probdata)->negativeFlowBinvars = NULL;
   (*probdata)->pressurevars        = NULL;
   (*probdata)->PIvars              = NULL;
   (*probdata)->pSlackVars          = NULL;
   (*probdata)->PIdiffvars          = NULL;
   (*probdata)->VALVE_binvars       = NULL;
   (*probdata)->CV_binvars          = NULL;
   (*probdata)->NLR_DeltaVars       = NULL;
   (*probdata)->NLR_AbsDeltaVars    = NULL;
   (*probdata)->NLR_AbsFlowVars     = NULL;
   (*probdata)->LR_smoothingFlow    = NULL;
   (*probdata)->LR_posFlowDir       = NULL;
   (*probdata)->LR_negFlowDir       = NULL;
   (*probdata)->CS_binvars          = NULL;
   (*probdata)->CS_resistorvars     = NULL;
   (*probdata)->BCM_flowCS          = NULL;
   (*probdata)->BCM_pressureIn      = NULL;
   (*probdata)->BCM_pressureOut     = NULL;
   (*probdata)->Lambda              = NULL;
   (*probdata)->mixingRatio         = NULL;
   (*probdata)->mixingRatioMass     = NULL;
   (*probdata)->nodeMixingRatio     = NULL;
   (*probdata)->nodeMixingRatioMol  = NULL;

   /* model options */
   (*probdata)->mixing              = FALSE;
   (*probdata)->nodeMu              = FALSE;
   (*probdata)->mol                 = FALSE;
   (*probdata)->linearMixing        = FALSE;
   (*probdata)->flowConsMixing      = FALSE;
   (*probdata)->flowConsMixingNode  = FALSE;
   (*probdata)->approx              = FALSE;
   (*probdata)->noFlowBinvars       = FALSE;
   (*probdata)->binaryFlowCons      = TRUE;
   (*probdata)->noCycles            = FALSE;
   (*probdata)->allCycles           = FALSE;
   (*probdata)->algebraic           = FALSE;
   (*probdata)->binVarEQone         = FALSE;
   (*probdata)->papay               = FALSE;
   (*probdata)->boxConstraintModel  = FALSE;
   (*probdata)->additionalFacets    = FALSE;
   (*probdata)->charDiagram         = FALSE;
   (*probdata)->relaxLowerBounds    = FALSE;
   (*probdata)->relaxUpperBounds    = FALSE;
   (*probdata)->relaxCVbounds       = FALSE;
   (*probdata)->relaxLBvalue        = 0.0;
   (*probdata)->relaxUBvalue        = 0.0;
   (*probdata)->relaxCVvalue        = 0.0;
   (*probdata)->resistorCVInternal  = FALSE;
   (*probdata)->idealCSminIncrease  = 1.0;
   (*probdata)->idealCSmaxIncrease  = 50.0;
   (*probdata)->csusemaxincrease    = TRUE;

   /* objective functions */
   (*probdata)->minCompSum          = FALSE;
   (*probdata)->minCompInc          = FALSE;
   (*probdata)->minPressSum         = FALSE;
   (*probdata)->noObjective         = FALSE;
   (*probdata)->powerLoss           = FALSE;
   (*probdata)->minSlack            = FALSE;
   (*probdata)->minSlackPerBound    = FALSE;
   (*probdata)->maxFlow             = FALSE;
   (*probdata)->minLambda           = FALSE;
   (*probdata)->addComponentCuts    = TRUE;
   (*probdata)->addParallelCSCuts   = FALSE;
   (*probdata)->reduceFlowConservation = FALSE;

   /* add parameters to set model options via a settings file
    * default values are those set above
    */
   SCIP_CALL( SCIPaddBoolParam(scip, "noFlowBinvars",
         "whether no binary variables for flowdirections should be used",
         &((*probdata)->noFlowBinvars), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "binaryFlowCons",
         "whether binary flow conservation constraints should be used",
         &((*probdata)->binaryFlowCons), FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "noCycles",
         "whether no cycles should be used for noCircularFlowConstraints",
         &((*probdata)->noCycles), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "allCycles",
         "whether all cycles in the graph should be used for noCircularFlowConstraints",
         &((*probdata)->allCycles), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "algebraic",
         "whether algebraic model should be used",
         &((*probdata)->algebraic), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "binVarEQone",
         "whether flowDirBinvar cons should be equal to one",
         &((*probdata)->binVarEQone), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "mixing",
         "whether mixing model should be used",
         &((*probdata)->mixing), FALSE, FALSE, NULL, NULL ) );

   SCIP_CALL( SCIPaddBoolParam(scip, "linearMixing",
         "whether linear mixing model should be used",
         &((*probdata)->linearMixing), FALSE, FALSE, NULL, NULL ) );

   SCIP_CALL( SCIPaddBoolParam(scip, "flowConsMixing",
         "whether the mixing flow cons model should be used",
         &((*probdata)->flowConsMixing), FALSE, FALSE, NULL, NULL ) );

   SCIP_CALL( SCIPaddBoolParam(scip, "flowConsMixingNode",
         "whether the mixing flow cons model should be used",
         &((*probdata)->flowConsMixingNode), FALSE, FALSE, NULL, NULL ) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nodeMu",
         "whether node mumixing model should be used",
         &((*probdata)->nodeMu), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "mol",
         "whether only mol& without auxilary variable mu_mass should be used",
         &((*probdata)->mol), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "papay",
         "whether formula of papay is used for compressibility factor",
         &((*probdata)->papay), FALSE, FALSE, NULL, NULL ) );

   SCIP_CALL( SCIPaddBoolParam(scip, "additionalFacets",
         "whether additional facets of box constraint model should be used",
         &((*probdata)->additionalFacets), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxLowerBounds",
         "Should lower pressure bounds be relaxed?",
         &((*probdata)->relaxLowerBounds), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxUpperBounds",
         "Should upper pressure bounds be relaxed?",
         &((*probdata)->relaxUpperBounds), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxCVbounds",
         "whether pressure bounds of controlvalves should be relaxed by relaxCVvalue",
         &((*probdata)->relaxCVbounds), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxLBvalue",
         "relax lower pressure bounds by this much",
         &((*probdata)->relaxLBvalue), FALSE, 0.0, 0.0, 100.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxUBvalue",
         "relax upper pressure bounds by this much",
         &((*probdata)->relaxUBvalue), FALSE, 0.0, 0.0, 100.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxCVvalue",
         "relax pressure bounds of controlvalves by this much",
         &((*probdata)->relaxCVvalue), FALSE, 0.0, 0.0, 100.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "resistorCVInternal",
         "whether resistors of a controlvave are treated to be internal",
         &((*probdata)->resistorCVInternal), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "idealCSminIncrease",
         "minimal pressure increase by an ideal compressor",
         &((*probdata)->idealCSminIncrease), FALSE, 1.0, 0.0, 100.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "idealCSmaxIncrease",
         "maximal pressure increase by an ideal compressor",
         &((*probdata)->idealCSmaxIncrease), FALSE, 50.0, 0.0, 100.0, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "csusemaxincrease",
         "whether a maximal increase for compressors should be used",
         &((*probdata)->csusemaxincrease), FALSE, TRUE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "noObjective",
         "whether only the feasibility problem should be solved",
         &((*probdata)->noObjective), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "minLambda",
         "whether the overall resistance in the pipes should be minimized (only for mixing)",
         &((*probdata)->minLambda), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "powerLoss",
         "whether objective powerLoss should be used",
         &((*probdata)->powerLoss), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "minCompSum",
         "minimize sum of active compressor stations",
         &((*probdata)->minCompSum), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "minPressSum",
         "minimize the pressure of exits",
         &((*probdata)->minPressSum), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "minCompInc",
         "minimize the pressure increase by compressor stations",
         &((*probdata)->minCompInc), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "addComponentCuts",
         "whether component cuts should be added",
         &((*probdata)->addComponentCuts), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "preprocessPassiveComponents",
         "preprocess passive components and fix flow",
         &((*probdata)->preprocessPassiveComponents), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "reduceFlowConservation",
         "skip one flow conservation for one node, to make matrix full rank",
         &((*probdata)->reduceFlowConservation), FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "aggregateParallelFlowdirVars",
         "whether flow direction variables of parallel arcs should be aggregated",
         &(*probdata)->aggregateParallelFlowdirVars, FALSE, FALSE, NULL, NULL ) );

   SCIP_CALL( SCIPaddRealParam(scip, "scalescn",
         "factor to scale scenario flow values with",
         &(*probdata)->scalescn, FALSE, 1.0, 0.0, SCIP_REAL_MAX, NULL, NULL ) );

   return SCIP_OKAY;
}


/** create stationary gas transport instance */
SCIP_RETCODE GAScreateProb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   const char*           netfilename,        /**< name of network file to read */
   const char*           scnfilename,        /**< name of scenario file to read */
   const char*           csfilename          /**< name of compressor file to read */
   )
{
   char probname[SCIP_MAXSTRLEN];
   char netname[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( probdata != NULL );

   /* read network */
   SCIPinfoMessage(scip, NULL, "Network file name:\t%s\n", netfilename);
   SCIP_CALL( GASLIBreadFile(scip, netfilename, probdata, probdata->boxConstraintModel, probdata->algebraic) );

   /* read scenario */
   SCIPinfoMessage(scip, NULL, "\nScenario file name:\t%s\n", scnfilename);
   SCIP_CALL( SCNreadFile(scip, scnfilename, probdata) );

   /* read scenario */
   if ( csfilename != NULL )
   {
      SCIPinfoMessage(scip, NULL, "\nCompressor file name:\t%s\n", csfilename);
      SCIP_CALL( CSreadFile(scip, csfilename, probdata) );
   }

   /* compute spanning tree */
   if ( ! probdata->noCycles )
   {
      SCIP_CALL( computeSpanningTree(scip, probdata->network) );
      if ( probdata->network->spanningTree->foundCycles )
      {
         SCIP_CALL( findFundamentalCycles(scip, probdata->network) );
         if ( (probdata->network->numBasisCycles > 1) && probdata->allCycles )
         {
            SCIP_CALL( combineCycles(scip, probdata->network) );
         }
      }
   }

   /* allocate memory for filenames */
   SCIP_CALL( getProblemName(netfilename, netname, SCIP_MAXSTRLEN) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(probdata->netname), netname, strlen(netname)+1) );
   SCIP_CALL( getProblemName(scnfilename, probname, SCIP_MAXSTRLEN) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(probdata->scnname), probname, strlen(probname)+1) );

   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, probname, probdelorigGAS, NULL, NULL, NULL, NULL, probcopyGAS, probdata) );

   return SCIP_OKAY;
}

/** create problem instance with existing network and determine neighbors */
SCIP_RETCODE GAScreateProbNetwork(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   struct GAS_network*   network             /**< network */
   )
{
   /* create problem */
   SCIP_CALL( SCIPcreateProb(scip, probname, probdelorigGAS, NULL, NULL, NULL, NULL, probcopyGAS, probdata) );

   SCIP_CALL( determineNeighbors(scip, network) );

   return SCIP_OKAY;
}

/** write network in xml file */
SCIP_RETCODE GASwriteNetwork(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   const char*           filename,           /**< name of network file to write */
   SCIP_Bool             withframework       /**< use key "framework"? */
   )
{
   GAS_Network* network;
   struct tm* tm;
   FILE* file;
   time_t t;
   int i;

   file = fopen(filename, "w");
   if ( file == NULL )
   {
      SCIPerrorMessage("Could not open file <%s>.\n", filename);
      return SCIP_WRITEERROR;
   }

   t = time(NULL);
   if ( t == (time_t)(-1))
   {
      SCIPerrorMessage("Error with obtaining time.\n");
      return SCIP_READERROR;
   }
   tm = gmtime(&t);

   network = probdata->network;
   assert( network != NULL );

   if ( withframework )
   {
      fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
      fprintf(file, "<network xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
      fprintf(file, "         xmlns=\"http://gaslib.zib.de/Gas\"\n");
      fprintf(file, "         xsi:schemaLocation=\"http://gaslib.zib.de/Gas Gas.xsd\"\n");
      fprintf(file, "         xmlns:framework=\"http://gaslib.zib.de/Framework\">\n");
      fprintf(file, "  <framework:information>\n");
      fprintf(file, "    <framework:title>%s</framework:title>\n", probdata->netname);
      fprintf(file, "    <framework:type>gas</framework:type>\n");
      fprintf(file, "    <framework:date>%.4d-%.2d-%.2d</framework:date>\n", 1900 + tm->tm_year, tm->tm_mon + 1, tm->tm_mday);
      fprintf(file, "  </framework:information>\n");

      fprintf(file, "  <framework:nodes>\n");
   }
   else
   {
      fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
      fprintf(file, "<network xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n");
      fprintf(file, "         xmlns=\"http://gaslib.zib.de/Gas\"\n");
      fprintf(file, "         xsi:schemaLocation=\"http://gaslib.zib.de/Gas Gas.xsd\"\n");
      fprintf(file, "         xmlns:framework=\"http://gaslib.zib.de/Framework\">\n");
      fprintf(file, "  <information>\n");
      fprintf(file, "    <title>%s</title>\n", probdata->netname);
      fprintf(file, "    <type>gas</type>\n");
      fprintf(file, "    <date>%.4d-%.2d-%.2d</date>\n", 1900 + tm->tm_year, tm->tm_mon + 1, tm->tm_mday);
      fprintf(file, "  </information>\n");

      fprintf(file, "  <nodes>\n");
   }

   for (i = 0; i < network->numnodes; ++i)
   {
      GAS_Node* node;
      Node_type type;

      node = &(network->nodes_ptr[i]);

#if 1
      /* determine type based on scenario */
      if ( SCIPisPositive(scip, node->flow) )
         type = ENTRY;
      else if ( SCIPisNegative(scip, node->flow) )
         type = EXIT;
      else
         type = INNODE;
#else
      type = node->type;
#endif

      if ( type == ENTRY )
         fprintf(file, "    <source id=\"%s\" x=\"%g\" y=\"%g\">\n", node->id, node->x, node->y);
      else if ( type == EXIT )
         fprintf(file, "    <sink id=\"%s\" x=\"%g\" y=\"%g\">\n", node->id, node->x, node->y);
      else
         fprintf(file, "    <innode id=\"%s\" x=\"%g\" y=\"%g\">\n", node->id, node->x, node->y);

      fprintf(file, "      <height value=\"%g\" unit=\"m\"/>\n", node->height);
      fprintf(file, "      <pressureMin unit=\"bar\" value=\"%g\"/>\n", node->netPressureMin);
      fprintf(file, "      <pressureMax unit=\"bar\" value=\"%g\"/>\n", node->netPressureMax);

      if ( type == ENTRY || type == EXIT )
      {
         fprintf(file, "      <flowMin unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", node->flow);
         fprintf(file, "      <flowMax unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", node->flow);
      }

      /* write out data for first node */
      if ( type == ENTRY )
      {
         fprintf(file, "      <gasTemperature unit=\"Celsius\" value=\"%g\"/>\n", network->gasTemperature - 273.15);
         fprintf(file, "      <calorificValue unit=\"MJ_per_m_cube\" value=\"36.4543670654\"/>\n");
         fprintf(file, "      <normDensity unit=\"kg_per_m_cube\" value=\"%g\"/>\n", network->normDensity);
         fprintf(file, "      <coefficient-A-heatCapacity value=\"31.8251781464\"/>\n");
         fprintf(file, "      <coefficient-B-heatCapacity value=\"-0.00846800766885\"/>\n");
         fprintf(file, "      <coefficient-C-heatCapacity value=\"7.44647331885e-05\"/>\n");
         fprintf( file, "      <molarMass unit=\"kg_per_kmol\" value=\"%g\"/>\n", network->molarMass1 * 1000.0 );
         fprintf(file, "      <pseudocriticalPressure unit=\"bar\" value=\"%.11g\"/>\n", network->pseudocriticalPressure);
         fprintf(file, "      <pseudocriticalTemperature unit=\"K\" value=\"%.11g\"/>\n", network->pseudocriticalTemperature);
         /*
           <coefficient-A-heatCapacity value="31.8251781464"/>
           <coefficient-B-heatCapacity value="-0.00846800766885"/>
           <coefficient-C-heatCapacity value="7.44647331885e-05"/>
         */
      }

      if ( type == ENTRY )
         fprintf(file, "    </source>\n");
      else if ( type == EXIT )
         fprintf(file, "    </sink>\n");
      else
         fprintf(file, "    </innode>\n");
      fprintf(file, "\n");
   }

   if ( withframework )
   {
      fprintf(file, "  </framework:nodes>\n");
      fprintf(file, "  <framework:connections>\n");
   }
   else
   {
      fprintf(file, "  </nodes>\n");
      fprintf(file, "  <connections>\n");
   }

   for (i = 0; i < network->numarcs; ++i)
   {
      GAS_Arc* arc;

      arc = &(network->arcs_ptr[i]);

      if ( arc->type == PIPE )
      {
         GAS_Pipe* pipe;

         pipe = (GAS_Pipe*) arc->detailed_info;

         fprintf(file, "    <pipe id=\"%s\" from=\"%s\" to=\"%s\">\n", arc->id, arc->sourcenode->id, arc->targetnode->id);
         fprintf(file, "      <flowMin unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMin);
         fprintf(file, "      <flowMax unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMax);
         fprintf(file, "      <length unit=\"km\" value=\"%g\"/>\n", pipe->length /1000.0);
         fprintf(file, "      <diameter unit=\"mm\" value=\"%g\"/>\n", pipe->diameter * 1000.0);
         fprintf(file, "      <roughness unit=\"mm\" value=\"%g\"/>\n", pipe->roughness * 1000.0);
         fprintf(file, "      <pressureMax unit=\"bar\" value=\"%g\"/>\n", pipe->pressureMax);
         fprintf(file, "      <heatTransferCoefficient unit=\"W_per_m_square_per_K\" value=\"2\"/>\n");
         fprintf(file, "    </pipe>\n\n");
      }
      else if ( arc->type == SHORTPIPE )
      {
         fprintf(file, "    <shortPipe id=\"%s\" from=\"%s\" to=\"%s\">\n", arc->id, arc->sourcenode->id, arc->targetnode->id);
         fprintf(file, "      <flowMin unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMin);
         fprintf(file, "      <flowMax unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMax);
         fprintf(file, "    </shortPipe>\n\n");
      }
      else if ( arc->type == VALVE )
      {
         GAS_Valve* valve;

         valve = (GAS_Valve*) arc->detailed_info;

         fprintf(file, "    <valve id=\"%s\" from=\"%s\" to=\"%s\">\n", arc->id, arc->sourcenode->id, arc->targetnode->id);
         fprintf(file, "      <flowMin unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMin);
         fprintf(file, "      <flowMax unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMax);
         fprintf(file, "      <pressureDifferentialMax unit=\"bar\" value=\"%g\"/>\n", valve->pressureDifferentialMax);
         fprintf(file, "    </valve>\n\n");
      }
      else if ( arc->type == CONTROLVALVE )
      {
         GAS_Controlvalve* cv;

         cv = (GAS_Controlvalve*) arc->detailed_info;

         fprintf(file, "    <controlValve id=\"%s\" from=\"%s\" to=\"%s\" internalBypassRequired=\"%d\">\n",
            arc->id, arc->sourcenode->id, arc->targetnode->id, cv->internalBypassRequired);
         fprintf(file, "      <flowMin unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMin);
         fprintf(file, "      <flowMax unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMax);
         fprintf(file, "      <pressureDifferentialMin unit=\"bar\" value=\"%g\"/>\n", cv->pressureDifferentialMin);
         fprintf(file, "      <pressureDifferentialMax unit=\"bar\" value=\"%g\"/>\n", cv->pressureDifferentialMax);
         fprintf(file, "      <pressureInMin unit=\"bar\" value=\"%g\"/>\n", cv->pressureInMin);
         fprintf(file, "      <pressureOutMax unit=\"bar\" value=\"%g\"/>\n", cv->pressureOutMax);
         fprintf(file, "      <pressureLossIn unit=\"bar\" value=\"%g\"/>\n", cv->pressureLossIn);
         fprintf(file, "      <pressureLossOut unit=\"bar\" value=\"%g\"/>\n", cv->pressureLossOut);
         fprintf(file, "    </controlValve>\n\n");
      }
      else if ( arc->type == CS )
      {
         GAS_CS* cs;

         cs = (GAS_CS*) arc->detailed_info;

         fprintf(file, "    <compressorStation id=\"%s\" from=\"%s\" to=\"%s\">\n", arc->id, arc->sourcenode->id, arc->targetnode->id);
         fprintf(file, "      <flowMin unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMin);
         fprintf(file, "      <flowMax unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMax);
         fprintf(file, "      <pressureLossIn unit=\"bar\" value=\"%g\"/>\n", cs->pressureLossIn);
         fprintf(file, "      <pressureLossOut unit=\"bar\" value=\"%g\"/>\n", cs->pressureLossOut);
         fprintf(file, "      <pressureInMin unit=\"bar\" value=\"%g\"/>\n", cs->pressureInMin);
         fprintf(file, "      <pressureOutMax unit=\"bar\" value=\"%g\"/>\n", cs->pressureOutMax);
         fprintf(file, "    </compressorStation>\n\n");
      }
      else if ( arc->type == RESISTOR )
      {
         GAS_Resistor* resistor;

         resistor = (GAS_Resistor*) arc->detailed_info;

         fprintf(file, "    <resistor id=\"%s\" from=\"%s\" to=\"%s\">\n", arc->id, arc->sourcenode->id, arc->targetnode->id);
         fprintf(file, "      <flowMin unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMin);
         fprintf(file, "      <flowMax unit=\"1000m_cube_per_hour\" value=\"%g\"/>\n", arc->flowMax);
         fprintf(file, "      <dragFactor value=\"%g\"/>\n", resistor->dragFactor);
         fprintf(file, "      <diameter unit=\"mm\" value=\"%g\"/>\n", resistor->diameter * 1000.0);
         fprintf(file, "    </resistor>\n\n");
      }
      else
      {
         SCIPerrorMessage("Cannot yet write out type %d.\n", arc->type);
      }
   }

   if ( withframework )
   {
      fprintf(file, "  </framework:connections>\n");
      fprintf(file, "</network>\n");
   }
   else
   {
      fprintf(file, "  </connections>\n");
      fprintf(file, "</network>\n");
   }

   fclose(file);

   return SCIP_OKAY;
}
