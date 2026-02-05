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

/**@file   generateModel_gas.c
 * @brief  generate constraints for describing stationary gas transport networks
 * @author Imke Joormann
 * @author Alexandra Stern
 * @author Oliver Habeck
 * @author Dennis Gabriel
 * @author Pascal Boerner
 */

#include <scip/scipdefplugins.h>
#include "generateModel_gas.h"
#include "characteristics_gas.h"
#include "scip/scip.h"
#include "cons_pipe_ode.h"
#include "probdata_gas.h"
#include <string.h>

#if SCIP_VERSION >= 800
/* direct access to the expression struct is a work-around and should be replaced */
#include "scip/struct_expr.h"
#endif

/* #define SCIP_DEBUG */

/*------------------local functions------------------------------------- */

#if SCIP_VERSION >= 800
/** frees an expression */
static
SCIP_RETCODE freeExpression(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_EXPR**           expr                /**< pointer to free the expression */
   )
{
   assert( scip != NULL );
   assert( expr != NULL );
   assert( (*expr)->quaddata == NULL );
   assert( (*expr)->ownerdata == NULL );

   if ( (*expr)->exprdata != NULL )
   {
      assert( (*expr)->exprhdlr->freedata != NULL );
      SCIP_CALL( (*expr)->exprhdlr->freedata(scip, *expr) );
   }

   /* free children array, if any */
   SCIPfreeBlockMemoryArray(scip, &(*expr)->children, (*expr)->childrensize);

   SCIPfreeBlockMemory(scip, expr);
   assert( *expr == NULL );

   return SCIP_OKAY;
}

/** recurively frees an expression */
static
SCIP_RETCODE GASexprFreeRecurisve(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_EXPR**           rootexpr            /**< pointer to expression */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPR* expr;

   assert( scip != NULL );
   assert( rootexpr != NULL );
   assert( *rootexpr != NULL );

   /* create iterator through subtree */
   SCIP_CALL( SCIPcreateExpriter(scip, &it) );

   /* use DFS and do not revisit nodes */
   SCIP_CALL( SCIPexpriterInit(it, *rootexpr, SCIP_EXPRITER_DFS, FALSE) );

   /* we need to free expressions when we leave them */
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);
   for (expr = SCIPexpriterGetCurrent(it); ! SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it)) /*lint !e2840*//*lint !e440*/
   {
      assert( SCIPexpriterGetStageDFS(it) == SCIP_EXPRITER_LEAVEEXPR );

      SCIP_CALL( freeExpression(scip, &expr) );
   }

   SCIPfreeExpriter(&it);

   return SCIP_OKAY;
}
#endif /* SCIP_VERSION >= 800 */

/** create continuous and binary flowvariables */
static
SCIP_RETCODE createFlowVariables(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Arc* arc;
   GAS_Node* node;
   char varname[SCIP_MAXSTRLEN];
   SCIP_Real lb;
   SCIP_Real ub;
   int k;
   int j;
   int flowdirposition = 0;

   assert( scip != NULL );

   probdata = SCIPgetProbData(scip);

   /* allocate memory */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->flowvars, probdata->network->numarcs) );

   if ( probdata->mixing || probdata->approx )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->Lambda, probdata->network->numarcs) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->nodeMixingRatio, probdata->network->numnodes) );

      if ( ! probdata->nodeMu )
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->mixingRatio, probdata->network->numarcs) );
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->mixingRatioMass, probdata->network->numarcs) );
      }

      if ( probdata->nodeMu )
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->nodeMixingRatioMol, probdata->network->numnodes) );
      }
   }

   if ( probdata->linearMixing )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->Lambda, probdata->network->numarcs) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->nodeMixingRatio, probdata->network->numnodes) );

      if ( ! probdata->nodeMu && ! probdata->flowConsMixingNode )
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->mixingRatioMass, probdata->network->numarcs) );
      }
   }

   /* compressors and controlvalves have a fixed flow direction, nonlinear resistors already have one in their model */
   if ( probdata->network->numflowdirvars > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->positiveFlowBinvars, probdata->network->numflowdirvars) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->negativeFlowBinvars, probdata->network->numflowdirvars) );
   }

   /* loop through all arcs and generate the variables */
   for (k = 0; k < probdata->network->numarcs; ++k)
   {
      arc = &probdata->network->arcs_ptr[k];
      assert( arc->type != UNKNOWNARCTYPE );

      /*
       *  first create the continuous flow variable
       */

      /* name of the flowvariable */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "q#%s#", arc->id);

      /* set the right flow bounds */
      ub = arc->flowMax;
      lb = arc->flowMin;

      /* if the arc is a compressor or a controlvalve with a bypass we have to be able to turn it off/close it
       * therefore we set the minimal flow to 0.0
       * if there has to be flow (positive) over the compressor station, we have either minflow on
       * the compressor or the bypass: 0.0 < flowMin \leq q_cs + q_va
       * generate compressor constraints
       */
      if ( arc->type == CS )
      {
         if ( ((GAS_CS*)arc->detailed_info)->internalBypassRequired == 1 )
         {
            lb = 0.0;
         }
      }
      else if ( arc->type == CONTROLVALVE )
      {
         if ( ((GAS_Controlvalve*)arc->detailed_info)->internalBypassRequired == 1 )
         {
            lb = 0.0;
         }
      }

      /* now create variable and add it */
      SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->flowvars[k]), varname, lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS) );

      /* do not multi-aggregate flow variables for pipes */
      if ( arc->type == PIPE )
      {
	 assert( probdata->flowvars[k] != NULL );
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, probdata->flowvars[k]) );
      }

      SCIP_CALL( SCIPaddVar(scip, probdata->flowvars[k]) );

      /* link flowvariable in arc */
      arc->flowvar = probdata->flowvars[k];

      /*
       * initialize the binary direction variables
       */

      /* only create flow binvars if option noFlowBinvars is false */
      if ( probdata->noFlowBinvars )
         continue;

      /* positive flow binvar */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "posFlowBinvar#%s#", arc->id);

      /* now create variable and add it */
      SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->positiveFlowBinvars[flowdirposition]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, probdata->positiveFlowBinvars[flowdirposition]) );

      /* build link such that the flowvariable can be adressed via Gas_Arc */
      arc->positiveFlowBinvar = probdata->positiveFlowBinvars[flowdirposition];

      /* negative flow binvar */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "negFlowBinvar#%s#", arc->id);
      /* now create variable and add it */
      SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->negativeFlowBinvars[flowdirposition]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
      SCIP_CALL( SCIPaddVar(scip, probdata->negativeFlowBinvars[flowdirposition]) );

      /* multi aggregation of these variables leads to some problem .... */
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, probdata->negativeFlowBinvars[flowdirposition]) );
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, probdata->positiveFlowBinvars[flowdirposition]) );

      /* build link such that the flowvariable can be adressed via Gas_Arc */
      arc->negativeFlowBinvar = probdata->negativeFlowBinvars[flowdirposition];
      ++flowdirposition;
      assert( flowdirposition <= probdata->network->numflowdirvars );

      /* Compressor stations and control valves have fixed flow direction.
       *
       * Here only the negative flow binvar is set to zero. The positive flow bin var is coupled to the compressor bin
       * var later in: generateIdealCSConstraints and generateBoxConstraintModel. 
       * For control valves the coupling happens in
       * generateControlValveConstraintsExternalResistors() and generateControlValveConstraintsInternalResistors(). */
      if ( arc->type == CS || arc->type == CONTROLVALVE )
      {
         SCIP_CALL( SCIPchgVarLb(scip, arc->negativeFlowBinvar, 0.0) );
         SCIP_CALL( SCIPchgVarUb(scip, arc->negativeFlowBinvar, 0.0) );
      }
   }

   if ( probdata->mixing || probdata->approx )
   {
      /* loop through all arcs and generate the variables */
      for (k = 0; k < probdata->network->numarcs; ++k)
      {
         arc = &probdata->network->arcs_ptr[k];
         assert( arc->type != UNKNOWNARCTYPE );

         if ( ! probdata->nodeMu )
         {
            /* name of the ratio variable */
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "mu#%s#", arc->id);

            /* now create variable and add it */
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->mixingRatio[k]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

            SCIP_CALL( SCIPaddVar(scip, probdata->mixingRatio[k]) );

            /* link ratio variable in arc */
            arc->mixingRatio = probdata->mixingRatio[k];

            /* name of the mass ratio variable */
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "mu_mass#%s#", arc->id);

            /* now create variable and add it */
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->mixingRatioMass[k]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

            SCIP_CALL( SCIPaddVar(scip, probdata->mixingRatioMass[k]) );

            /* link ratio variable in arc */
            arc->mixingRatioMass = probdata->mixingRatioMass[k];
         }

         /* name of the variable */
         (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "LAMBDA#%s#", arc->id);

         if ( probdata->minLambda )
         {
            /* now create variable and add it */
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->Lambda[k]), varname, 0, 100, 1.0, SCIP_VARTYPE_CONTINUOUS) );
         }
         else
         {
            /* now create variable and add it */
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->Lambda[k]), varname, 0, 100, 0.0, SCIP_VARTYPE_CONTINUOUS) );
         }

         SCIP_CALL( SCIPaddVar(scip, probdata->Lambda[k]) );

         /* link mixed speed of sound variable in arc */
         arc->Lambda = probdata->Lambda[k];
      }

      /* loop through all nodes and generate the variables */
      for (k = 0; k < probdata->network->numnodes; ++k)
      {
         node = &probdata->network->nodes_ptr[k];

         /* name of the node mixing ratio variable */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "mu_mass_node#%s#", node->id);

         /* now create variable and add it */
         SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->nodeMixingRatio[k]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

         SCIP_CALL( SCIPaddVar(scip, probdata->nodeMixingRatio[k]) );

         /* link ratio variable in node */
         node->nodeMixingRatio = probdata->nodeMixingRatio[k];

         if ( probdata->nodeMu )
         {
            /* name of the node mixing ratio variable */
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "mu_mol_node#%s#", node->id);

            /* now create variable and add it */
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->nodeMixingRatioMol[k]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

            SCIP_CALL( SCIPaddVar(scip, probdata->nodeMixingRatioMol[k]) );

            /* link ratio variable in node */
            node->nodeMixingRatioMol = probdata->nodeMixingRatioMol[k];
         }
      }
   }

   if ( probdata->linearMixing )
   {
      /* loop through all arcs and generate the variables */
      for (k = 0; k < probdata->network->numarcs; ++k)
      {
         arc = &probdata->network->arcs_ptr[k];
         assert( arc->type != UNKNOWNARCTYPE );

         if ( ! probdata->nodeMu && ! probdata->flowConsMixingNode )
         {
            /* name of the mass ratio variable */
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "mu_mass#%s#", arc->id);

            /* now create variable and add it */
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->mixingRatioMass[k]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

            SCIP_CALL( SCIPaddVar(scip, probdata->mixingRatioMass[k]) );

            /* link ratio variable in arc */
            arc->mixingRatioMass = probdata->mixingRatioMass[k];
         }

         /* name of the variable */
         (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "LAMBDA#%s#", arc->id);

         if ( probdata->minLambda )
         {
            /* now create variable and add it */
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->Lambda[k]), varname, 0, 100, 1.0, SCIP_VARTYPE_CONTINUOUS) );
         }
         else
         {
            /* now create variable and add it */
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->Lambda[k]), varname, 0, 100, 0.0, SCIP_VARTYPE_CONTINUOUS) );
         }

         SCIP_CALL( SCIPaddVar(scip, probdata->Lambda[k]) );

         /* link mixed speed of sound variable in arc */
         arc->Lambda = probdata->Lambda[k];
      }

      /* loop through all nodes and generate the variables */
      for (k = 0; k < probdata->network->numnodes; ++k)
      {
         node = &probdata->network->nodes_ptr[k];

         /* name of the ndoe mixing ratio variable */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "mu_mass_node#%s#", node->id);

         /* now create variable and add it */
         SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->nodeMixingRatio[k]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

         SCIP_CALL( SCIPaddVar(scip, probdata->nodeMixingRatio[k]) );

         /* link ratio variable in node */
         node->nodeMixingRatio = probdata->nodeMixingRatio[k];
      }
   }

   if ( ! probdata->noFlowBinvars )
   {
      if ( flowdirposition != probdata->network->numflowdirvars )
      {
         SCIPerrorMessage("Something went wrong while creating the binary flow direction variables.\n");
         return SCIP_ERROR;
      }
   }

   if ( probdata->maxFlow )
   {
      assert( probdata->netflowvars == NULL );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->netflowvars, probdata->network->numnodes) );

      /* loop through all nodes and set netflow variable */
      for (j = 0; j < probdata->network->numnodes; j++)
      {
         GAS_Node* Node = &(probdata->network->nodes_ptr[j]);
         assert( Node->type != UNKNOWNNODETYPE );

         if ( Node->type == ENTRY ) /* true if node is entry */
         {
            /* name of the netflowvariable */
            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "netflow#%s#", Node->id);

            lb = Node->scnFlowMin;
            ub = Node->scnFlowMax;

            /* create variable */
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->netflowvars[j]), varname, lb, ub, -1.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, probdata->netflowvars[j]) );
         }
      }
   }

   return SCIP_OKAY;
}

/** create pressure, pressure squared and pressure slack variables */
static
SCIP_RETCODE createPressureVariables(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Node* N;
   char varname[SCIP_MAXSTRLEN];
   SCIP_Real lb, ub;
   SCIP_Real maxSlack;
   SCIP_Real obj;
   int k;
   int nPIvarsCreated = 0;

   probdata = SCIPgetProbData(scip);

   /* set the right upper bound for slack variables */
   maxSlack = MAX( MAX( probdata->relaxUBvalue, probdata->relaxLBvalue ), probdata->relaxCVvalue );

   /*
    * allocate memory
    */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->pressurevars, probdata->network->numnodes) );
   if ( probdata->network->numpivars > 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->PIvars, probdata->network->numpivars) );
   }

   if ( probdata->minSlackPerBound )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->pSlackVars, probdata->network->numnodes) );
   }

   /*
    * set the objective coefficient of pressure variables ( only the case powerloss depends on the specific node )
    */

   /* default objective is to maximize the pressure */
   obj = -1.0;

   if ( probdata->minPressSum )
   {
      /* the objective function minimizes the sum of all pressure values. */
      obj = 1.0;
   }
   else if ( probdata->noObjective || probdata->minCompSum || probdata->minCompInc || probdata->minSlack || probdata->minSlackPerBound || probdata->maxFlow || probdata->minLambda)
   {
      /* in all this cases the pressure variables do not appear in the objective function */
      obj = 0.0;
   }

   /* loop through all nodes and create corresponding variables */
   for (k = 0; k < probdata->network->numnodes; ++k)
   {
      N = &(probdata->network->nodes_ptr[k]);

      /* determine maximal pressure at node */
      N->pressureMax = MIN3(N->netPressureMax, N->scnPressureMax, N->arcPressureMax);

      /* add slack */
      if ( probdata->relaxUpperBounds )
      {
         N->pressureMax += probdata->relaxUBvalue;
      }

      /* determine minimal pressure at node (there is no minimal pressure for arcs) */
      N->pressureMin = MAX(N->netPressureMin, N->scnPressureMin);

      /* add slack */
      if ( probdata->relaxLowerBounds )
      {
         N->pressureMin = MAX(N->pressureMin - probdata->relaxLBvalue, 1.01325);
      }

      lb = N->pressureMin;
      ub = N->pressureMax;

      /* set the objective coefficient in the case powerloss */
      if ( probdata->powerLoss )
      {
         /* objective is to minimize the power loss \sum_v q_v * p_v */
         obj = N->flow;
      }

      /* set variable name */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "p#%s#", N->id);

      SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->pressurevars[k]), varname, lb, ub, obj, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, probdata->pressurevars[k]) );

      SCIP_CALL( SCIPaddVar(scip, probdata->pressurevars[k]) );

      /* Create a link such that variable can be adressed through the corresponding node*/
      N->pressurevar = probdata->pressurevars[k];

      /* create pressure squared variable if need. e.g., when using the algebraic model */
      if ( N->needPIvar )
      {
         /* set variable name */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "PI#%s#", N->id);

         SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->PIvars[nPIvarsCreated]), varname, lb * lb, ub * ub, 0.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, probdata->PIvars[nPIvarsCreated]) );

         /* link PI variable to corresponding node*/
         N->PIvar = probdata->PIvars[nPIvarsCreated];
         ++nPIvarsCreated;
      }

      /* create slack variable for minimizing the relaxation */
      if ( probdata->minSlackPerBound )
      {
         /* set variable name */
         (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "pSlackVar#%s#", N->id);

         SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->pSlackVars[k]), varname, 0.0, maxSlack, 1.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPmarkDoNotMultaggrVar(scip, probdata->pSlackVars[k]) );

         SCIP_CALL( SCIPaddVar(scip, probdata->pSlackVars[k]) );

         /* Create a link such that variable can be adressed through the corresponding node*/
         N->pSlackVar = probdata->pSlackVars[k];
      }
   }

   if ( nPIvarsCreated != probdata->network->numpivars )
   {
      SCIPerrorMessage("Something went wrong when creating the pressure squared variables.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** create variables for algebraic model */
static
SCIP_RETCODE createVariablesAlgebraicModel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Arc* arc;
   GAS_Pipe* pipe;
   SCIP_Real  lb;
   SCIP_Real  ub;
   GAS_Node* source;
   GAS_Node* target;
   int varsCreated = 0;
   int k;
   char varname[SCIP_MAXSTRLEN];

   probdata = SCIPgetProbData(scip);

   /* allocate memory */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->PIdiffvars, probdata->network->numpipes) );

   for (k = 0; k < probdata->network->numarcs; ++k)
   {
      arc = &(probdata->network->arcs_ptr[k]);
      assert( arc->type != UNKNOWNARCTYPE );
      if ( arc->type != PIPE )
         continue;

      pipe = arc->detailed_info;
      source = arc->sourcenode;
      target = arc->targetnode;

      /* calculate upper and lower bounds for the squared pressure difference */
      lb = (source->pressureMin) * (source->pressureMin) - (target->pressureMax) * (target->pressureMax) ;
      ub = (source->pressureMax) * (source->pressureMax) - (target->pressureMin) * (target->pressureMin);

      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "PI_diff#%s#", arc->id);
      SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->PIdiffvars[varsCreated]), varname, lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(scip, probdata->PIdiffvars[varsCreated]) );

      pipe->PIdiffVar = probdata->PIdiffvars[varsCreated];
      ++varsCreated;
   }

   if ( varsCreated != probdata->network->numpipes )
   {
      SCIPerrorMessage("Something went wrong when creating variables for algebraic model. Did not create as many variables as there are pipes.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** create variables for valves */
static
SCIP_RETCODE createVariablesForValves(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Arc* arc;
   GAS_Valve* valve;
   int varsCreated = 0;
   int k;
   char varname[SCIP_MAXSTRLEN];

   assert( scip != NULL );

   probdata = SCIPgetProbData(scip);

   /* allocate memory */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->VALVE_binvars, probdata->network->numvalves) );

   /* loop through all arcs */
   for (k = 0; k < (probdata->network->numarcs); k++)
   {
      arc = &(probdata->network->arcs_ptr[k]);
      if ( arc->type != VALVE )
         continue;

      valve = arc->detailed_info;

      /* name */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "VALVE_binvar#%s#", arc->id);

      /* creating valve_binvar*/
      SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->VALVE_binvars[varsCreated]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );

      SCIP_CALL( SCIPaddVar(scip, probdata->VALVE_binvars[varsCreated]) );
      valve->valve_binvar = probdata->VALVE_binvars[varsCreated];
      ++varsCreated;
   }

   if ( varsCreated != probdata->network->numvalves )
   {
      SCIPerrorMessage("Something went wrong when creating variables for valves. Did not create as many variables as there are valves.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** create variables for valves */
static
SCIP_RETCODE createVariablesForControlvalves(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Arc* arc;
   GAS_Controlvalve* cv;
   int varsCreated = 0;
   int k;
   char varname[SCIP_MAXSTRLEN];

   assert( scip != NULL );

   probdata = SCIPgetProbData(scip);

   /* allocate memory */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->CV_binvars, probdata->network->numcontrolvalves) );

   /* loop through all arcs */
   for (k = 0; k < probdata->network->numarcs; k++)
   {
      arc = &(probdata->network->arcs_ptr[k]);
      if ( arc->type != CONTROLVALVE )
         continue;

      cv = arc->detailed_info;

      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "CV_binvar#%s#", arc->id);

      SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->CV_binvars[varsCreated]),varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );

      SCIP_CALL( SCIPaddVar(scip, probdata->CV_binvars[varsCreated]) );

      cv->binvar = probdata->CV_binvars[varsCreated];

      ++varsCreated;
   }

   if ( varsCreated != probdata->network->numcontrolvalves )
   {
      SCIPerrorMessage("Something went wrong when creating variables for controlvalves. Did not create as many variables as there are controlvalves.\n");
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** create variables needed for nonlinear resistor constraint */
static
SCIP_RETCODE createVariablesforNonLinResistor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   arcposition,        /**< position of arc */
   int                   resistorposition    /**< position of resitor */
   )
{
   GAS_Arc* arc;
   GAS_Node* Sourcenode;
   GAS_Node* Targetnode;
   SCIP_Real pressureDiffMin;
   SCIP_Real pressureDiffMax;
   SCIP_Real absDeltaMin;
   SCIP_Real absDeltaMax;
   SCIP_Real absFlowMin;
   SCIP_Real absFlowMax;
   char name[SCIP_MAXSTRLEN];

   Sourcenode = probdata->network->arcs_ptr[arcposition].sourcenode;
   Targetnode = probdata->network->arcs_ptr[arcposition].targetnode;
   arc = &probdata->network->arcs_ptr[arcposition];

   /* calculate upper and lower bound for variable Delta, absDelta and absFlow */
   pressureDiffMin = (Sourcenode->pressureMin)  - (Targetnode->pressureMax);
   pressureDiffMax = (Sourcenode->pressureMax)  - (Targetnode->pressureMin);
   absDeltaMin = ABS(pressureDiffMin) * pressureDiffMin;
   absDeltaMax = ABS(pressureDiffMax) * pressureDiffMax;
   absFlowMin = ABS(arc->flowMin) * arc->flowMin;
   absFlowMax = ABS(arc->flowMax) * arc->flowMax;

   /* for each resistor we need a binary variable resistor_binvars and support
      variables Delta, AbsDeltaVars, AbsFlowVars */

   /* for each nonlinear resistor generate the continuous variable Delta=p_u-p_v */

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "NLR_DeltaVar#%s#", arc->id);

   SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->NLR_DeltaVars[resistorposition]), name,
         pressureDiffMin, pressureDiffMax, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, probdata->NLR_DeltaVars[resistorposition]) );

   /* continuous variable AbsDeltaVar = |Delta|*Delta */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "NLR_AbsDeltaVar#%s#", arc->id);

   SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->NLR_AbsDeltaVars[resistorposition]),
         name, absDeltaMin, absDeltaMax, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, probdata->NLR_AbsDeltaVars[resistorposition]) );

   /* continuous variable AbsFlowVar = |flowvar|*flowvar */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "NLR_AbsFlowVar#%s#", arc->id);

   SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->NLR_AbsFlowVars[resistorposition]),
         name, absFlowMin, absFlowMax, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, probdata->NLR_AbsFlowVars[resistorposition]) );

   return SCIP_OKAY;
}

/** create three variables needed for linear resistor constraint */
static
SCIP_RETCODE createVariablesforLinResistor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   arcposition,        /**< position of arc */
   int                   resistorposition    /**< position of resistor */
   )
{
   char varname[SCIP_MAXSTRLEN];

   /* Creates binary variable LR_posFlowDir*/
   (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "LR_posFlowDir#%s#", probdata->network->arcs_ptr[arcposition].id);
   SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->LR_posFlowDir[resistorposition]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, probdata->LR_posFlowDir[resistorposition]) );

   /* Creates binary variable LR_negFlowDir */
   (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "LR_negFlowDir#%s#", probdata->network->arcs_ptr[arcposition].id);
   SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->LR_negFlowDir[resistorposition]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
   SCIP_CALL( SCIPaddVar(scip, probdata->LR_negFlowDir[resistorposition]) );

   /* Creates continuous variable LR_smoothingFlow */
   (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "LR_smoothingFlow#%s#", probdata->network->arcs_ptr[arcposition].id);
   SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->LR_smoothingFlow[resistorposition]), varname, -1.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, probdata->LR_smoothingFlow[resistorposition]) );

   return SCIP_OKAY;
}

/** create variables for resistors */
static
SCIP_RETCODE createResistorVars(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Arc* arc;
   GAS_Resistor* res;
   int k;

   assert( scip != NULL );

   probdata = SCIPgetProbData(scip);

   /* allocate memory for nonlinear resistors */
   if ( probdata->network->numnonlinresistor != 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->NLR_DeltaVars, probdata->network->numnonlinresistor) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->NLR_AbsDeltaVars, probdata->network->numnonlinresistor) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->NLR_AbsFlowVars, probdata->network->numnonlinresistor) );
   }

   /* allocate memory for linear resistors */
   if ( probdata->network->numlinresistor != 0 )
   {
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->LR_posFlowDir, probdata->network->numlinresistor) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->LR_smoothingFlow, probdata->network->numlinresistor) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->LR_negFlowDir, probdata->network->numlinresistor) );
   }

   for (k = 0; k < probdata->network->numarcs; ++k)
   {
      arc = &(probdata->network->arcs_ptr[k]);

      if ( arc->type != RESISTOR )
         continue;

      res = arc->detailed_info;

      if ( res->type == NONLINEAR )
      {
         SCIP_CALL( createVariablesforNonLinResistor(scip, probdata, k, res->NLR_position) );
      }
      else if ( res->type == LINEAR)
      {
         SCIP_CALL( createVariablesforLinResistor(scip, probdata, k, res->LR_position) );
      }
      else
      {
         SCIPerrorMessage(" Unknown resistor type.\n");
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}


/** create variables for compressor stations (except flowvariable) */
static
SCIP_RETCODE createVariablesForCompressors(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Arc*       arc;
   GAS_CS*        cs;
   GAS_CSBoxCons* boxcons;
   SCIP_Real      obj;
   SCIP_Real      lb;
   SCIP_Real      ub;
   int            k,l;
   int            binvarsCreated = 0;
   int            bcmvarsCreated = 0;
   int            resvarsCreated = 0;
   char           varname[SCIP_MAXSTRLEN];

   assert( scip != NULL );

   probdata = SCIPgetProbData(scip);

   /* if the sum of active compressors shall be minimized */
   if ( probdata->minCompSum )
      obj = 1.0;
   else
      obj = 0.0;

   /* allocate memory */
   if ( ! probdata->boxConstraintModel )
   {
      /* so far only boxconstraintmodel and idealized model are implemented, i.e.,
       * this is the case of idealized compressor model */
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->CS_binvars, probdata->network->numcompressor) );
   }
   else
   {
      /* we need a binary variable for every station and every configuration */
      k = probdata->network->numcompressor + probdata->network->numcsconfigurations;
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->CS_binvars, k) );

      /* Memory for nonlinear resistor variables in compressor stations */
      if ( probdata->network->numcsvarsfornlrs != 0 )
      {
         SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->CS_resistorvars, probdata->network->numcsvarsfornlrs) );
      }

      /* we need a flow, a pressure_in and a pressure_out variable for each config */
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->BCM_flowCS, probdata->network->numcsconfigurations) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->BCM_pressureIn, probdata->network->numcsconfigurations) );
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &probdata->BCM_pressureOut, probdata->network->numcsconfigurations) );
   }

   for (k = 0; k < (probdata->network->numarcs); k++)
   {
      arc = &(probdata->network->arcs_ptr[k]);

      if ( arc->type != CS )
         continue;

      cs = arc->detailed_info;

      /* create binary variable for compressor station */
      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "CS_binvar#%s#", arc->id);
      SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->CS_binvars[binvarsCreated]), varname, 0.0, 1.0, obj, SCIP_VARTYPE_BINARY) );

      SCIP_CALL( SCIPaddVar(scip, probdata->CS_binvars[binvarsCreated]) );
      cs->compressor_binvar = probdata->CS_binvars[binvarsCreated];

      ++binvarsCreated;

      /*
       *  create variables for the box constraint model
       */
      if ( probdata->boxConstraintModel )
      {
         SCIP_Real objIn = 0.0;
         SCIP_Real objOut = 0.0;
         if (probdata->minCompInc)
         {
            objIn = 0.0;
            objOut = 1.0;
         }

         /*
          *  create variables for configurations
          */
         for (l = 0; l < cs->numconfigurations; ++l)
         {
            boxcons = &(cs->boxcons[l]);

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "CS_ConfigBinvar#%s_%s#", arc->id, boxcons->id);

            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->CS_binvars[binvarsCreated]), varname, 0.0, 1.0, 0.0, SCIP_VARTYPE_BINARY) );
            SCIP_CALL( SCIPaddVar(scip, probdata->CS_binvars[binvarsCreated]) );
            boxcons->config_binvar = probdata->CS_binvars[binvarsCreated];
            ++binvarsCreated;

            /* create variable BCMflowCS  */
            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "BCM_flowCS#%s_%s#", arc->id, boxcons->id);
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->BCM_flowCS[bcmvarsCreated]), varname, boxcons->massFlowMin, boxcons->massFlowMax,
                  0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, probdata->BCM_flowCS[bcmvarsCreated]) );
            boxcons->BCM_flow = probdata->BCM_flowCS[bcmvarsCreated];

            /* create variable BCMpressureIn */
            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "BCM_Pin#%s_%s#", arc->id, boxcons->id);
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->BCM_pressureIn[bcmvarsCreated]), varname, boxcons->pressureInMin, boxcons->pressureInMax,
                  objIn, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, probdata->BCM_pressureIn[bcmvarsCreated]) );
            boxcons->BCM_pin = probdata->BCM_pressureIn[bcmvarsCreated];

            /* create variable BCMpressureOut */
            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "BCM_Pout#%s_%s#", arc->id, boxcons->id);
            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->BCM_pressureOut[bcmvarsCreated]), varname, boxcons->pressureOutMin, boxcons->pressureOutMax,
                  objOut, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, probdata->BCM_pressureOut[bcmvarsCreated]) );
            boxcons->BCM_pout = probdata->BCM_pressureOut[bcmvarsCreated];

            ++bcmvarsCreated;
         }

         /*
          *  create variables for nonlinear resistors
          */
         /* variables for in resistor */
         if ( ! SCIPisZero(scip, cs->dragFactorIn) )
         {
            /* create pressure variable NLRin_pressure */
            lb = arc->sourcenode->pressureMin;
            ub = arc->sourcenode->pressureMax;

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "CS_ResInP#%s#", arc->id);

            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->CS_resistorvars[resvarsCreated]), varname, lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS) );
            /* add variable */
            SCIP_CALL( SCIPaddVar(scip, probdata->CS_resistorvars[resvarsCreated]) );
            /* Create a link */
            cs->NLRin_pressure = probdata->CS_resistorvars[resvarsCreated];

            ++resvarsCreated;

            /* create pressure difference variable NLRin_pressureDiff */
            lb = 0.0;
            ub = arc->sourcenode->pressureMax - arc->sourcenode->pressureMin;

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "CS_ResInPDiff#%s#", arc->id);

            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->CS_resistorvars[resvarsCreated]), varname, lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS) );
            /* add variable */
            SCIP_CALL( SCIPaddVar(scip, probdata->CS_resistorvars[resvarsCreated]) );
            /* Create a link */
            cs->NLRin_pressureDiff = probdata->CS_resistorvars[resvarsCreated];

            ++resvarsCreated;
         }

         /* variables for out resistor */
         if ( ! SCIPisZero(scip, cs->dragFactorOut) )
         {
            /* create pressure variable NLRout_pressure */
            lb = arc->targetnode->pressureMin;
            ub = arc->targetnode->pressureMax;

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "CS_ResOutP#%s#", arc->id);

            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->CS_resistorvars[resvarsCreated]), varname, lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS) );
            /* add variable */
            SCIP_CALL( SCIPaddVar(scip, probdata->CS_resistorvars[resvarsCreated]) );
            /* Create a link */
            cs->NLRout_pressure = probdata->CS_resistorvars[resvarsCreated];

            ++resvarsCreated;

            /* create pressure difference variable NLRin_pressureDiff */
            lb = 0.0;
            ub = arc->targetnode->pressureMax - arc->targetnode->pressureMin;

            (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "CS_ResOutPDiff#%s#", arc->id);

            SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->CS_resistorvars[resvarsCreated]), varname, lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS) );
            /* add variable */
            SCIP_CALL( SCIPaddVar(scip, probdata->CS_resistorvars[resvarsCreated]) );
            /* Create a link */
            cs->NLRout_pressureDiff = probdata->CS_resistorvars[resvarsCreated];

            ++resvarsCreated;
         }

         /* end if boxconstraintmodel */
      }

      /* end loop through all constraints */
   }

   /* consistency check */
   if ( probdata->boxConstraintModel )
   {
      if ( binvarsCreated != (probdata->network->numcompressor +  probdata->network->numcsconfigurations) )
      {
         SCIPerrorMessage("Did not create as many binary variables for compressor stations as there should be!\n");
         return SCIP_ERROR;
      }
      if ( bcmvarsCreated != probdata->network->numcsconfigurations )
      {
         SCIPerrorMessage("Did not create as many binary variables for box constraint model as there should be!\n");
         return SCIP_ERROR;
      }
      if ( resvarsCreated != probdata->network->numcsvarsfornlrs )
      {
         SCIPerrorMessage("Did not create as many variables for nonlinear resistors in compressor stations as there should be!\n");
         return SCIP_ERROR;
      }
   }
   else
   {
      if ( binvarsCreated != probdata->network->numcompressor )
      {
         SCIPerrorMessage("Did not create as many binary variables for compressor stations as there should be!\n");
         return SCIP_ERROR;
      }
   }

   return SCIP_OKAY;
}

/** generate variables */
static
SCIP_RETCODE generateVariables(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   char varname[SCIP_MAXSTRLEN];
   SCIP_Real ub = 0.0;

   probdata = SCIPgetProbData(scip);
   assert ( probdata != NULL );
   assert ( probdata->network != NULL );

   /* create variable for objective function */
   if ( probdata->minSlack )
   {
      ub = MAX( MAX(probdata->relaxLBvalue, probdata->relaxUBvalue), probdata->relaxCVvalue);

      (void) SCIPsnprintf(varname, SCIP_MAXSTRLEN, "objective");
      SCIP_CALL( SCIPcreateVarBasic(scip, &(probdata->objective), varname, 0.0, ub, 1.0, SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(scip, probdata->objective) );
   }

   /*
    *  generate flow and flow direction variables
    */
   SCIP_CALL( createFlowVariables(scip) );

   /*
    * generate pressure, pressure squared and pressure slack variables
    */
   SCIP_CALL( createPressureVariables(scip) );

   /*
    * create vars for algebraic model
    */
   if ( probdata->algebraic && probdata->network->numpipes > 0 )
   {
      SCIP_CALL( createVariablesAlgebraicModel(scip) );
   }

   /*
    * create vars for valves
    */
   if ( probdata->network->numvalves > 0 )
   {
      SCIP_CALL( createVariablesForValves(scip) );
   }

   /*
    * create vars for controlvalves
    */
   if ( probdata->network->numcontrolvalves > 0 )
   {
      SCIP_CALL( createVariablesForControlvalves(scip) );
   }

   /*
    * create variables for resistors
    */
   if ( (probdata->network->numlinresistor > 0) || (probdata->network->numnonlinresistor > 0) )
   {
      SCIP_CALL( createResistorVars(scip) );
   }

   /*
    * create variables for compressors
    */
   if ( probdata->network->numcompressor > 0 )
   {
      SCIP_CALL( createVariablesForCompressors(scip) );
   }

   return SCIP_OKAY;
}

/** generate flow constraint on arcs */
static
SCIP_RETCODE generateFlowConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char consname[SCIP_MAXSTRLEN];
   GAS_Network* Network = probdata->network;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_CONS* cons;
   SCIP_Bool removedcons = FALSE;
   int j;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   /* get memory for the variables and the coefficients */
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, probdata->network->numarcs + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, probdata->network->numarcs + 1) );

   /* generate constraint for every node */
   for (j = 0; j < Network->numnodes; j++)
   {
      GAS_Node*  Node  = &(Network->nodes_ptr[j]);
      int        pos   = 0;
      GAS_Arc*   arc;

      assert( Node != NULL );

      /* fill in incoming arcs */
      arc = Node->inarcs;
      while ( arc != NULL )
      {
         vals[pos] = -1.0;
         vars[pos++] = arc->flowvar;
         assert( arc->flowvar != NULL );
         arc = arc->next_inarc;
      }

      /* fill in outgoing arcs */
      arc = Node->outarcs;
      while ( arc != NULL )
      {
         vals[pos] = 1.0;
         vars[pos++] = arc->flowvar;
         assert( arc->flowvar != NULL );
         arc = arc->next_outarc;
      }
      assert ( pos <= probdata->network->numarcs );

      /* if maximizing the flow, create different constraints */
      if ( probdata->maxFlow )
      {
         /* if node is entry we need one additional variable (netflow) */
         if ( Node->type == ENTRY )
         {
            vals[pos] = -1.0;
            vars[pos++] = probdata->netflowvars[j];
            assert( probdata->netflowvars[j] );
         }
         /* should have found as much incident arcs as before */
         assert ( pos <= probdata->network->numarcs + 1 );

         /* create the linear constraint name */
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "flowConservation#%s#", Node->id);
         if ( Node->type == INNODE )
         {
            /* create the linear constraint and add it to SCIP */
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, pos, vars, vals, Node->flow, Node->flow) ); /* Node->flow=0 if node is innode */
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
         else if ( Node->type == EXIT )
         {
            assert( Node->scnFlowMax <= Node->scnFlowMin && Node->scnFlowMin <= 0.0 );

            /* Create the linear constraint and add it to scip */
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, pos, vars, vals, Node->scnFlowMax, Node->scnFlowMin) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
         else if ( Node->type == ENTRY )
         {
            /* create the linear constraint and add it to SCIP */
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, pos, vars, vals, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
         else
         {
            SCIPerrorMessage("Wrong node type.\n");
            SCIPABORT();
         }
      }
      /* if not maximizing the flow, create standard flow conservation constraints */
      else
      {
         /* if there is only one incident arc fix flow directly, instead of adding a constraint */
         if ( pos == 1 )
         {
            SCIP_CALL( SCIPchgVarLb(scip, vars[0], vals[0] * Node->flow) );
            SCIP_CALL( SCIPchgVarUb(scip, vars[0], vals[0] * Node->flow) );
         }
         else
         {
            /* Create the linear constraint name */
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "flowConservation#%s#", Node->id);

            /* create the linear constraint and add it to SCIP */
            if ( probdata->reduceFlowConservation && ! removedcons )
            {
               /* add constraint as non-initial, non-separable, non-enforable, but checked and propagated */
               SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, pos, vars, vals, Node->flow, Node->flow,
                     FALSE, FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
               removedcons = TRUE;
            }
            else
            {
               SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, pos, vars, vals, Node->flow, Node->flow) );
            }

            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }
   }
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &vals);

   return SCIP_OKAY;
}

/** generate flow conservation with binary variables */
static
SCIP_RETCODE generateBinaryFlowConservationConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char consname[SCIP_MAXSTRLEN];
   GAS_Network* Network;
   GAS_Node* node;
   GAS_Arc* inarc;
   GAS_Arc* outarc;
   GAS_CS* cs;
   GAS_Controlvalve* cv;
   SCIP_VAR* vars2[2] = {NULL, NULL};
   SCIP_Real vals2[2];
   SCIP_VAR** bfinvars;               /* ingoing flow direction variables (bf = binary flow) */
   SCIP_VAR** bfoutvars;              /* outgoing flow direction variables (bf = binary flow) */
   SCIP_VAR** owinvars;               /* binary variables for ingoing cs or cv (ow = one way) */
   SCIP_VAR** owoutvars;              /* binary variables for outgoing cs or cv (ow = one way) */
   SCIP_VAR** invars;                 /* array for generating linear constraints */
   SCIP_Real* invals;                 /* coefficient array */
   SCIP_VAR** outvars;                /* array for generating linear constraints */
   SCIP_Real* outvals;                /* coefficient array */
   SCIP_CONS* flowconservation;
   int numbfdirvars;
   int numowinvars;
   int numowoutvars;
   int numoutvars;
   int numinvars;
   int length;
   int i;
   int j;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   Network = probdata->network;

   length = Network->maxnumincidentarcs + 1;
   SCIP_CALL( SCIPallocClearBufferArray(scip, &bfinvars, length) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &bfoutvars, length) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &owinvars, length) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &owoutvars, length) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &invars, length) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &invals, length) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &outvars, length) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &outvals, length) );

   /* generate constraints for every node */
   for (j = 0; j < Network->numnodes; ++j)
   {
      node = &(Network->nodes_ptr[j]);
      assert( node != NULL );

      if ( node->numincidentarcs == 0 )
      {
         if ( ! SCIPisFeasZero(scip, node->flow) )
         {
            SCIPerrorMessage("Model infeasible: Found an isolated node <%s> with nonzero flow %g.\n", node->id, node->flow);
            return SCIP_INVALIDDATA;
         }
         else
            SCIPinfoMessage(scip, NULL, "Warning: Node <%s> is isolated with 0 flow requirement.\n", node->id);
      }

      /* reset variables */
      numbfdirvars = 0;
      numowinvars = 0;
      numowoutvars = 0;

      /* find flow direction variables/cs-/cv-binvars of incident arcs */
      inarc = node->inarcs;
      while ( inarc != NULL )
      {
         if ( inarc->type == CONTROLVALVE )
         {
            cv = inarc->detailed_info;
            owinvars[numowinvars++] = cv->binvar;
         }
         else if ( inarc->type == CS )
         {
            cs = inarc->detailed_info;
            owinvars[numowinvars++] = cs->compressor_binvar;
         }
         else
         {
            bfinvars[numbfdirvars] = inarc->positiveFlowBinvar;
            bfoutvars[numbfdirvars] = inarc->negativeFlowBinvar;
            ++numbfdirvars;
         }
         inarc = inarc->next_inarc;
      }

      outarc = node->outarcs;
      while ( outarc != NULL )
      {
         if ( outarc->type == CONTROLVALVE )
         {
            cv = outarc->detailed_info;
            owoutvars[numowoutvars++] = cv->binvar;
         }
         else if ( outarc->type == CS )
         {
            cs = outarc->detailed_info;
            owoutvars[numowoutvars++] = cs->compressor_binvar;
         }
         else
         {
            bfinvars[numbfdirvars] = outarc->negativeFlowBinvar;
            bfoutvars[numbfdirvars] = outarc->positiveFlowBinvar;
            ++numbfdirvars;
         }
         outarc = outarc->next_outarc;
      }

      if ( (numbfdirvars + numowinvars + numowoutvars) != node->numincidentarcs )
      {
         SCIPerrorMessage("Somewhere counting went wrong: numincidentarcs of node <%s> does not match number of incident arcs.\n", node->id);
         return SCIP_ERROR;
      }

      /* REMARK: Correctness for maxFlow objective
       *         For all objectives besides maxFlow node->flow is set to lowerflowbound from
       *         the scn file. Since in for all those objectives lower and upper bound in the
       *         scn-file are required to be equal, it actually doesn't matter which bound we
       *         choose. For the maxFlow objective we choose the upper bound instead. Due to
       *         this, case 4 does not detect entries or exits with in- or outflow equal zero.
       *         (TO DO) Rework at least case 4 and add some if case "if ( node->type == INNODE )"
       *         similar to case 1.
       */

      /*
       *  CASE 1: node->numincidentarcs == 1
       *          The node has only one neighbor. Then the gas can only flow to the node
       *          (if it is an exit) or away from the node (if it is an entry). Otherwise,
       *          it would flow in a cycle. If the node is an innode, there can be no flow.
       *          Since there is only one incident arc, we can set the corresponding outgoing
       *          binary variable to 1 if the node is an entry with positive inflow. Analog,
       *          for the ingoing binary variable, if the node is an exit.
       */
      if ( ! probdata->binVarEQone && node->numincidentarcs == 1 )
      {
         assert( numbfdirvars + numowinvars + numowoutvars == 1 );

         /* exit node: */
         if ( SCIPisNegative(scip, node->flow) )
         {
            if ( numbfdirvars == 1 )
            {
               SCIP_CALL( SCIPchgVarLb(scip, bfinvars[0], 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, bfoutvars[0], 0.0) );
            }
            else if ( numowinvars == 1 )
            {
               SCIP_CALL( SCIPchgVarLb(scip, owinvars[0], 1.0) );
            }
            else if ( numowoutvars == 1 )
            {
               SCIPerrorMessage("Wrong input data: Found an exit <%s> with only outgoing arcs.\n", node->id);
               return SCIP_INVALIDDATA;
            }
         }
         /* entry node: */
         else if ( SCIPisPositive(scip, node->flow) )
         {
            if ( numbfdirvars == 1 )
            {
               SCIP_CALL( SCIPchgVarLb(scip, bfoutvars[0], 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, bfinvars[0], 0.0) );
            }
            else if ( numowoutvars == 1 )
            {
               SCIP_CALL( SCIPchgVarLb(scip, owoutvars[0], 1.0) );
            }
            else if ( numowinvars == 1 )
            {
               SCIPerrorMessage("Wrong input data: Found an entry <%s> with only ingoing arcs.\n", node->id);
               return SCIP_INVALIDDATA;
            }
         }
         /* otherwise we have a node with zero flow: */
         else
         {
            assert( SCIPisZero(scip, node->flow) );
            if ( numbfdirvars == 1 )
            {
               SCIP_CALL( SCIPchgVarUb(scip, bfoutvars[0], 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, bfinvars[0], 0.0) );
            }
            else if ( numowoutvars == 1 )
            {
               SCIP_CALL( SCIPchgVarUb(scip, owoutvars[0], 0.0) );
            }
            else if ( numowinvars == 1 )
            {
               SCIP_CALL( SCIPchgVarUb(scip, owinvars[0], 0.0) );
            }
         }
      }
      /*
       *  CASE 2: node->numneighbors == 1
       *          The node has only one neighbor. Then the gas can only flow to the node
       *          (if it is an exit) or away from the node (if it is an entry). Otherwise,
       *          it would flow in a cycle. If the node is an innode, there can be no flow.
       */
      else if ( ! probdata->binVarEQone && node->numneighbors == 1 )
      {
         assert( numbfdirvars + numowinvars + numowoutvars == node->numneighbors );

         /* if node is an entry or innode, no flow can go into the node */
         if ( SCIPisGE(scip, node->flow, 0.0) )
         {
            for (i = 0; i < numbfdirvars; ++i)
            {
               SCIP_CALL( SCIPchgVarUb(scip, bfinvars[i], 0.0) );
            }
            for (i = 0; i < numowinvars; ++i)
            {
               SCIP_CALL( SCIPchgVarUb(scip, owinvars[i], 0.0) );
            }
         }

         /* if node is an exit or innode, no flow can go out from the node */
         if ( SCIPisLE(scip, node->flow, 0.0) )
         {
            for (i = 0; i < numbfdirvars; ++i)
            {
               SCIP_CALL( SCIPchgVarUb(scip, bfoutvars[i], 0.0) );
            }
            for (i = 0; i < numowoutvars; ++i)
            {
               SCIP_CALL( SCIPchgVarUb(scip, owoutvars[i], 0.0) );
            }
         }
      }
      /*
       *  CASE 3: node->numincidentarcs == 2 && SCIPisZero(scip, node->flow)
       *          The "innode" has exactly two neighbors, since one neighbor was the
       *          case before. Hence, we can couple the binary (direction) variables.
       *
       *          The case of entries and exits with two incident arcs will be treated
       *          in the general case later on.
       */
      else if ( node->numincidentarcs == 2 && SCIPisZero(scip, node->flow) )
      {
         if ( numowinvars == 2 )
         {
            SCIP_CALL( SCIPchgVarUb(scip, owinvars[0], 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, owinvars[1], 0.0) );
         }
         else if ( numowoutvars == 2 )
         {
            SCIP_CALL( SCIPchgVarUb(scip, owoutvars[0], 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, owoutvars[1], 0.0) );
         }
         else if ( numbfdirvars == 2 )
         {
            vars2[0] = bfinvars[0];
            vars2[1] = bfoutvars[1];
            vals2[0] = 1.0;
            vals2[1] = - 1.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowConsA#%s#", node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, 2, vars2, vals2, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, flowconservation) );
            SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );

            vars2[0] = bfinvars[1];
            vars2[1] = bfoutvars[0];
            vals2[0] = 1.0;
            vals2[1] = - 1.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowConsB#%s#", node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, 2, vars2, vals2, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, flowconservation) );
            SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );
         }
         else if ( numowinvars > 0 && numowoutvars > 0 )
         {
            assert( numowinvars == 1 && numowoutvars == 1 );

            vars2[0] = owinvars[0];
            vars2[1] = owoutvars[0];
            vals2[0] = 1.0;
            vals2[1] = - 1.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowCons#%s#", node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, 2, vars2, vals2, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, flowconservation) );
            SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );
         }
         else if ( numowinvars > 0 && numbfdirvars > 0 )
         {
            assert( numowinvars == 1 && numbfdirvars == 1 );

            vars2[0] = owinvars[0];
            vars2[1] = bfoutvars[0];
            vals2[0] = 1.0;
            vals2[1] = - 1.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowCons#%s#", node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, 2, vars2, vals2, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, flowconservation) );
            SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );

            SCIP_CALL( SCIPchgVarUb(scip, bfinvars[0], 0.0) );
         }
         else if ( numbfdirvars > 0 && numowoutvars > 0 )
         {
            assert( numowoutvars == 1 && numbfdirvars == 1 );

            vars2[0] = bfinvars[0];
            vars2[1] = owoutvars[0];
            vals2[0] = 1.0;
            vals2[1] = - 1.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowCons#%s#", node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, 2, vars2, vals2, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, flowconservation) );
            SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );

            SCIP_CALL( SCIPchgVarUb(scip, bfoutvars[0], 0.0) );
         }
         else
         {
            SCIPerrorMessage("This case in generateBinaryFlowConservationConstraints should not occurr!\n");
            return SCIP_ERROR;
         }
      }
      /*
       *  CASE 4: SCIPisZero(scip, node->flow)
       */
      else if ( ! probdata->binVarEQone && SCIPisZero(scip, node->flow) )
      {
         int cntin = 1;
         int cntout = 1;

         numinvars = numbfdirvars + numowinvars + 1;
         numoutvars = numbfdirvars + numowoutvars + 1;

         /* initialize */
         outvals[0] = -1.0;
         invals[0] = -1.0;

         /* set variables for constraints of the form: out <= in1 + in2 + ... and in <= out1 + out2 + ... */
         for (i = 0; i < numbfdirvars; ++i)
         {
            invars[cntin] = bfinvars[i];
            invals[cntin] = 1.0;
            ++cntin;

            outvars[cntout] = bfoutvars[i];
            outvals[cntout] = 1.0;
            ++cntout;
         }
         for (i = 0; i < numowoutvars; ++i)
         {
            outvars[cntout] = owoutvars[i];
            outvals[cntout] = 1.0;
            ++cntout;
         }
         for (i = 0; i < numowinvars; ++i)
         {
            invars[cntin] = owinvars[i];
            invals[cntin] = 1.0;
            ++cntin;
         }
         assert( cntin == numinvars );
         assert( cntout == numoutvars );
         assert( cntin <= Network->maxnumincidentarcs + 1 );
         assert( cntout <= Network->maxnumincidentarcs + 1 );

         /* generate constraints */
         for (i = 0; i < numbfdirvars; ++i)
         {
            invars[0] = bfoutvars[i];
            outvars[0] = bfinvars[i];

            /* temporarily overwrite coefficient and restore later */
            invals[i + 1]  = 0.0;
            outvals[i + 1] = 0.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowConsID#%d#%s#", i, node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, numinvars, invars, invals, 0.0, SCIPinfinity(scip)) );
            SCIP_CALL( SCIPaddCons(scip, flowconservation) );
            SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowConsOD#%d#%s#", i, node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, numoutvars, outvars, outvals, 0.0, SCIPinfinity(scip)) );
            SCIP_CALL( SCIPaddCons(scip, flowconservation) );
            SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );

            /* restore value */
            invals[i+1]  = 1.0;
            outvals[i+1] = 1.0;
         }
         for (i = 0; i < numowoutvars; ++i)
         {
            invars[0]  = owoutvars[i];

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowConsO#%d#%s#", i, node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, numinvars, invars, invals, 0.0, SCIPinfinity(scip)) );
            SCIP_CALL( SCIPaddCons(scip, flowconservation) );
            SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );
         }
         for (i = 0; i < numowinvars; ++i)
         {
            outvars[0]  = owinvars[i];

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowConsI#%d#%s#", i, node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, numoutvars, outvars, outvals, 0.0, SCIPinfinity(scip)) );
            SCIP_CALL( SCIPaddCons(scip, flowconservation) );
            SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );
         }
      }
      /*
       *  CASE 5: exit: SCIPisNegative(scip, node->flow)
       *          The node is an exit and there is outflow. Then at least one arc has to be directed
       *          to the exit.
       */
      else if ( SCIPisNegative(scip, node->flow) )
      {
         int cnt = 0;
         for (i = 0; i < numbfdirvars; ++i)
         {
            invals[cnt] = 1.0;
            invars[cnt++] = bfinvars[i];
         }
         for (i = 0; i < numowinvars; ++i)
         {
            invals[cnt] = 1.0;
            invars[cnt++] = owinvars[i];
         }
         assert( cnt == numbfdirvars + numowinvars );

         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowCons#%s#", node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, cnt, invars, invals, 1.0, SCIPinfinity(scip)) );
         SCIP_CALL( SCIPaddCons(scip, flowconservation) );
         SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );
      }
      /*
       *  CASE 6: entry: SCIPisPositive(scip, node->flow)
       *          The node is an entry and there is inflow. Then at least one arc has to be directed
       *          away from the node.
       */
      else if ( SCIPisPositive(scip, node->flow) )
      {
         int cnt = 0;
         for (i = 0; i < numbfdirvars; ++i)
         {
            outvals[cnt] = 1.0;
            outvars[cnt++] = bfoutvars[i];
         }
         for (i = 0; i < numowoutvars; ++i)
         {
            outvals[cnt] = 1.0;
            outvars[cnt++] = owoutvars[i];
         }
         assert( cnt == numbfdirvars + numowoutvars );

         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "binaryFlowCons#%s#", node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowconservation, consname, cnt, outvars, outvals, 1.0, SCIPinfinity(scip)) );
         SCIP_CALL( SCIPaddCons(scip, flowconservation) );
         SCIP_CALL( SCIPreleaseCons(scip, &flowconservation) );
      }
      else if ( ! probdata->binVarEQone )
      {
         /* SCIPerrorMessage("Either something is wrong with the data or I forgot to implement some case!\n");
          * return SCIP_ERROR;
         */
      }
   }   /* end for all nodes */

   SCIPfreeBufferArray(scip, &outvals);
   SCIPfreeBufferArray(scip, &outvars);
   SCIPfreeBufferArray(scip, &invals);
   SCIPfreeBufferArray(scip, &invars);
   SCIPfreeBufferArray(scip, &owoutvars);
   SCIPfreeBufferArray(scip, &owinvars);
   SCIPfreeBufferArray(scip, &bfoutvars);
   SCIPfreeBufferArray(scip, &bfinvars);

   return SCIP_OKAY;
}

/** Generate the constraints that couple the flow with the binary variable for the direction. */
static
SCIP_RETCODE generateFlowDirectionConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   GAS_Network* Network;
   GAS_Arc*     arc;
   SCIP_VAR*    vars[2];
   SCIP_Real    vals[2];
   SCIP_VAR*    vars3[3];
   SCIP_Real    vals3[3];
   SCIP_Real    lhs;
   SCIP_Real    rhs;
   SCIP_Real    qlb;
   SCIP_Real    qub;
   SCIP_CONS*   binvarcoupling;
   SCIP_CONS*   negativeflowcons;
   SCIP_CONS*   positiveflowcons;
   SCIP_CONS*   pressureDiffCons;
   char         consname[SCIP_MAXSTRLEN];
   int          j;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->network != NULL );

   Network = probdata->network;

   for (j = 0; j < probdata->network->numarcs; ++j)
   {
      arc = &Network->arcs_ptr[j];

      /* we do not have a binary variable if the arc is a CS or control valve */
      if ( arc->type == CS || arc->type == CONTROLVALVE )
         continue;

      qlb = arc->flowMin;
      qub = arc->flowMax;

      /* add coupling of the binary variables
       * 0.0 <= positiveFlowBinvar + negativeFlowBinvar <= 1.0
         if option binVarEQone is set then: positiveFlowBinvar + negativeFlowBinvar = 1.0 */
      if (SCIPisFeasNegative(scip, qlb) && SCIPisFeasPositive(scip, qub))
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "flowDirBinvars#%s#", arc->id);

         vars[0] = arc->positiveFlowBinvar;
         vars[1] = arc->negativeFlowBinvar;
         vals[0] = 1.0;
         vals[1] = 1.0;
         rhs     = 1.0;

         if ( probdata->binVarEQone )
            lhs = 1.0;
         else
            lhs = 0.0;

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &binvarcoupling, consname, 2, vars, vals, lhs, rhs) );
         SCIP_CALL( SCIPaddCons(scip, binvarcoupling) );
         SCIP_CALL( SCIPreleaseCons(scip, &binvarcoupling) );
      }

      /* add negativeFlowCons:
       * 0.0 <= q_a - q_lb * negativeFlowBinvar <= + SCIPinfinity
       * or fix binary variables if flow bounds are already nonpositive or nonnegative */
      if ( SCIPisFeasGE(scip, qlb, 0.0) )
      {
         SCIP_CALL( SCIPchgVarUb(scip, arc->negativeFlowBinvar, 0.0) );
      }
      else if ( SCIPisFeasNegative(scip, qub) )
      {
         SCIP_CALL( SCIPchgVarLb(scip, arc->negativeFlowBinvar, 1.0) );
      }
      else
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "negFlowDir#%s#", arc->id);

         lhs     = 0.0;
         rhs     = SCIPinfinity(scip);

         SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &negativeflowcons, consname, arc->flowvar, arc->negativeFlowBinvar, - qlb, lhs, rhs) );
         SCIP_CALL( SCIPaddCons(scip, negativeflowcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &negativeflowcons) );
      }

      /* add positiveFlowCons:
       * - SCIPinfinity <= q_a - q_ub * positiveFlowBinvar <= 0.0
       * or fix binary variables if flow bounds are already nonpositive or nonnegative */
      if ( SCIPisFeasLE(scip, qub, 0.0) )
      {
         SCIP_CALL( SCIPchgVarUb(scip, arc->positiveFlowBinvar, 0.0) );
      }
      else if ( SCIPisFeasPositive(scip, qlb) )
      {
         SCIP_CALL( SCIPchgVarLb(scip, arc->positiveFlowBinvar, 1.0) );
      }
      else
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "posFlowDir#%s#", arc->id);

         lhs     = - SCIPinfinity(scip);
         rhs     = 0.0;

         SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &positiveflowcons, consname, arc->flowvar, arc->positiveFlowBinvar, - qub, lhs, rhs) );
         SCIP_CALL( SCIPaddCons(scip, positiveflowcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &positiveflowcons) );
      }

      /* couple the flow direction binvar with the pressure variables for pipes */
      if ( arc->type == PIPE )
      {
         if ( ! probdata->relaxLowerBounds && ! probdata->relaxUpperBounds )
         {
            vars3[0] = Network->arcs_ptr[j].sourcenode->pressurevar;
            vars3[1] = Network->arcs_ptr[j].targetnode->pressurevar;
            vars3[2] = arc->positiveFlowBinvar;
            vals3[0] = 1.0;
            vals3[1] = - 1.0;

            /* if flow is positive, then the pressure at the target has to be lower:
             * p_source - p_target - positiveFlowBinvar * (p_source^ub - p_target^lb) <= 0.0
             */
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PIPE_posPressureDiff#%s#", arc->id);

            vals3[2] = - (Network->arcs_ptr[j].sourcenode->pressureMax - Network->arcs_ptr[j].targetnode->pressureMin);
            lhs = - SCIPinfinity(scip);
            rhs = 0.0;

            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &pressureDiffCons, consname, 3, vars3, vals3, lhs, rhs) );
            SCIP_CALL( SCIPaddCons(scip, pressureDiffCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &pressureDiffCons) );

            /* if flow is negative, then the pressure at the target has to be bigger:
             * 0.0 <= p_source - p_target - negativeFlowBinvar * (p_source^lb - p_target^ub)
             */
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PIPE_negPressureDiff#%s#", arc->id);

            vars3[2] = arc->negativeFlowBinvar;
            vals3[2] = - (Network->arcs_ptr[j].sourcenode->pressureMin - Network->arcs_ptr[j].targetnode->pressureMax);
            lhs = 0.0;
            rhs = SCIPinfinity(scip);

            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &pressureDiffCons, consname, 3, vars3, vals3, lhs, rhs) );
            SCIP_CALL( SCIPaddCons(scip, pressureDiffCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &pressureDiffCons) );
         }
      }
      else if ( arc->type == VALVE )
      {
         GAS_Valve* valve;
         SCIP_CONS* cons;

         valve = arc->detailed_info;

         if ( ! probdata->binVarEQone )
         {
            /* couple valve and flow directions: positiveFlowBinvar + negativeFlowBinvar - valve_binvar == 0.0 */
            vars3[0] = arc->positiveFlowBinvar;
            vars3[1] = arc->negativeFlowBinvar;
            vars3[2] = valve->valve_binvar;

            vals3[0] = 1.0;
            vals3[1] = 1.0;
            vals3[2] = -1.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "VALVE_BinvarCoupling#%s#", arc->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 3, vars3, vals3, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
         else
         {
            qlb = arc->flowMin;
            qub = arc->flowMax;

            if ( SCIPisPositive(scip, qub) )
            {
               /* q_a - q_ub * valveBinvar <= 0.0 */
               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "VALVE_BinvarCoupling#%s#ub", arc->id);
               SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, arc->flowvar, valve->valve_binvar, -qub, - SCIPinfinity(scip), 0.0) );
               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }

            if ( SCIPisNegative(scip, qlb) )
            {
               /* 0 <= q_a - q_lb * valveBinvar */
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "VALVE_BinvarCoupling#%s#lb", arc->id);
               SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, arc->flowvar, valve->valve_binvar, -qlb, 0.0, SCIPinfinity(scip)) );
               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }
         }
      }
      else if ( arc->type == RESISTOR )
      {
         GAS_Resistor* resistor;
         SCIP_CONS* resistorcons;

         resistor = arc->detailed_info;

         if ( resistor->type == LINEAR )
         {
            /* LR_posFlowDir (for flow greater than some q_epsilon) can only be 1 if positiveFlowBinvar is 1 */
            vars[0] = arc->positiveFlowBinvar;
            vars[1] = probdata->LR_posFlowDir[resistor->LR_position];

            vals[0] = - 1.0;
            vals[1] = 1.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LR_BinvarCoupling1#%s#", arc->id);

            lhs = - 1.0;
            rhs = 0.0;

            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &resistorcons, consname, 2, vars, vals, lhs, rhs) );
            SCIP_CALL( SCIPaddCons(scip, resistorcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &resistorcons) );

            /* LR_negFlowDir (for flow smaller than some -q_epsilon) can only be 1 if negativeFlowBinvar is 1 */
            vars[0] = arc->negativeFlowBinvar;
            vars[1] = probdata->LR_negFlowDir[resistor->LR_position];

            vals[0] = - 1.0;
            vals[1] = 1.0;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LR_BinvarCoupling2#%s#", arc->id);

            lhs = -1.0;
            rhs = 0.0;

            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &resistorcons, consname, 2, vars, vals, lhs, rhs) );
            SCIP_CALL( SCIPaddCons(scip, resistorcons) );
            SCIP_CALL( SCIPreleaseCons(scip, &resistorcons) );
         }
         /* eventually add something for the nonlinear resistors */
      }
   }

   return SCIP_OKAY;
}

/** generate circular flow constraints
 *
 *  Gas cannot flow in cycles!
 *  For each cycle we found, we introduce two cuts:
 *
 *  positive Cycle: cut forbids a circular flow with positive direction
 *  negative Cycle: cut forbids a circular flow with negative direction
 *
 *  Note that controlvalves and compressores can only be passed in one
 *  direction. Hence we might only add one or no cut at all, when cvs or css
 *  are part of the cycle.
 */
static
SCIP_RETCODE generateNoCircularFlowConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   GAS_Network*      network;
   GAS_Valve*        valve;
   GAS_Controlvalve* cv;
   GAS_CS*           cs;
   GAS_Cycle*        cycle;
   SCIP_VAR**        posvars;                               /* variable array for positive cycle */
   SCIP_VAR**        negvars;                               /* variable array for negative cycle */
   SCIP_Real*        vals;                                  /* coefficient array */
   SCIP_CONS*        nocircularflowcons;
   char              consname[SCIP_MAXSTRLEN];
   int               nposvars;                              /* variables in the positive cycle cut */
   int               nnegvars;                              /* variables in the negative cycle cut */
   int               prhs;                                  /* rhs for positive flow direction cut */
   int               nrhs;                                  /* rhs for negative flow direction cut */
   int               i;
   int               ncons;                                 /* total number of generated noCircularFlow constraints */
   SCIP_Bool         posCycle;                              /* wheter or not the positive cirlce cut should be created */
   SCIP_Bool         negCycle;                              /* wheter or not the negative cirlce cut should be created */

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->network != NULL );
   assert( probdata->network->firstCycle != NULL );
   assert( probdata->network->lastCycle != NULL );

   network = probdata->network;
   cycle   = network->firstCycle;
   ncons   = 0;

   SCIP_CALL( SCIPallocClearBufferArray(scip, &posvars, network->numarcs) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &negvars, network->numarcs) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &vals, network->numarcs) );

   /* array coefficients are always 1.0 */
   for (i = 0; i < network->numarcs; ++i)
      vals[i] = 1.0;

   while ( cycle != NULL )
   {
      posCycle = TRUE;
      negCycle = TRUE;
      prhs      = -1;
      nrhs      = -1;
      nposvars  = 0;
      nnegvars  = 0;

      /* set variable arrays! */
      for (i = 0; i < network->numarcs; ++i)
      {
         if ( cycle->arcsInCycle[i] == 0 )
         {
            continue;
         }
         else
         {
            if ( network->arcs_ptr[i].type == CS )
            {
               cs = network->arcs_ptr[i].detailed_info;
               if ( cycle->arcsInCycle[i] == 1 )
               {
                  posvars[nposvars] = cs->compressor_binvar;
                  ++prhs;
                  ++nposvars;

                  /* the compressor cannot be passed in negative direction */
                  negCycle = FALSE;
               }
               else                                         /* cycle->arcsInCycle[i] == -1 */
               {
                  negvars[nnegvars] = cs->compressor_binvar;
                  ++nrhs;
                  ++nnegvars;

                  /* the compressor cannot be passed in negative direction */
                  posCycle = FALSE;
               }
            }
            else if ( network->arcs_ptr[i].type == CONTROLVALVE )
            {
               cv = network->arcs_ptr[i].detailed_info;
               if ( cycle->arcsInCycle[i] == 1 )
               {
                  posvars[nposvars] = cv->binvar;
                  ++prhs;
                  ++nposvars;

                  /* the controlvalve cannot be passed in negative direction */
                  negCycle = FALSE;
               }
               else
               {
                  negvars[nnegvars] = cv->binvar;
                  ++nrhs;
                  ++nnegvars;

                  /* the controlvalve cannot be passed in negative direction */
                  posCycle = FALSE;
               }
            }
            else if ( network->arcs_ptr[i].type == VALVE )
            {
               valve = network->arcs_ptr[i].detailed_info;
               if ( cycle->arcsInCycle[i] == 1 )
               {
                  posvars[nposvars] = network->arcs_ptr[i].positiveFlowBinvar;
                  ++nposvars;
                  ++prhs;
                  negvars[nnegvars] = network->arcs_ptr[i].negativeFlowBinvar;
                  ++nnegvars;
                  ++nrhs;

                  /* if we bypassed a cv or cs we can add the corresponding variable.
                   * we do not need to alter the rhs since either valve or cv/cs can
                   * be open or on
                   */
                  if ( valve->bypassed_cv != NULL )
                  {
                     cv = valve->bypassed_cv->detailed_info;
                     posvars[nposvars] = cv->binvar;
                     ++nposvars;
                  }
                  else if ( valve->bypassed_cs != NULL )
                  {
                     cs = valve->bypassed_cs->detailed_info;
                     posvars[nposvars] = cs->compressor_binvar;
                     ++nposvars;
                  }
               }
               else
               {
                  posvars[nposvars] = network->arcs_ptr[i].negativeFlowBinvar;
                  ++nposvars;
                  ++prhs;
                  negvars[nnegvars] = network->arcs_ptr[i].positiveFlowBinvar;
                  ++nnegvars;
                  ++nrhs;

                  /* if we bypassed a cv or cs we can add the corresponding variable.
                   * we do not need to alter the rhs since either valve or cv/cs can
                   * be open or on
                   */
                  if ( valve->bypassed_cv != NULL )
                  {
                     cv = valve->bypassed_cv->detailed_info;
                     negvars[nnegvars] = cv->binvar;
                     ++nnegvars;
                  }
                  else if ( valve->bypassed_cs != NULL )
                  {
                     cs = valve->bypassed_cs->detailed_info;
                     negvars[nnegvars] = cs->compressor_binvar;
                     ++nnegvars;
                  }
               }
            }
            else
            {
              if ( cycle->arcsInCycle[i] == 1 )
               {
                  posvars[nposvars] = network->arcs_ptr[i].positiveFlowBinvar;
                  ++nposvars;
                  ++prhs;
                  negvars[nnegvars] = network->arcs_ptr[i].negativeFlowBinvar;
                  ++nnegvars;
                  ++nrhs;
               }
               else
               {
                  posvars[nposvars] = network->arcs_ptr[i].negativeFlowBinvar;
                  ++nposvars;
                  ++prhs;
                  negvars[nnegvars] = network->arcs_ptr[i].positiveFlowBinvar;
                  ++nnegvars;
                  ++nrhs;
               }
            }
         }

         if ( !posCycle && !negCycle )
         {
            break;
         }
      }

      if ( posCycle )
      {
         ++ncons;
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "noCircularFlowCons#%d", ncons);

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &nocircularflowcons, consname, nposvars, posvars, vals, - SCIPinfinity(scip), (SCIP_Real) prhs) );
         SCIP_CALL( SCIPaddCons(scip, nocircularflowcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &nocircularflowcons) );
      }
      if ( negCycle )
      {
         ++ncons;
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "noCircularFlowCons#%d", ncons);

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &nocircularflowcons, consname, nnegvars, negvars, vals, - SCIPinfinity(scip), (SCIP_Real) nrhs) );
         SCIP_CALL( SCIPaddCons(scip, nocircularflowcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &nocircularflowcons) );
      }

      cycle = cycle->next_cycle;
   }

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &negvars);
   SCIPfreeBufferArray(scip, &posvars);

   return SCIP_OKAY;
}

/** generate constraints of the form PI = p^2 for nonlinear resistors and weymouth equation */
static
SCIP_RETCODE generatePressureSquaredConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   GAS_Node* node;
   SCIP_VAR* vars1[1];
   SCIP_VAR* vars2[1];
   SCIP_Real vals1[1];
   SCIP_Real vals2[1];
   SCIP_CONS* constraint;
   char consname[SCIP_MAXSTRLEN];
   int i;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->network != NULL );

   for (i = 0; i < (probdata->network->numnodes); ++i)
   {
      node = &(probdata->network->nodes_ptr[i]);

      if ( node->needPIvar )
      {
         vars1[0] = node->PIvar;
         vars2[0] = node->pressurevar;
         vals1[0] = 1;
         vals2[0] = -1;

         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "SquaredPressureCons#%s#", node->id);
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
         SCIP_CALL( SCIPcreateConsBasicQuadraticNonlinear(scip, &constraint, consname, 1,  vars1, vals1, 1, vars2, vars2, vals2, 0.0, 0.0) );
#else
         SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &constraint, consname, 1,  vars1, vals1, 1, vars2, vars2, vals2, 0.0, 0.0) );
#endif
         SCIP_CALL( SCIPaddCons(scip, constraint) );
         SCIP_CALL( SCIPreleaseCons(scip, &constraint) );
      }
   }

   return SCIP_OKAY;
}

/** generates the constraints for the ODE gas flow model (ISO1) */
static
SCIP_RETCODE generateODEModel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   arcposition         /**< position in arc array */
   )
{
   GAS_Arc* arc;
   GAS_Pipe* pipe;
   GAS_Node* sourcenode;
   GAS_Node* targetnode;
   SCIP_Real vals[2];
   SCIP_VAR* vars[2];
   SCIP_CONS* machbound;
   SCIP_CONS* PipeODE_Cons;
   char consname[SCIP_MAXSTRLEN];
   SCIP_Real A;                              /* cross sectional are of pipe*/
   int N = -1;                               /* Number of discretization steps necessary to fulfill bound on step size */
   SCIP_Real lambda;                         /* friction coefficient of the pipe */
   SCIP_Real c;                              /* speed of sound */
   int mach_numerator;
   SCIP_Real p_m;
   SCIP_Real z_m;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   arc = &(probdata->network->arcs_ptr[arcposition]);
   pipe = arc->detailed_info;
   sourcenode = arc->sourcenode;
   targetnode = arc->targetnode;

   A = pi * pow((pipe->diameter)/2.0, 2.0);

   p_m = ComputeMeanPressure(sourcenode->pressureMin, sourcenode->pressureMax, targetnode->pressureMin, targetnode->pressureMax);
   z_m = MeanCompressibilityFactor(probdata->network, p_m, probdata->papay);
   c = computeSpeedOfSound( probdata->network->gasTemperature, z_m, probdata->network->molarMass1 );

   /*
    *  Generate the MachBoundConstraints:
    *  (1) mach_numerator * A * p_in  * 1e5 => - 5 * c * q  ---> 0.0 <= mach_numerator * A * p_in  * 1e5 + 5 * c * q
    *  (2) mach_numerator * A * p_out * 1e5 =>   5 * c * q  ---> 0.0 <= mach_numerator * A * p_out * 1e5 - 5 * c * q
    */
   SCIP_CALL( SCIPgetIntParam(scip, "machbound", &mach_numerator) );

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "MachBound_negFlow#%s#", arc->id);

   vars[0] = sourcenode->pressurevar;
   vars[1] = arc->flowvar;
   vals[0] = mach_numerator * A * 1e5;
   vals[1] = 5 * c;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &machbound, consname, 2, vars, vals, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, machbound) );
   SCIP_CALL( SCIPreleaseCons(scip, &machbound) );

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "MachBound_posFlow#%s#", arc->id);

   vars[0] = targetnode->pressurevar;
   vals[1] = - 5 * c;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &machbound, consname, 2, vars, vals, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, machbound) );
   SCIP_CALL( SCIPreleaseCons(scip, &machbound) );

   /*
    *  generate the pipeODE constraint
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PipeODE#%s#", probdata->network->arcs_ptr[arcposition].id);

   lambda = NikuradzesEquation( pipe->diameter, pipe->roughness);

   SCIPdebugMessage("lambda: <%f>\n", lambda);
   if ( mach_numerator == 1 )
   {
      N = (int) ( lambda * pipe->length / ( 29.15 * pipe->diameter ) ) + 1;
   }
   else if ( mach_numerator == 2 )
   {
      N = (int) ( lambda * pipe->length / ( 4.925 * pipe->diameter ) ) + 1;
   }
   else if ( mach_numerator == 3 )
   {
      SCIPerrorMessage("TODO: Determine maximal step size for numerical methods in case v/c <= 0.6\n");

      return SCIP_ERROR;
   }
   else if ( mach_numerator == 4 )
   {
      N = (int) (100 * lambda * pipe->length / ( 16 * pipe->diameter ) ) + 1;
   }
   assert( N >= 0 );

   SCIP_CALL( SCIPcreateConsBasicPipeODE(scip, &PipeODE_Cons, consname, sourcenode->pressurevar, targetnode->pressurevar, probdata->flowvars[arcposition],
         probdata->network->arcs_ptr[arcposition].positiveFlowBinvar,
         probdata->network->arcs_ptr[arcposition].negativeFlowBinvar,
         N, A, pipe->length, pipe->diameter, lambda, c) );
   SCIP_CALL( SCIPaddCons(scip, PipeODE_Cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &PipeODE_Cons) );

   return SCIP_OKAY;
}


/** generates the constraints for the algebraic gas flow model (ISO4 / Weymouth) */
static
SCIP_RETCODE generateAlgebraicModel(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   arcposition         /**< position in arc array */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Arc* arc;
   GAS_Pipe* pipe;
   GAS_Node* sourcenode;
   GAS_Node* targetnode;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];
   SCIP_CONS* PI_diff_cons;
   SCIP_CONS* WeymouthEquality_Cons;
   char consname[SCIP_MAXSTRLEN];
   SCIP_Real LAMBDA;
   SCIP_Real meanpressure;

   assert( scip != NULL );

   probdata =  SCIPgetProbData(scip);

   arc = &(probdata->network->arcs_ptr[arcposition]);
   pipe = arc->detailed_info;
   sourcenode = arc->sourcenode;
   targetnode = arc->targetnode;

   assert( sourcenode->PIvar != NULL );
   assert( targetnode->PIvar != NULL );

   /*
    *  generate the coupling between the squared pressure variables and the pressure diff variable
    *  PI_diff - (PI_source - PI_target) = 0.0
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PIdiffCons#%s#", arc->id);

   vars[0] = pipe->PIdiffVar;
   vars[1] = sourcenode->PIvar;
   vars[2] = targetnode->PIvar;

   vals[0] =  1;
   vals[1] = -1;
   vals[2] =  1;

   /* Create linear constraint and add it to scip */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &PI_diff_cons, consname, 3, vars, vals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, PI_diff_cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &PI_diff_cons) );

   /*
    *  Generate the Weymouth equation
    */

   meanpressure = ComputeMeanPressure(sourcenode->pressureMin, sourcenode->pressureMax, targetnode->pressureMin, targetnode->pressureMax);

   if ( probdata->mixing || probdata->linearMixing )
   {
#if SCIP_VERSION >= 800
      char name[SCIP_MAXSTRLEN];
      SCIP_CONS* newMixWeymouth;
      SCIP_EXPR* exprMu;
      SCIP_EXPR* exprMol;
      SCIP_EXPR* exprFrac = NULL;
      SCIP_EXPR* exprCsqrdQ;
      SCIP_EXPR* exprFIN;
      SCIP_EXPR* expr[3];
      SCIP_Real  coef[2];
      SCIP_CONS* Lambda;
      SCIP_EXPR* exprQsqrd;
      SCIP_EXPR* exprQ;
      SCIP_EXPR* exprPIdiff;
      SCIP_EXPR* expCsqrd;

      /* set coefficients */
      SCIP_Real LambdaMix = computeWeymouthConstantMixing(probdata->network, pipe->length, pipe->diameter, pipe->roughness, meanpressure, probdata->papay);
      SCIP_Real molarMass1 = probdata->molarMass1;
      SCIP_Real molarMass2 = probdata->molarMass2;
      SCIP_Real coefM = molarMass1 - molarMass2;

#ifdef USE_ADIABATIC_CAPACITY
      /* if the adiabatic/specific heat capacity should be taken into account for speed of sound calculation uncommend this */
      SCIP_Real specificHeatCap1 = probdata->specificHeatCap1;
      SCIP_Real specificHeatCap2 = probdata->specificHeatCap2;
      SCIP_Real coefG = specificHeatCap1 - specificHeatCap2;
      SCIP_EXPR* exprGamma;
#endif

      /* create expr variables */
      assert( probdata->flowvars[arcposition] != NULL );
      SCIP_CALL( SCIPcreateExprVar(scip, &exprQ, probdata->flowvars[arcposition], NULL, NULL) );
      SCIP_CALL( SCIPcreateExprSignpower(scip, &exprQsqrd, exprQ, 2.0, NULL, NULL) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exprQ) );
      assert( pipe->PIdiffVar != NULL );
      SCIP_CALL( SCIPcreateExprVar(scip, &exprPIdiff, pipe->PIdiffVar, NULL, NULL) );

      if ( ! probdata->linearMixing )
      {
         /* in the following the big Lambda coefficient of the Weymouth Eqn is built and then the Weym. Eqn. c= RTz/M(mu) */
         assert( arc->mixingRatio != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &exprMu, arc->mixingRatio, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &exprMu, &coefM, molarMass2, NULL, NULL) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprMu) );
         SCIP_CALL( SCIPcreateExprPow(scip, &exprFrac, exprMol, -1.0, NULL, NULL) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprMol) );
      }

#ifdef USE_ADIABATIC_CAPACITY
      SCIP_CALL( SCIPcreateExprSum(scip, &exprGamma, 1, &exprMu, &coefG, specificHeatCap2, NULL, NULL) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exprMu) );
      expr[0] = exprGamma;
      expr[1] = exprFrac;
      SCIP_CALL( SCIPcreateExprProduct(scip, &expGMRTz, 2, expr, LambdaMix, NULL, NULL) ); /* Would be c^2 if LamdaMix is replaced by RTz */
#endif

      if ( ! probdata->linearMixing )
      {
         assert( exprFrac != NULL );
         assert( arc->Lambda != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &expCsqrd, arc->Lambda, NULL, NULL) );
         expr[0] = expCsqrd;
         expr[1] = exprFrac;
         coef[0] = 1.0;
         coef[1] = -LambdaMix;
         SCIP_CALL( SCIPcreateExprSum(scip, &exprFIN, 2, expr, coef, 0.0, NULL, NULL) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expCsqrd) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprFrac) );

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "WeymouthConstant%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &Lambda, name, exprFIN, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, Lambda) );
         SCIP_CALL( SCIPreleaseCons(scip, &Lambda) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprFIN) );
      }
      else
      {

         /* Build the linear mass% convex comb of the isolated speed of sounds */
         assert( arc->mixingRatioMass != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &exprMu, arc->mixingRatioMass, NULL, NULL) );

         /* (1/M1 -1/M2)* mu^mass + 1/M2 */
         coefM = 1.0 / molarMass1 - 1.0 / molarMass2;

         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &exprMu, &coefM, 1.0 / molarMass2, NULL, NULL) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprMu) );

         /* LambdaMix * [(1/M1 -1/M2)* mu^mass + 1/M2 ]   - LAMBDA =0*/
         assert( arc->Lambda != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &expCsqrd, arc->Lambda, NULL, NULL) );
         coef[0] = LambdaMix;
         coef[1] = -1.0;
         expr[0] = exprMol;
         expr[1] = expCsqrd;
         SCIP_CALL( SCIPcreateExprSum(scip, &exprFIN, 2, expr, coef, 0.0, NULL, NULL) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprMol) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expCsqrd) );

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "WeymouthConstant%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &Lambda, name, exprFIN, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, Lambda) );
         SCIP_CALL( SCIPreleaseCons(scip, &Lambda) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprFIN) );
      }

      /* here the Weymouth EQN is built */
      assert( arc->Lambda != NULL );
      SCIP_CALL( SCIPcreateExprVar(scip, &expCsqrd, arc->Lambda, NULL, NULL) );

      expr[0] = expCsqrd;
      expr[1] = exprQsqrd;

      /* if RTz is used above to calculate the real SoS then here Lamda instead of 1.0 needs to be used */
      SCIP_CALL( SCIPcreateExprProduct(scip, &exprCsqrdQ, 2, expr, 1.0, NULL, NULL) );
      SCIP_CALL( SCIPreleaseExpr(scip, &expCsqrd) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exprQsqrd) );

      expr[0] = exprPIdiff;
      expr[1] = exprCsqrdQ;
      coef[0] = 1.0;
      coef[1] = -1.0;
      SCIP_CALL( SCIPcreateExprSum(scip, &exprFIN, 2, expr, coef, 0.0, NULL, NULL) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exprPIdiff) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exprCsqrdQ) );

      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "MixWeymouth%s", arc->id);
      SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &newMixWeymouth, name, exprFIN, 0.0, 0.0) );
      SCIP_CALL( SCIPaddCons(scip, newMixWeymouth) );
      SCIP_CALL( SCIPreleaseCons(scip, &newMixWeymouth) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exprFIN) );

#ifdef USE_ADIABATIC_CAPACITY
      SCIP_CALL( SCIPreleaseExpr(scip, &expGMRTz) );
#endif

#ifdef VARIANT_WEYMOUTH
      {
         SCIP_CONS* lowerLambda;
         SCIP_EXPR* arcLambda;
         SCIP_EXPR* diff;
         SCIP_EXPR* exprZplus;
         SCIP_EXPR* exprZminus;
         SCIP_Real unmixedLambda = 0.65 * computeWeymouthConstant(probdata->network, pipe->length, pipe->diameter, pipe->roughness, meanpressure, probdata->papay);

         assert( arc->Lambda != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &arcLambda, arc->Lambda, NULL, NULL) );

         SCIP_CALL( SCIPcreateExprVar(scip, &exprZplus, arc->positiveFlowBinvar, NULL, NULL) );
         assert( arc->negativeFlowBinvar != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &exprZminus, arc->negativeFlowBinvar, NULL, NULL) );

         expr[0] = exprZplus;
         expr[1] = exprZminus;
         SCIP_CALL( SCIPcreateExprSum(scip, &diff, 2, expr, NULL, 0.0, NULL, NULL) );
         expr[0] = diff;
         expr[1] = arcLambda;
         SCIP_CALL( SCIPcreateExprProduct(scip, &exprFIN, 2, expr, 1.0, NULL, NULL) );

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "LambdaMin%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &lowerLambda, name, exprFIN, unmixedLambda, 100.0) );
         SCIP_CALL( SCIPaddCons(scip, lowerLambda) );
         SCIP_CALL( SCIPprintCons(scip, lowerLambda, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &lowerLambda) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprFIN) );

         SCIP_CALL( SCIPreleaseExpr(scip, &arcLambda) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprZplus) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprZminus) );
         SCIP_CALL( SCIPreleaseExpr(scip, &diff) );
      }
#endif

#ifdef USE_ADIABATIC_CAPACITY
      SCIP_CALL( SCIPreleaseExpr(scip, &exprGamma) );
#endif
#endif  /* SCIP_VERSION >= 800 */
   }
   else
   {
      LAMBDA = computeWeymouthConstant(probdata->network, pipe->length, pipe->diameter, pipe->roughness, meanpressure, probdata->papay);
      SCIPdebugMessage("Lambda of pipe <%s>: <%f>\n", arc->id, LAMBDA);

      /* use cons_abspower.h : 1/Lambda *(PI_diff) = flowvars * |flowvars| */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "WeymouthEq#%s#", arc->id);

      /* Create constraint and add it to scip */
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
      SCIP_CALL( SCIPcreateConsBasicSignpowerNonlinear(scip, &WeymouthEquality_Cons, consname, probdata->flowvars[arcposition], pipe->PIdiffVar,
            2.0, 0.0, - (1.0 / LAMBDA), 0.0, 0.0) );
#else
      SCIP_CALL( SCIPcreateConsBasicAbspower(scip, &WeymouthEquality_Cons, consname, probdata->flowvars[arcposition], pipe->PIdiffVar,
            2.0, 0.0, - (1.0 / LAMBDA), 0.0, 0.0) );
#endif
      SCIP_CALL( SCIPaddCons(scip, WeymouthEquality_Cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &WeymouthEquality_Cons) );
   }

   return SCIP_OKAY;
}

/** generates the constraints for the algebraic gas flow model (ISO4 / Weymouth) */
static
SCIP_RETCODE generateAlgebraicModelNodeMu(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   arcposition         /**< position in arc array */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Arc* arc;
   GAS_Pipe* pipe;
   GAS_Node* sourcenode;
   GAS_Node* targetnode;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];
   SCIP_CONS* PI_diff_cons;
   SCIP_CONS* WeymouthEquality_Cons;
   char consname[SCIP_MAXSTRLEN];
   SCIP_Real LAMBDA;
   SCIP_Real meanpressure;

   assert( scip != NULL );

   probdata =  SCIPgetProbData(scip);

   arc = &(probdata->network->arcs_ptr[arcposition]);
   pipe = arc->detailed_info;
   sourcenode = arc->sourcenode;
   targetnode = arc->targetnode;

   assert( sourcenode->PIvar != NULL );
   assert( targetnode->PIvar != NULL );

   /*
    *  generate the coupling between the squared pressure variables and the pressure diff variable
    *  PI_diff - (PI_source - PI_target) = 0.0
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "PIdiffCons#%s#", arc->id);

   vars[0] = pipe->PIdiffVar;
   vars[1] = sourcenode->PIvar;
   vars[2] = targetnode->PIvar;

   vals[0] =  1.0;
   vals[1] = -1.0;
   vals[2] =  1.0;

   /* Create linear constraint and add it to SCIP */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &PI_diff_cons, consname, 3, vars, vals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, PI_diff_cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &PI_diff_cons) );

   /*
    *  Generate the Weymouth equation
    */

   meanpressure = ComputeMeanPressure(sourcenode->pressureMin, sourcenode->pressureMax, targetnode->pressureMin, targetnode->pressureMax);

   if ( ( probdata->mixing || probdata->linearMixing ) && ! probdata->approx )
   {
#if SCIP_VERSION >= 800
      char name[SCIP_MAXSTRLEN];
      SCIP_CONS* newMixWeymouth;
      SCIP_EXPR* exprMol;
      SCIP_EXPR* exprFrac;
      SCIP_EXPR* exprFIN;
      SCIP_EXPR* expr[3];
      SCIP_Real  coef[2];
      SCIP_EXPR* exprQsqrd;
      SCIP_EXPR* exprQ;
      SCIP_EXPR* exprPIdiff;
      SCIP_EXPR* exprZplus;
      SCIP_EXPR* exprZminus;
      SCIP_EXPR* MU;
      SCIP_EXPR* LAMBDAsource;
      SCIP_EXPR* LAMBDAsink;
      SCIP_EXPR* diff;
      SCIP_EXPR* diffQsqrd;

      /* set coefficients */
      SCIP_Real LambdaMix = computeWeymouthConstantMixing(probdata->network, pipe->length, pipe->diameter, pipe->roughness, meanpressure, probdata->papay);
      SCIP_Real molarMass1 = probdata->molarMass1;
      SCIP_Real molarMass2 = probdata->molarMass2;
      SCIP_Real coefM = molarMass1 - molarMass2;

      /* create expr variables */
      assert( arc->positiveFlowBinvar != NULL );
      SCIP_CALL( SCIPcreateExprVar(scip, &exprZplus, arc->positiveFlowBinvar, NULL, NULL) );
      assert( arc->negativeFlowBinvar != NULL );
      SCIP_CALL( SCIPcreateExprVar(scip, &exprZminus, arc->negativeFlowBinvar, NULL, NULL) );

      if ( ! probdata->linearMixing )
      {
         /* Lambda at source */
         assert( arc->sourcenode->nodeMixingRatioMol != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &MU, arc->sourcenode->nodeMixingRatioMol, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &MU, &coefM, molarMass2, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprPow(scip, &exprFrac, exprMol, -1.0, NULL, NULL) );
         expr[0] = exprFrac;
         expr[1] = exprZplus;
         SCIP_CALL( SCIPcreateExprProduct(scip, &LAMBDAsource, 2, expr, LambdaMix, NULL, NULL) );
         /* free expr */
         SCIP_CALL( SCIPreleaseExpr(scip, &MU) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprMol) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprFrac) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprZplus) );

         /* Lambda at sink */
         assert( arc->targetnode->nodeMixingRatioMol != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &MU, arc->targetnode->nodeMixingRatioMol, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &MU, &coefM, molarMass2, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprPow(scip, &exprFrac, exprMol, -1.0, NULL, NULL) );
         expr[0] = exprFrac;
         expr[1] = exprZminus;
         SCIP_CALL( SCIPcreateExprProduct(scip, &LAMBDAsink, 2, expr, LambdaMix, NULL, NULL) );
         /* free expr */
         SCIP_CALL( SCIPreleaseExpr(scip, &MU) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprMol) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprFrac) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprZminus) );
      }
      else
      {
         /* Lambda at source */
         /* (1/M1 -1/M2) * muMass + 1/M2 */
         assert( arc->sourcenode->nodeMixingRatio != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &MU, arc->sourcenode->nodeMixingRatio, NULL, NULL) );
         coefM = 1.0 / molarMass1 - 1.0 / molarMass2;
         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &MU, &coefM, 1.0 / molarMass2, NULL, NULL) );
         SCIP_CALL( SCIPreleaseExpr(scip, &MU) );
         expr[0] = exprMol;
         expr[1] = exprZplus;

         /* LambdaMix* z^+ * ( (1/M1 -1/M2) * muMass + 1/M2 ) */
         SCIP_CALL( SCIPcreateExprProduct(scip, &LAMBDAsource, 2, expr, LambdaMix, NULL, NULL) );

         /* free expr */
         SCIP_CALL( SCIPreleaseExpr(scip, &exprZplus) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprMol) );

         /* Lambda at sink */
         assert( arc->targetnode->nodeMixingRatio != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &MU, arc->targetnode->nodeMixingRatio, NULL, NULL) );
         coefM = 1.0 / molarMass1 - 1.0 / molarMass2;
         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &MU, &coefM, 1.0 / molarMass2, NULL, NULL) );
         SCIP_CALL( SCIPreleaseExpr(scip, &MU) );
         expr[0] = exprMol;
         expr[1] = exprZminus;
         SCIP_CALL( SCIPcreateExprProduct(scip, &LAMBDAsink, 2, expr, LambdaMix, NULL, NULL) );

         /* free expr */
         SCIP_CALL( SCIPreleaseExpr(scip, &exprMol) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprZminus) );
      }

#ifdef VARIANT_WEYMOUTH_NODE1
      /* create expr variables */
      assert( probdata->flowvars[arcposition] != NULL );
      SCIP_CALL( SCIPcreateExprVar(scip, &exprQ, probdata->flowvars[arcposition], NULL, NULL) );
      SCIP_CALL( SCIPcreateExprPow(scip, &exprQsqrd, exprQ, 2.0, NULL, NULL) );
      assert( pipe->PIdiffVar != NULL );
      SCIP_CALL( SCIPcreateExprVar(scip, &exprPIdiff, pipe->PIdiffVar, NULL, NULL) );

      /* difference of the big Lambdas */
      coef[0] = 1.0;
      coef[1] = -1.0;
      expr[0] = LAMBDAsource;
      expr[1] = LAMBDAsink;
      SCIP_CALL( SCIPcreateExprSum(scip, &diff, 2, expr, coef, 0.0, NULL, NULL) );

      /* difference muliplied with q^2 in order to get the rhs of the Weymouth eqn */
      expr[0] = diff;
      expr[1] = exprQsqrd;
      SCIP_CALL( SCIPcreateExprProduct(scip, &diffQsqrd, 2, expr, 1.0, NULL, NULL) );

      /* Weymouth eqn solved for zero */
      expr[0] = exprPIdiff;
      expr[1] = diffQsqrd;
      SCIP_CALL( SCIPcreateExprSum(scip, &exprFIN, 2, expr, coef, 0.0, NULL, NULL) );

      /* EDIT: faster progress in the dual objective in Gaslib135 with this eqn */
      /* This type of equation is more problematic for SCIP it needs 4 times more time to solve the same problems with this type of constraint than the below one. */

      (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "MixWeymouth%s", arc->id);
      SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &newMixWeymouth, name, exprFIN, 0.0, 0.0) );
      SCIP_CALL( SCIPaddCons(scip, newMixWeymouth) );
      SCIP_CALL( SCIPreleaseCons(scip, &newMixWeymouth) );
#endif

      /* Generate LAMBDA via its own constraint and then the WeyEqn */
      /* This is fastest by now! It seems like the signpower expr which we use to genereate the weymout eqn is really good for SCIP. */
      {
         SCIP_CONS* LambdaCONS;
         SCIP_EXPR* arcLambda;

         /* sum of the big Lambdas */
         expr[0] = LAMBDAsource;
         expr[1] = LAMBDAsink;
         SCIP_CALL( SCIPcreateExprSum(scip, &diff, 2, expr, NULL, 0.0, NULL, NULL) );

         assert( arc->Lambda != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &arcLambda, arc->Lambda, NULL, NULL) );

         /* LambdaSource + LambdaSink - arcLambda =0*/
         expr[0] = arcLambda;
         expr[1] = diff;
         coef[0] = 1;
         coef[1] = -1;
         SCIP_CALL( SCIPcreateExprSum(scip, &exprFIN, 2, expr, coef, 0.0, NULL, NULL) );

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "Lambda%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &LambdaCONS, name, exprFIN, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, LambdaCONS) );
         SCIP_CALL( SCIPreleaseCons(scip, &LambdaCONS) );

         /* free expr */
         SCIP_CALL( SCIPreleaseExpr(scip, &LAMBDAsource) );
         SCIP_CALL( SCIPreleaseExpr(scip, &LAMBDAsink) );
         SCIP_CALL( SCIPreleaseExpr(scip, &arcLambda) );
         SCIP_CALL( SCIPreleaseExpr(scip, &diff) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprFIN) );

         /* WEYMOUTH EQN */

         /* generate expr needed again since GASexprFree did free them above*/
         /* no assert needed since the assert is already applied when first generating the epxr */
         SCIP_CALL( SCIPcreateExprVar(scip, &arcLambda, arc->Lambda, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprVar(scip, &exprQ, probdata->flowvars[arcposition], NULL, NULL) );
         SCIP_CALL( SCIPcreateExprVar(scip, &exprPIdiff, pipe->PIdiffVar, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSignpower(scip, &exprQsqrd, exprQ, 2.0, NULL, NULL) );

         /* now the WeyEqn is generated */
         expr[0] = arcLambda;
         expr[1] = exprQsqrd;
         SCIP_CALL( SCIPcreateExprProduct(scip, &diffQsqrd, 2, expr, 1.0, NULL, NULL) );

         /* Weymouth eqn solved for zero */
         expr[0] = exprPIdiff;
         expr[1] = diffQsqrd;
         SCIP_CALL( SCIPcreateExprSum(scip, &exprFIN, 2, expr, coef, 0.0, NULL, NULL) );

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "MixWeymouth%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &newMixWeymouth, name, exprFIN, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, newMixWeymouth) );
         SCIP_CALL( SCIPreleaseCons(scip, &newMixWeymouth) );

         SCIP_CALL( SCIPreleaseExpr(scip, &exprQ) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprPIdiff) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprQsqrd) );
         SCIP_CALL( SCIPreleaseExpr(scip, &arcLambda) );
         SCIP_CALL( SCIPreleaseExpr(scip, &diffQsqrd) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprFIN) );
      }

#ifdef VARIANT_WEYMOUTH_NODE2
      /* makes scip run slower by a factor of 3... same in the arc mu case. */
      {
         SCIP_CONS* lowerLambda;
         SCIP_EXPR* arcLambda;
         SCIP_Real unmixedLambda = 0.65* computeWeymouthConstant( probdata->network, pipe->length, pipe->diameter, pipe->roughness, meanpressure, probdata->papay);

         assert( arc->Lambda != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &arcLambda, arc->Lambda, NULL, NULL) );

         expr[0] = exprZplus;
         expr[1] = exprZminus;
         SCIP_CALL( SCIPcreateExprSum(scip, &diff, 2, expr, NULL, 0.0, NULL, NULL) );
         expr[0] = diff;
         expr[1] = arcLambda;
         SCIP_CALL( SCIPcreateExprProduct(scip, &exprFIN, 2, expr, 1.0, NULL, NULL) );

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "LambdaMin%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &lowerLambda, name, exprFIN, unmixedLambda, 100.0) );
         SCIP_CALL( SCIPaddCons(scip, lowerLambda) );
         SCIP_CALL( SCIPprintCons(scip, lowerLambda, NULL) );
         SCIP_CALL( SCIPreleaseCons(scip, &lowerLambda) );

         SCIP_CALL( SCIPreleaseExpr(scip, &arcLambda) );
      }
#endif

      /* second fastest by now */
      {
#ifdef VARIANT_WEYMOUTH_NODE3
         SCIP_CONS* posWeymouth;
         SCIP_CONS* negWeymouth;

         /* Lambda source */
         assert( arc->sourcenode->nodeMixingRatioMol != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &MU, arc->sourcenode->nodeMixingRatioMol, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &MU, &coefM, molarMass2, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprPow(scip, &exprFrac, exprMol, -1.0, NULL, NULL) );
         expr[0] = exprFrac;
         expr[1] = exprQsqrd;
         SCIP_CALL( SCIPcreateExprProduct(scip, &LAMBDAsource, 2, expr, LambdaMix, NULL, NULL) );

         expr[0] = exprPIdiff;
         expr[1] = LAMBDAsource;
         SCIP_CALL( SCIPcreateExprSum(scip, &diff, 2, expr, coef, 0.0, NULL, NULL) );

         expr[0] = diff;
         expr[1] = exprZplus;
         SCIP_CALL( SCIPcreateExprProduct(scip, &exprFIN, 2, expr, 1.0, NULL, NULL) );

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "posWeymouth%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &posWeymouth, name, exprFIN, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, posWeymouth) );
         SCIP_CALL( SCIPreleaseCons(scip, &posWeymouth) );

         /* Lambda at sink */
         assert( arc->targetnode->nodeMixingRatioMol != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &MU, arc->targetnode->nodeMixingRatioMol, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &MU, &coefM, molarMass2, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprPow(scip, &exprFrac, exprMol, -1.0, NULL, NULL) );
         expr[0] = exprFrac;
         expr[1] = exprQsqrd;
         SCIP_CALL( SCIPcreateExprProduct(scip, &LAMBDAsink, 2, expr, LambdaMix, NULL, NULL) );

         expr[0] = exprPIdiff;
         expr[1] = LAMBDAsink;
         SCIP_CALL( SCIPcreateExprSum(scip, &diff, 2, expr, NULL, 0.0, NULL, NULL) ); /* added in this case i.e. coef =NULL = all are one */

         expr[0] = diff;
         expr[1] = exprZminus;
         SCIP_CALL( SCIPcreateExprProduct(scip, &exprFIN, 2, expr, 1.0, NULL, NULL) );

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "negWeymouth%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &negWeymouth, name, exprFIN, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, negWeymouth) );
         SCIP_CALL( SCIPreleaseCons(scip, &negWeymouth) );
#endif
      }

#ifdef VARIANT_WEYMOUTH_NODE4
      /* the following is worse */
      {
         SCIP_CONS* posWeymouth;
         SCIP_CONS* negWeymouth;

         /* Lambda at source */
         assert( arc->sourcenode->nodeMixingRatioMol != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &MU, arc->sourcenode->nodeMixingRatioMol, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &MU, &coefM, molarMass2, NULL, NULL) );
         expr[0] = exprMol;
         expr[1] = exprPIdiff;
         SCIP_CALL( SCIPcreateExprProduct(scip, &LAMBDAsource, 2, expr, (1.0 / LambdaMix), NULL, NULL) );

         expr[0] = LAMBDAsource;
         expr[1] = exprQsqrd;
         SCIP_CALL( SCIPcreateExprSum(scip, &diff, 2, expr, coef, 0.0, NULL, NULL) );

         expr[0] = diff;
         expr[1] = exprZplus;
         SCIP_CALL( SCIPcreateExprProduct(scip, &exprFIN, 2, expr, 1.0, NULL, NULL) );

         (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "posWeymouth%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &posWeymouth, name, exprFIN, 0.0, 0.0) );
         SCIP_CALL( SCIPprintCons(scip, posWeymouth, NULL) );
         SCIP_CALL( SCIPaddCons(scip, posWeymouth) );
         SCIP_CALL( SCIPreleaseCons(scip, &posWeymouth) );

         /* Lambda at sink */
         assert( arc->sourcenode->nodeMixingRatioMol != NULL );
         SCIP_CALL( SCIPcreateExprVar(scip, &MU, arc->targetnode->nodeMixingRatioMol, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSum(scip, &exprMol, 1, &MU, &coefM, molarMass2, NULL, NULL) );
         expr[0] = exprMol;
         expr[1] = exprPIdiff;
         SCIP_CALL( SCIPcreateExprProduct(scip, &LAMBDAsink, 2, expr, (1.0 / LambdaMix), NULL, NULL) );

         expr[0] = LAMBDAsink;
         expr[1] = exprQsqrd;
         SCIP_CALL( SCIPcreateExprSum(scip, &diff, 2, expr, NULL, 0.0, NULL, NULL) ); /*adding*/

         expr[0] = diff;
         expr[1] = exprZminus;
         SCIP_CALL( SCIPcreateExprProduct(scip, &exprFIN, 2, expr, 1.0, NULL, NULL) );

          (void)SCIPsnprintf(name, SCIP_MAXSTRLEN, "negWeymouth%s", arc->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &negWeymouth, name, exprFIN, 0.0, 0.0) );
         SCIP_CALL( SCIPprintCons(scip, negWeymouth, NULL) );
         SCIP_CALL( SCIPaddCons(scip, negWeymouth) );
         SCIP_CALL( SCIPreleaseCons(scip, &negWeymouth) );
      }
#endif
#endif /* SCIP_VERSION >= 800 */
   }
   else
   {
      LAMBDA = computeWeymouthConstant(probdata->network, pipe->length, pipe->diameter, pipe->roughness, meanpressure, probdata->papay);
      SCIPdebugMessage("Lambda of pipe <%s>: <%f>\n", arc->id, LAMBDA);

      /* Use cons_abspower.h : 1/Lambda *(PI_diff) = flowvars * |flowvars| */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "WeymouthEq#%s#", arc->id);

      /* Create constraint and add it to scip */
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
      SCIP_CALL( SCIPcreateConsBasicSignpowerNonlinear(scip, &WeymouthEquality_Cons, consname, probdata->flowvars[arcposition], pipe->PIdiffVar, 2.0, 0.0, -(1/LAMBDA), 0.0, 0.0) );
#else
      SCIP_CALL( SCIPcreateConsBasicAbspower(scip, &WeymouthEquality_Cons, consname, probdata->flowvars[arcposition], pipe->PIdiffVar, 2.0, 0.0, -(1/LAMBDA), 0.0, 0.0) );
#endif
      SCIP_CALL( SCIPaddCons(scip, WeymouthEquality_Cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &WeymouthEquality_Cons) );
   }

   return SCIP_OKAY;
}


/** Adjust mixing ratio lower bounds based on gas properties from scenario file */
static
SCIP_RETCODE adjustMixigRatioLowerBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   GAS_Network* Network = probdata->network;
   SCIP_Real muLow = 1.0;
   SCIP_Real epsilon = 1e-5;

   /* the lower bound is reduced by epsilon since for some reason GL40 with mixing is infeasible without this */

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   for (int j = 0; j < Network->numnodes; j++)
   {
      GAS_Node* Node = &(Network->nodes_ptr[j]);
      assert( Node != NULL );

      if ( Node->type == ENTRY && Node->flow > 0 )
      {
         /* calculate the mass% mu */
         SCIP_Real mu = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);

         if ( mu < muLow )
            muLow = mu - epsilon;
      }
   }

   if (muLow < 0.0)
      muLow = 0.0;

   assert( muLow < 1.0 );
   assert( muLow >= 0.0 );

   /* loop through nodes and just lower bound */
   for (int j = 0; j < Network->numnodes; j++)
   {
      GAS_Node* Node = &(Network->nodes_ptr[j]);
      assert( Node != NULL );

      SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, muLow) );
   }

   return SCIP_OKAY;
}

/** generate mixing ratio caclulation "constraint" for each node */
static
SCIP_RETCODE generateMixingRatioNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
#if SCIP_VERSION >= 800
   char consname[SCIP_MAXSTRLEN];
   GAS_Network* Network = probdata->network;
   SCIP_CONS* mixingCons;
   SCIP_CONS* mixingConsMol;
   int j;
   int c = 0;
   int sum = 0;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   for (j = 0; j < Network->numnodes; j++)
   {
      GAS_Node* Node = &(Network->nodes_ptr[j]);
      GAS_Arc* arc;
      int numInarcs = 0;
      int numOutarcs = 0;

      assert( Node != NULL );

      /* Determine the number of in and outarcs */
      arc = Node->inarcs;
      while ( arc != NULL )
      {
         numInarcs++;
         arc = arc->next_inarc;
      }

      arc = Node->outarcs;
      while ( arc != NULL )
      {
         numOutarcs++;
         arc = arc->next_outarc;
      }

      sum = numInarcs + numOutarcs;

      if ( sum == 1 && Node->type == ENTRY )
      {

         assert( Node->nodeMixingRatio != NULL );

         if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
         {
            SCIPdebugMsg(scip, "Heavy Gas, i.e., fix mu = 1.\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 1.0) );

            if ( ! probdata->linearMixing )
            {
               assert( Node->nodeMixingRatioMol != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatioMol, 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatioMol, 1.0) );
            }
         }
         else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
         {
            SCIPdebugMsg(scip, "Light Gas i.e. fix mu = 0\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 0.0) );

            if ( ! probdata->linearMixing )
            {
               assert( Node->nodeMixingRatioMol != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatioMol, 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatioMol, 0.0) );
            }
         }
         /* otherwise the gas is a mixture of the gaes with the highest and lowest molar mass and we calculate mu */
         else
         {
            SCIP_Real mu;
            SCIP_Real muMol;
            SCIP_Real M1 = 1.0 / probdata->molarMass1;
            SCIP_Real M2 = 1.0 / probdata->molarMass2;

            /* calculate the mass% mu */
            mu = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);

            SCIPdebugMsg(scip, "Mixed Gas, molar mass = %f, mu = %f", Node->molarMass, mu);
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, mu) );

            if ( ! probdata->linearMixing )
            {
               /* calculate mol% mu */
               muMol = (mu * M1) / (mu * (M1 - M2) + M2);
               SCIPdebugMsg(scip, "mu mol = %f\n", muMol);

               SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, muMol) );
               SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, muMol) );
            }
         }
      }

      if ( sum == 1 && Node->type == EXIT )
      {
         SCIP_CONS* fixNodeMuCons;
         SCIP_VAR* vars[2];
         SCIP_Real vals[2];

         vals[0] = 1.0;
         vals[1] = -1.0;

         if ( Node->inarcs != NULL )
         {
            arc = Node->inarcs;
            vars[0] = arc->sourcenode->nodeMixingRatio;
            assert( arc->sourcenode->nodeMixingRatio != NULL );
            SCIPdebugMsg(scip, "Arc is inarc\n");
         }
         else
         {
            arc = Node->outarcs;
            vars[0] = arc->targetnode->nodeMixingRatio;
            assert( arc->targetnode->nodeMixingRatio != NULL );
            SCIPdebugMsg(scip, "Arc is outarc\n");
         }

         assert( Node->nodeMixingRatio != NULL );
         vars[1] = Node->nodeMixingRatio;

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-3-EXIT-fixNodeMixRatio%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixNodeMuCons, consname, 2, vars, vals, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, fixNodeMuCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &fixNodeMuCons) );
      }
      /* node has 2 or more neighbours and mu cannot be fixed */
      else if ( sum >= 2 )
      {
         SCIP_EXPR* exprQZprod[2];
         SCIP_EXPR** exprQZprodResult = NULL;
         SCIP_EXPR* exprQZMuprod[2];
         SCIP_EXPR** exprQZMuprodResult = NULL;
         SCIP_EXPR* exprSUMin = NULL;
         SCIP_EXPR* exprSUMinMu = NULL;
         SCIP_EXPR* exprQaOut = NULL;
         SCIP_EXPR* exprZaOut = NULL;
         SCIP_EXPR* exprMUaOut = NULL;
         SCIP_EXPR** exprQZprodResultO = NULL;
         SCIP_EXPR** exprQZMuprodResultO = NULL;
         SCIP_EXPR* exprSUMout = NULL;
         SCIP_EXPR* exprSUMoutMu = NULL;
         SCIP_EXPR* exprZsum = NULL;
         SCIP_EXPR* exprZarr[2];
         SCIP_Real coefPlusOneMinusOne[2];
         SCIP_EXPR* exprN[2];
         SCIP_EXPR* exprMuNode = NULL;
         SCIP_EXPR* exprProdN[2];
         SCIP_EXPR* exprMuSumProd = NULL;
         SCIP_EXPR* exprFsumm[2];
         SCIP_EXPR* exprFinMu;
         SCIP_EXPR* expFin;
         SCIP_EXPR* expMuMolNode;

         /* used in the loops below as offset i.e if the node is an entry the inflow is added to the flow sums */
         SCIP_Real inflow = 0.0;
         SCIP_Real inflowMU = 0.0;
         SCIP_Real copyInflow = 0.0;
         SCIP_Real copyInflowMU = 0.0;

         coefPlusOneMinusOne[0] = 1.0;
         coefPlusOneMinusOne[1] = -1.0;

         /* Adding the inflow at an as entry offset is complicated and happens as follows: if either ONLY inarcs or
          * outarcs exist, the offset is added in one of the next two if loops.*/
         /* Otherwise the offset is added in the 3rd if loop (i.e. if numInarcs && numoutarcs >0) */
         /* The complication is caused by the fact that the flow value has to be added exactly once. */
         /* The flow is always added in the denominator but only in the numerator if mu=1. */
         if ( Node->type == ENTRY )
         {
            if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
            {
               /* heavy gas => mu=1 */
               inflow = Node->flow;
               inflowMU = 1;
            }
            else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
            {
               /* light gas => mu=0 */
               inflow = Node->flow;
               inflowMU = 0;
            }
            /* gas is a mixture and we first need to calculate mu */
            else
            {
               /* calculate the mass% mu */
               inflowMU = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);
               inflow = Node->flow;
               SCIPdebugMsg(scip, "inflowMU = %f, molar mass = %f, node = %s\n", inflowMU, Node->molarMass, Node->id);
            }
         }

         if ( numInarcs > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZprodResult, numInarcs) );
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResult, numInarcs) );

            /* The loop below builds the pairwise producs q_a*z_a and q_a*z_a*mu_a for the inarcs. */
            /* NOTE: MASS% mixing Ratio is used i.e. mixingRatioMass */

            /* reset counter */
            c = 0;

            /* fill in incoming arcs */
            arc = Node->inarcs;
            while ( arc != NULL )
            {
               SCIP_EXPR* exprQaIn = NULL;
               SCIP_EXPR* exprZaIn = NULL;
               SCIP_EXPR* exprMUaIn = NULL;

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaIn, arc->flowvar, NULL, NULL) );
               assert( arc->positiveFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaIn, arc->positiveFlowBinvar, NULL, NULL) );

               exprQZprod[0] = exprQaIn;
               exprQZprod[1] = exprZaIn;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZprodResult[c], 2, exprQZprod, 1.0, NULL, NULL) ); /* Pairwise q_a*z_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaIn) );

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaIn, arc->flowvar, NULL, NULL) );
               assert( arc->positiveFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaIn, arc->positiveFlowBinvar, NULL, NULL) );
               assert( arc->sourcenode->nodeMixingRatio != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaIn, arc->sourcenode->nodeMixingRatio, NULL, NULL) );

               exprQZMuprod[0] = exprQZprodResult[c];
               exprQZMuprod[1] = exprMUaIn;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResult[c], 2, exprQZMuprod, 1.0, NULL, NULL) ); /* Pairwise q_a*z_a*mu_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaIn) );

               arc = arc->next_inarc;
               c++;
            }

            assert( c == numInarcs );
            c = 0;

            /* need to add the inflow value as offset - see details above */
            if ( numOutarcs == 0 )
            {
               copyInflow = inflow;
               copyInflowMU = inflowMU;
            }

            /* now the sums over the above created producs are build */

            /* Sum over inarcs q_a*z_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMin, numInarcs, exprQZprodResult, NULL, copyInflow, NULL, NULL) );

            /* Sum over inarcs q_a*z_a*mu_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMinMu, numInarcs, exprQZMuprodResult, NULL, copyInflowMU * copyInflow, NULL, NULL) );

            /*free arrays*/
            for (int i = 0; i < numInarcs; i++)
            {
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZprodResult[i]) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZMuprodResult[i]) );
            }

            SCIPfreeBufferArray(scip, &exprQZMuprodResult );
            SCIPfreeBufferArray(scip, &exprQZprodResult);
         }

         if ( numOutarcs > 0 )
         {
            /* Analog calculation now for the out arcs. */
            /* The loop below builds the pairwise producs q_a*z_a and q_a*z_a*mu_a for the outarcs. */
            /* NOTE: MASS% mixing Ratio is used i.e. mixingRatioMass. */
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZprodResultO, numOutarcs) );
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResultO, numOutarcs) );

            c = 0;
            arc = Node->outarcs;

            /* fill in outgoing arcs */
            while ( arc != NULL )
            {
               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaOut, arc->flowvar, NULL, NULL) );
               assert( arc->negativeFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaOut, arc->negativeFlowBinvar, NULL, NULL) );
               exprQZprod[0] = exprQaOut;
               exprQZprod[1] = exprZaOut;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZprodResultO[c], 2, exprQZprod, 1.0, NULL, NULL) ); /* pairwise q_a * z_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaOut) );

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaOut, arc->flowvar, NULL, NULL) );
               assert( arc->negativeFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaOut, arc->negativeFlowBinvar, NULL, NULL) );
               assert( arc->targetnode->nodeMixingRatio != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaOut, arc->targetnode->nodeMixingRatio, NULL, NULL) );
               exprQZMuprod[0] = exprQZprodResultO[c];
               exprQZMuprod[1] = exprMUaOut;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResultO[c], 2, exprQZMuprod, 1.0, NULL, NULL) ); /* pairwise q_a * z_a * mu_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaOut) );

               arc = arc->next_outarc;
               c++;
            }

            assert( c == numOutarcs );

            /* In the case that the node has only inflow on outarcs i.e. the flow is negative we multiply it by -1. */
            SCIP_Real minusOne[numOutarcs];

            if ( numInarcs == 0 )
            {
               copyInflow = inflow;
               copyInflowMU = inflowMU;

               for (int i = 0; i < numOutarcs; i++)
                  minusOne[i] = -1.0;
            }
            else
            {
               for (int i = 0; i < numOutarcs; i++)
                  minusOne[i] = 1.0;
            }

            /* now the sums over the above created producs are built */

            /* Sum over inarcs q_a*z_a prod + inflow if needed. see above*/
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMout, numOutarcs, exprQZprodResultO, minusOne, copyInflow, NULL, NULL) );

            /* Sum over inarcs q_a*z_a*mu_a prod + inflow if needed. see above*/
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMoutMu, numOutarcs, exprQZMuprodResultO, minusOne, copyInflowMU * copyInflow, NULL, NULL) );

            /*free arrays*/
            for (int i = 0; i < numOutarcs; i++)
            {
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZprodResultO[i]) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZMuprodResultO[i]) );
            }

            SCIPfreeBufferArray(scip, &exprQZMuprodResultO);
            SCIPfreeBufferArray(scip, &exprQZprodResultO);
         }

         /* In the following the mixing equation - Equation (7) of the mixing paper- is build */
         if ( numInarcs > 0 && numOutarcs > 0 )
         {
            SCIP_EXPR* exprNsum = NULL;

            /* Now we first build the numerator for the mixing eqaution and then the denominator down below */
            /* Note that the equation will be written in non fractional form (no div by 0) */

            /* Numerator sum of the mixing equation - EQ (7) in the paper */
            exprZarr[0] = exprSUMinMu;
            exprZarr[1] = exprSUMoutMu;
            SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 2, exprZarr, coefPlusOneMinusOne, inflowMU * inflow, NULL, NULL) );

            /* Mu_Node * NominatorSum - Denominator sum muliplied by mu */
            exprN[0] = exprSUMin;
            exprN[1] = exprSUMout;

            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );

            SCIP_CALL( SCIPcreateExprSum(scip, &exprNsum, 2, exprN, coefPlusOneMinusOne, inflow, NULL, NULL) ); /* inflow at node */
            exprProdN[0] = exprMuNode;
            exprProdN[1] = exprNsum;
            SCIP_CALL( SCIPcreateExprProduct(scip, &exprMuSumProd, 2, exprProdN, 1.0, NULL, NULL) ); /* mu * inflow or C1 * M1 */

            exprFsumm[0] = exprMuSumProd;
            exprFsumm[1] = exprZsum;

            SCIP_CALL( SCIPreleaseExpr(scip, &exprNsum) );
         }
         else if ( numInarcs > 0 && numOutarcs == 0)
         {
            /* Denominator is only sum over inarcs = exprSUMinMu */
            /* Mu_Node * NominatorSum */
            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );
            exprProdN[0] = exprMuNode;
            exprProdN[1] = exprSUMin;
            SCIP_CALL( SCIPcreateExprProduct(scip, &exprMuSumProd, 2, exprProdN, 1.0, NULL, NULL) ); /* mu * inflow or C1 * M1 */

            exprFsumm[0] = exprMuSumProd;
            exprFsumm[1] = exprSUMinMu;
         }
         else if ( numInarcs == 0 && numOutarcs > 0 )
         {
            /* Denominator is only sum over outarcs = exprSUMoutMu */
            /* Mu_Node * NominatorSum. Note that the inflow here is negative which is why it is multiplied with -1 above !*/
            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );
            exprProdN[0] = exprMuNode;
            exprProdN[1] = exprSUMout;
            SCIP_CALL( SCIPcreateExprProduct(scip, &exprMuSumProd, 2, exprProdN, 1.0, NULL, NULL) ); /* mu * inflow or C1 * M1 */

            exprFsumm[0] = exprMuSumProd;
            exprFsumm[1] = exprSUMoutMu;
         }

         /* final equation */
         /* Mu_Node * NominatorSum - DenominatorSum = 0 */
         SCIP_CALL( SCIPcreateExprSum(scip, &exprFinMu, 2, exprFsumm, coefPlusOneMinusOne, 0.0, NULL, NULL) );

         /* Constraint */
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-NodeMixingRatioCalculation_%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprFinMu, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, mixingCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );

         /* free all expr */
         if ( numInarcs > 0 && numOutarcs > 0 )
         {
            /*release expr*/
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMinMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMoutMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMin) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMout) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprZsum) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuSumProd) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprFinMu) );
         }
         else if ( numInarcs > 0 && numOutarcs == 0 )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMinMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMin) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuSumProd) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprFinMu) );
         }
         else if ( numInarcs == 0 && numOutarcs > 0 )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMout) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMoutMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuSumProd) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprFinMu) );
         }

         /* here the mole% mixr. is calculated */
         if ( ! probdata->linearMixing )
         {
            SCIP_EXPR* expr[2];
            SCIP_Real coef[1];
            SCIP_EXPR* C2inv;
            SCIP_EXPR* C1fracC2;
            SCIP_Real M2 = 1.0 / probdata->molarMass2;
            SCIP_EXPR* exprC2;

            /* need to genereate the expr again since it is freeded above.*/
            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );

            coef[0] = (1.0 / probdata->molarMass1) - (1.0 / probdata->molarMass2); /* this coef is negative */
            expr[0] = exprMuNode;

            SCIP_CALL( SCIPcreateExprSum(scip, &exprC2, 1, expr, coef, M2, NULL, NULL) );
            SCIP_CALL( SCIPcreateExprPow(scip, &C2inv, exprC2, -1.0, NULL, NULL) );

            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );

            expr[0] = exprMuNode;
            expr[1] = C2inv;
            SCIP_CALL( SCIPcreateExprProduct(scip, &C1fracC2, 2, expr, 1.0 / probdata->molarMass1, NULL, NULL) );

            assert( Node->nodeMixingRatioMol != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &expMuMolNode, Node->nodeMixingRatioMol, NULL, NULL) );

            expr[0] = expMuMolNode;
            expr[1] = C1fracC2;
            SCIP_CALL( SCIPcreateExprSum(scip, &expFin, 2, expr, coefPlusOneMinusOne, 0.0, NULL, NULL) );

            /* Here we calculate the MOL% Node mixing ratio. */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-7-NodeMixingRatioMOLCalculation_%s", Node->id);
            SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingConsMol, consname, expFin, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, mixingConsMol) );
            SCIP_CALL( SCIPreleaseCons(scip, &mixingConsMol) );

            /* free all expressions used to form root expression */
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprC2) );
            SCIP_CALL( SCIPreleaseExpr(scip, &C2inv) );
            SCIP_CALL( SCIPreleaseExpr(scip, &C1fracC2) );
            SCIP_CALL( SCIPreleaseExpr(scip, &expMuMolNode) );
            SCIP_CALL( SCIPreleaseExpr(scip, &expFin) );
         }
      }
   }
#endif /* SCIP_VERSION >= 800 */

   return SCIP_OKAY;
}


/** generate mixing ratio calculation "constraint" for each node */
static
SCIP_RETCODE generateFlowConsMixingNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
#if SCIP_VERSION >= 800
   char consname[SCIP_MAXSTRLEN];
   GAS_Network* Network = probdata->network;
   SCIP_CONS* mixingCons;
   int j;
   int c = 0;
   int sum = 0;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   for (j = 0; j < Network->numnodes; j++)
   {
      GAS_Node* Node = &(Network->nodes_ptr[j]);
      GAS_Arc* arc;
      int numInarcs = 0;
      int numOutarcs = 0;

      assert( Node != NULL );

      /* Determine the number of in and outarcs */
      arc = Node->inarcs;
      while ( arc != NULL )
      {
         numInarcs++;
         arc = arc->next_inarc;
      }

      arc = Node->outarcs;
      while ( arc != NULL )
      {
         numOutarcs++;
         arc = arc->next_outarc;
      }

      sum = numInarcs + numOutarcs;

      if ( sum == 1 && Node->type == ENTRY )
      {

         assert( Node->nodeMixingRatio != NULL );

         if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
         {
            SCIPdebugMsg(scip, "Heavy Gas, i.e., fix mu = 1.\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 1.0) );

            if ( ! probdata->linearMixing )
            {
               assert( Node->nodeMixingRatioMol != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatioMol, 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatioMol, 1.0) );
            }
         }
         else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
         {
            SCIPdebugMsg(scip, "Light Gas i.e. fix mu = 0\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 0.0) );

            if ( ! probdata->linearMixing )
            {
               assert( Node->nodeMixingRatioMol != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatioMol, 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatioMol, 0.0) );
            }
         }
         /* otherwise the gas is a mixture of the gaes with the highest and lowest molar mass and we calculate mu */
         else
         {
            SCIP_Real mu;
            SCIP_Real muMol;
            SCIP_Real M1 = 1.0 / probdata->molarMass1;
            SCIP_Real M2 = 1.0 / probdata->molarMass2;

            /* calculate the mass% mu */
            mu = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);

            SCIPdebugMsg(scip, "Mixed Gas, molar mass = %f, mu = %f", Node->molarMass, mu);
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, mu) );

            if ( ! probdata->linearMixing )
            {
               /* calculate mol% mu */
               muMol = (mu * M1) / (mu * (M1 - M2) + M2);
               SCIPdebugMsg(scip, "mu mol = %f\n", muMol);

               SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, muMol) );
               SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, muMol) );
            }
         }
      }

      if ( sum == 1 && Node->type == EXIT )
      {
         SCIP_CONS* fixNodeMuCons;
         SCIP_VAR* vars[2];
         SCIP_Real vals[2];

         vals[0] = 1.0;
         vals[1] = -1.0;

         if ( Node->inarcs != NULL )
         {
            arc = Node->inarcs;
            vars[0] = arc->sourcenode->nodeMixingRatio;
            assert( arc->sourcenode->nodeMixingRatio != NULL );
            SCIPdebugMsg(scip, "Arc is inarc\n");
         }
         else
         {
            arc = Node->outarcs;
            vars[0] = arc->targetnode->nodeMixingRatio;
            assert( arc->targetnode->nodeMixingRatio != NULL );
            SCIPdebugMsg(scip, "Arc is outarc\n");
         }

         assert( Node->nodeMixingRatio != NULL );
         vars[1] = Node->nodeMixingRatio;

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-3-EXIT-fixNodeMixRatio#%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixNodeMuCons, consname, 2, vars, vals, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, fixNodeMuCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &fixNodeMuCons) );
      }
      /* node has 2 or more neighbours and mu cannot be fixed */
      else if ( sum >= 2 )
      {
         SCIP_EXPR* exprQZMuprod[2];
         SCIP_EXPR** exprQZMuprodResult = NULL;
         SCIP_EXPR* exprSUMinMu = NULL;
         SCIP_EXPR* exprQaOut = NULL;
         SCIP_EXPR** exprQZMuprodResultO = NULL;
         SCIP_EXPR* exprSUMoutMu = NULL;
         SCIP_EXPR* exprZsum = NULL;
         SCIP_EXPR* exprZarr[2];
         SCIP_EXPR* exprMuNode = NULL;
         SCIP_Real coefPlusOneMinusOne[2];

         /* used in the loops below as offset i.e if the node is an entry the inflow is added to the flow sums */
         SCIP_Real inflow = 0.0;
         SCIP_Real inflowMU = 0.0;
         SCIP_Real copyInflow = 0.0;
         SCIP_Real copyInflowMU = 0.0;

         coefPlusOneMinusOne[0] = 1.0;
         coefPlusOneMinusOne[1] = -1.0;

         /* Adding the inflow at an as entry offset is complicated and happens as follows: if either ONLY inarcs or
          * outarcs exist, the offset is added in one of the next two if loops. */
         /* Otherwise the offset is added in the 3rd if-loop (i.e. if numInarcs && numoutarcs > 0) */
         /* The complication is caused by the fact that the flow value has to be added exactly once. */
         /* The flow is always added in the denominator but only in the numerator if mu=1. */
         if ( Node->type == ENTRY )
         {
            if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
            {
               /* constituent 1 => mu=1 */
               inflow = Node->flow;
               inflowMU = 1;
            }
            else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
            {
               /* constituent 0 => mu=0 */
               inflow = Node->flow;
               inflowMU = 0;
            }
            /* gas is a mixture and we first need to calculate mu */
            else
            {
               /* calculate the mass% mu */
               inflowMU = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);
               inflow = Node->flow;
               SCIPdebugMsg(scip, "inflowMU = %f, molar mass = %f, node = %s\n", inflowMU, Node->molarMass, Node->id);
            }
         }

         if ( numInarcs > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResult, numInarcs) );

            /* reset counter */
            c = 0;

            /* fill in incoming arcs */
            arc = Node->inarcs;
            while ( arc != NULL )
            {
               SCIP_EXPR* exprQaIn = NULL;
               SCIP_EXPR* exprZaIn = NULL;
               SCIP_EXPR* exprZaNegIn = NULL;
               SCIP_EXPR* exprMUaIn = NULL;
               SCIP_EXPR* exprMUaTargIn = NULL;
               SCIP_EXPR* exprPosMUprod = NULL;
               SCIP_EXPR* exprNegMUprod = NULL;
               SCIP_EXPR* exprSumMU = NULL;
               SCIP_EXPR* exprMUTermArr[2];

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaIn, arc->flowvar, NULL, NULL) );
               assert( arc->positiveFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaIn, arc->positiveFlowBinvar, NULL, NULL) );
               assert( arc->negativeFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaNegIn, arc->negativeFlowBinvar, NULL, NULL) );
               assert( arc->sourcenode->nodeMixingRatio != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaIn, arc->sourcenode->nodeMixingRatio, NULL, NULL) );
               assert( arc->targetnode->nodeMixingRatio != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaTargIn, arc->targetnode->nodeMixingRatio, NULL, NULL) );

               /* mu_sourcenode * positiveFlowBinvar */
               exprMUTermArr[0] = exprMUaIn;
               exprMUTermArr[1] = exprZaIn;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprPosMUprod, 2, exprMUTermArr, 1.0, NULL, NULL) );

               /* mu_targetnode * negativeFlowBinvar */
               exprMUTermArr[0] = exprMUaTargIn;
               exprMUTermArr[1] = exprZaNegIn;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprNegMUprod, 2, exprMUTermArr, 1.0, NULL, NULL) );

               /* sum = mu_sourcenode*pos + mu_targetnode*neg */
               {
                  SCIP_EXPR* exprSumArr[2];
                  exprSumArr[0] = exprPosMUprod;
                  exprSumArr[1] = exprNegMUprod;
                  SCIP_CALL( SCIPcreateExprSum(scip, &exprSumMU, 2, exprSumArr, NULL, 0.0, NULL, NULL) );
               }

               /* final product q_a * (mu_s*pos + mu_t*neg) */
               exprQZMuprod[0] = exprQaIn;
               exprQZMuprod[1] = exprSumMU;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResult[c], 2, exprQZMuprod, 1.0, NULL, NULL) );

               /* release intermediates */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaNegIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaTargIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprPosMUprod) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprNegMUprod) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprSumMU) );

               arc = arc->next_inarc;
               c++;
            }

            assert( c == numInarcs );
            c = 0;

            /* need to add the inflow value as offset - see details above */
            if ( numOutarcs == 0 )
            {
               copyInflow = inflow;
               copyInflowMU = inflowMU;
            }

            /* now the sum over the above created producs is build */

            /* Sum over inarcs q_a*z_a*mu_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMinMu, numInarcs, exprQZMuprodResult, NULL, copyInflowMU * copyInflow, NULL, NULL) );

            /*free arrays*/
            for (int i = 0; i < numInarcs; i++)
            {
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZMuprodResult[i]) );
            }

            SCIPfreeBufferArray(scip, &exprQZMuprodResult);
         }

         if ( numOutarcs > 0 )
         {
            /* Analogous calculation now for the out arcs. */
            /* The loop below builds the pairwise producs q_a*z_a and q_a*z_a*mu_a for the outarcs. */
            /* NOTE: MASS% mixing Ratio is used i.e. mixingRatioMass. */
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResultO, numOutarcs) );

            c = 0;
            arc = Node->outarcs;

            /* fill in outgoing arcs */
            while ( arc != NULL )
            {
               /* Build expression q_a*( mu_sourcenode*positiveFlowBinvar + mu_targetnode*negativeFlowBinvar ) for outarcs */
               SCIP_EXPR* exprZaPosOut = NULL;
               SCIP_EXPR* exprZaNegOut = NULL;
               SCIP_EXPR* exprMUaSourOut = NULL;
               SCIP_EXPR* exprMUaTargOut = NULL;
               SCIP_EXPR* exprPosMUprodOut = NULL;
               SCIP_EXPR* exprNegMUprodOut = NULL;
               SCIP_EXPR* exprSumMUOut = NULL;
               SCIP_EXPR* exprMUTermArrOut[2];

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaOut, arc->flowvar, NULL, NULL) );
               assert( arc->positiveFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaPosOut, arc->positiveFlowBinvar, NULL, NULL) );
               assert( arc->negativeFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaNegOut, arc->negativeFlowBinvar, NULL, NULL) );
               assert( arc->sourcenode->nodeMixingRatio != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaSourOut, arc->sourcenode->nodeMixingRatio, NULL, NULL) );
               assert( arc->targetnode->nodeMixingRatio != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaTargOut, arc->targetnode->nodeMixingRatio, NULL, NULL) );

               /* mu_source * positiveFlowBinvar */
               exprMUTermArrOut[0] = exprMUaSourOut;
               exprMUTermArrOut[1] = exprZaPosOut;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprPosMUprodOut, 2, exprMUTermArrOut, 1.0, NULL, NULL) );

               /* mu_target * negativeFlowBinvar */
               exprMUTermArrOut[0] = exprMUaTargOut;
               exprMUTermArrOut[1] = exprZaNegOut;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprNegMUprodOut, 2, exprMUTermArrOut, 1.0, NULL, NULL) );

               /* sum = mu_source*pos + mu_target*neg */
               {
                  SCIP_EXPR* exprSumArrOut[2];
                  exprSumArrOut[0] = exprPosMUprodOut;
                  exprSumArrOut[1] = exprNegMUprodOut;
                  SCIP_CALL( SCIPcreateExprSum(scip, &exprSumMUOut, 2, exprSumArrOut, NULL, 0.0, NULL, NULL) );
               }

               /* final product q_a * (sum) */
               exprQZMuprod[0] = exprQaOut;
               exprQZMuprod[1] = exprSumMUOut;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResultO[c], 2, exprQZMuprod, 1.0, NULL, NULL) );

               /* release intermediates */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaPosOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaNegOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaSourOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaTargOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprPosMUprodOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprNegMUprodOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprSumMUOut) );

               arc = arc->next_outarc;
               c++;
            }

            assert( c == numOutarcs );

            /* In the case that the node has only inflow on outarcs i.e. the flow is negative we multiply it by -1. */
            SCIP_Real minusOne[numOutarcs];

            if ( numInarcs == 0 )
            {
               copyInflow = inflow;
               copyInflowMU = inflowMU;

               for (int i = 0; i < numOutarcs; i++)
                  minusOne[i] = -1.0;
            }
            else
            {
               for (int i = 0; i < numOutarcs; i++)
                  minusOne[i] = 1.0;
            }

            /* now the sums over the above created producs are built */

            /* sum over inarcs q_a * z_a * mu_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMoutMu, numOutarcs, exprQZMuprodResultO, minusOne, copyInflowMU * copyInflow, NULL, NULL) );

            /* free arrays */
            for (int i = 0; i < numOutarcs; i++)
            {
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZMuprodResultO[i]) );
            }

            SCIPfreeBufferArray(scip, &exprQZMuprodResultO);
         }

         if ( numInarcs > 0 && numOutarcs > 0 )
         {
            if( Node->type != EXIT )
            {
               exprZarr[0] = exprSUMinMu;
               exprZarr[1] = exprSUMoutMu;
               SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 2, exprZarr, coefPlusOneMinusOne, inflowMU * inflow, NULL, NULL) );

               /* constraint */
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprZsum, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );
            }
            else
            {
               SCIP_EXPR* exprSum3[3];
               SCIP_Real coef3[3];

               coef3[0] = 1.0;
               coef3[1] = -1.0;
               coef3[2] = Node->flow;
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );

               exprSum3[0] = exprSUMinMu;
               exprSum3[1] = exprSUMoutMu;
               exprSum3[2] = exprMuNode;
               SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 3, exprSum3, coef3, 0.0, NULL, NULL) );

               SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );

               (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id );
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprZsum, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );
            }
         }
         else if ( numInarcs > 0 && numOutarcs == 0 )
         {

            if ( Node->type != EXIT )
            {
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprSUMinMu, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );
            }
            else
            {
               SCIP_EXPR* exprSum3[3];
               SCIP_Real coef3[3];
               coef3[0] = 1.0;
               coef3[1] = Node->flow;
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );

               exprSum3[0] = exprSUMinMu;
               exprSum3[1] = exprMuNode;
               SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 2, exprSum3, coef3, 0.0, NULL, NULL) );

               SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprZsum, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );
            }
         }
         else if ( numInarcs == 0 && numOutarcs > 0 )
         {
            if ( Node->type != EXIT )
            {
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprSUMoutMu, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );
            }
            else
            {
               SCIP_EXPR* exprSum3[3];
               SCIP_Real coef3[3];

               coef3[0] = 1.0;
               coef3[1] = Node->flow;
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );

               exprSum3[0] = exprSUMoutMu;
               exprSum3[1] = exprMuNode;
               SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 2, exprSum3, coef3, 0.0, NULL, NULL) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprZsum, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );
            }
         }

         /* free all expr */
         if ( numInarcs > 0 && numOutarcs > 0 )
         {
            /*release expr*/
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMinMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMoutMu) );;
            SCIP_CALL( SCIPreleaseExpr(scip, &exprZsum) );
         }
         else if ( numInarcs > 0 && numOutarcs == 0 )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMinMu) );
            if( Node->type == EXIT )
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZsum) );
         }
         else if ( numInarcs == 0 && numOutarcs > 0 )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMoutMu) );
            if( Node->type == EXIT )
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZsum) );
         }
      }
   }
#endif /* SCIP_VERSION >= 800 */

   return SCIP_OKAY;
}


/** generate mixing ratio caclulation "constraint" for each arc and node */
static
SCIP_RETCODE generateMixingRatio(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
#if SCIP_VERSION >= 800
   GAS_Network* Network = probdata->network;
   SCIP_CONS* mixingCons;
   char consname[SCIP_MAXSTRLEN];
   int j;
   int c = 0;
   int sum = 0;
   SCIP_Bool foundCS = FALSE;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   for (j = 0; j < Network->numnodes; j++)
   {
      GAS_Node* Node = &(Network->nodes_ptr[j]);
      GAS_Arc* arc;
      int numInarcs = 0;
      int numOutarcs = 0;

      assert( Node != NULL );

      /* Determine the number of in and outarcs */
      arc = Node->inarcs;
      while ( arc != NULL )
      {
         numInarcs++;
         arc = arc->next_inarc;
      }

      arc = Node->outarcs;
      while ( arc != NULL )
      {
         numOutarcs++;
         arc = arc->next_outarc;
      }

      sum = numInarcs + numOutarcs;

      if ( sum == 1 && Node->type == ENTRY )
      {

         if ( Node->inarcs != NULL )
            arc = Node->inarcs;
         else
            arc = Node->outarcs;

         assert( Node->nodeMixingRatio != NULL );
         assert( arc->mixingRatioMass != NULL );

         if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
         {
            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 1, i.e., fix mu = 1\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 1.0) );

            /* Second the arc_mixingRatioMass can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatioMass, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatioMass, 1.0) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatio, 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatio, 1.0) );
            }
         }
         else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
         {
            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 0, i.e., fix mu = 0\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 0.0) );

            /* Second the arc_mixingRatioMass can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatioMass, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatioMass, 0.0) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatio, 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatio, 0.0) );
            }
         }
         /* otherwise the gas is a mixture and we need to calculate mu */
         else
         {
            SCIP_Real mu;
            SCIP_Real muMol;
            SCIP_Real M1 = 1.0 / probdata->molarMass1;
            SCIP_Real M2 = 1.0 / probdata->molarMass2;

            /* calculate the mass% mu */
            mu = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);
            /* calculate mol% mu */
            muMol = (mu * M1) / (mu * (M1 - M2) + M2);
            SCIPdebugMsg(scip, "Mixed Gas, molarMass = %f, mu = %f, muMol = %f\n", Node->molarMass, mu, muMol);

            /* First the node mixing ratio can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, mu) );

            /* Second the arc_mixingRatioMass can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatioMass, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatioMass, mu) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatio, muMol) );
               SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatio, muMol) );
            }
         }
      }

      /* If an entry only has a compressor as outarc then the above won't detect it because of the bypass arc. We have to check this here. */
      /* only check outarcs because there is no in-arc compressor */
      arc = Node->outarcs;

      while ( arc != NULL )
      {
         if ( arc->type == CS )
            foundCS = TRUE;
         arc = arc->next_outarc;
      }

      /* only continue if a CS is found */
      if ( sum == 2 && Node->type == ENTRY && foundCS )
      {
         GAS_Arc* arc1 = NULL;
         GAS_Arc* arc2 = NULL;
         SCIP_Bool secondRun = FALSE;

         arc = Node->outarcs;
         /* set arc1 and arc2 */

         while ( arc != NULL )
         {
            if ( ! secondRun )
            {
               arc1 = arc;
               secondRun = TRUE;
               arc = arc->next_outarc;
            }
            else if ( secondRun )
            {
               arc2 = arc;
               arc = arc->next_outarc;
            }
         }

         assert( Node->nodeMixingRatio != NULL );
         assert( arc->mixingRatioMass != NULL );

         /* Generate constraints to fix mu on Node, CS and Bypass */
         if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
         {
            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 1, i.e., fix mu = 1\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 1.0) );

            /* Second the arc_mixingRatioMass can be fixed */
            assert( arc1->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatioMass, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatioMass, 1.0) );
            assert( arc2->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatioMass, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatioMass, 1.0) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               assert( arc1->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatio, 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatio, 1.0) );
               assert( arc2->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatio, 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatio, 1.0) );
            }
         }
         else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2 ) )
         {

            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 0, i.e., fix mu = 0\n");
            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 1, i.e., fix mu = 1\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 0.0) );

            /* Second the arc_mixingRatioMass can be fixed */
            assert( arc1->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatioMass, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatioMass, 0.0) );
            assert( arc2->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatioMass, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatioMass, 0.0) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               assert( arc1->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatio, 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatio, 0.0) );
               assert( arc2->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatio, 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatio, 0.0) );
            }
         }
         /* otherwise the gas is a mixture and we need to calculate mu */
         else
         {
            SCIP_Real mu;
            SCIP_Real muMol;
            SCIP_Real M1 = 1.0 / probdata->molarMass1;
            SCIP_Real M2 = 1.0 / probdata->molarMass2;

            /* calculate the mass% mu */
            mu = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);
            /* calculate mol% mu */
            muMol = (mu * M1) / (mu * (M1 - M2) + M2);
            SCIPdebugMsg(scip, "Mixed Gas, molarMass = %f, mu = %f, muMol = %f\n", Node->molarMass, mu, muMol );

            /* First the node mixing ratio can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, mu) );

            /* Second the arc_mixingRatioMass can be fixed */
            assert( arc1->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatioMass, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatioMass, mu) );
            assert( arc2->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatioMass, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatioMass, mu) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               assert( arc1->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatio, muMol) );
               SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatio, muMol) );
               assert( arc2->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatio, muMol) );
               SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatio, muMol) );
            }
         }
      }
      /* if sum of in and outarcs = 2 and node is no entry OR exit, we can fix mu and the nodeMixingRatio */
      else if ( sum == 2 && Node->type == INNODE )
      {
         SCIP_CONS* fixMuCons;
         SCIP_CONS* fixMuMassCons;
         SCIP_CONS* fixMuNodeCons;
         SCIP_VAR* vars[2];
         SCIP_Bool secondRun = FALSE;
         SCIP_Real vals[2];
         GAS_Arc* arc1 = NULL;
         GAS_Arc* arc2 = NULL;
         GAS_Arc* arcs[2];

         arcs[0] = Node->outarcs;
         arcs[1] = Node->inarcs;
         vals[0] = 1;
         vals[1] = -1;

         /* set arc1 and arc2 */
         for (int i = 0; i < 2; i++)
         {
            arc = arcs[i];
            while ( arc != NULL )
            {
               /* Check if outarc */
               if ( i == 0 && !secondRun )
               {
                  arc1 = arc;
                  secondRun = TRUE;
                  arc = arc->next_outarc;
               }
               else if ( i == 0 && secondRun )
               {
                  arc2 = arc;
                  arc = arc->next_outarc;
               }

               /* Check if inarc */
               if ( i == 1 && !secondRun )
               {
                  arc1 = arc;
                  secondRun = TRUE;
                  arc = arc->next_inarc;
               }
               else if ( i == 1 && secondRun )
               {
                  arc2 = arc;
                  arc = arc->next_inarc;
               }
            }
         }

         /* Now that arc1 and arc2 are set, we generate the contraint */
         /* mu_arc1 = mu_arc2 */
         if ( ! probdata->linearMixing )
         {
            assert( arc1 != NULL );
            assert( arc2 != NULL );
            assert( arc1->mixingRatio != NULL );
            assert( arc2->mixingRatio != NULL );

            vars[0] = arc1->mixingRatio;
            vars[1] = arc2->mixingRatio;

            /* Create the first constraint for Mu and add it to SCIP */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1A-2Arcs-fixMuCons#%s", Node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuCons, consname, 2, vars, vals, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, fixMuCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &fixMuCons) );
         }

         assert( arc1->mixingRatioMass != NULL );
         assert( arc2->mixingRatioMass != NULL );

         vars[0] = arc1->mixingRatioMass;
         vars[1] = arc2->mixingRatioMass;

         /* Create the second constraint for MuMass and add it to SCIP */
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1A-2Arcs-fixMuMassCons#%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuMassCons, consname, 2, vars, vals, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, fixMuMassCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &fixMuMassCons) );

         /* fix mu node aka nodeMixingRatio */
         assert( Node->nodeMixingRatio != NULL );
         vars[0] = Node->nodeMixingRatio;
         vars[1] = arc1->mixingRatioMass; /* does not matter if arc1 or arc2 is used */
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1N-2Arcs-fixNodeMu#%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuNodeCons, consname, 2, vars, vals, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, fixMuNodeCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &fixMuNodeCons) );
      }
      /* node has 2 or more neighbours and mu cannot be fixed */
      else if ( sum >= 2 )
      {
         SCIP_EXPR* exprQZprod[2];
         SCIP_EXPR** exprQZprodResult = NULL;
         SCIP_EXPR* exprQZMuprod[2];
         SCIP_EXPR** exprQZMuprodResult = NULL;
         SCIP_EXPR* exprSUMin = NULL;
         SCIP_EXPR* exprSUMinMu = NULL;
         SCIP_EXPR** exprQZprodResultO = NULL;
         SCIP_EXPR** exprQZMuprodResultO = NULL;
         SCIP_EXPR* exprSUMout = NULL;
         SCIP_EXPR* exprSUMoutMu = NULL;
         SCIP_EXPR* exprZsum = NULL;
         SCIP_EXPR* exprZarr[2];
         SCIP_Real coefPlusOneMinusOne[2];
         SCIP_EXPR* exprNsum = NULL;
         SCIP_EXPR* exprN[2];
         SCIP_EXPR* exprMuNode = NULL;
         SCIP_EXPR* exprProdN[2];
         SCIP_EXPR* exprMuSumProd = NULL;
         SCIP_EXPR* exprFsumm[2];
         SCIP_EXPR* exprFinMu;
         /* used in the loops below as offset i.e if the node is an entry the inflow is added to the flow sums */
         SCIP_Real inflow = 0.0;
         SCIP_Real inflowMU = 0.0;
         SCIP_Real copyInflow = 0.0;
         SCIP_Real copyInflowMU = 0.0;

         coefPlusOneMinusOne[0] = 1;
         coefPlusOneMinusOne[1] = -1;

         /* Adding the inflow at an as entry offset is complicated and happens as follows: if either ONLY inarcs or
          * outarcs exist, the offset is added in one of the next two if loops.*/
         /* Otherwise the offset is added in the 3rd if loop (i.e. if numInarcs && numoutarcs >0). */
         /* The complication is cauesd by the fact that the flow value has to be added exactly once. */
         /* The flow is always added in the denominator but only in the numerator if mu=1. */
         if ( Node->type == ENTRY )
         {
            if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
            {
               /* heavy gas => mu=1 */
               inflow = Node->flow;
               inflowMU = 1;
            }
            else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
            {
               /* light gas => mu=0 */
               inflow = Node->flow;
               inflowMU = 0;
            }
            /* gas is a mixture and we first need to calculate mu */
            else
            {
               /* calculate the mass% mu */
               inflowMU = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);
               inflow = Node->flow;
               SCIPdebugMsg(scip, "inflow on node with >=2 arcs and mu= %f, molarMass = %f, node =%s\n", inflowMU, Node->molarMass, Node->id );
            }
         }

         if ( numInarcs > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZprodResult, numInarcs) );
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResult, numInarcs) );

            /* The loop below builds the pairwise producs q_a*z_a and q_a*z_a*mu_a for the inarcs. */
            /* NOTE: MASS% mixing Ratio is used i.e. mixingRatioMass */

            /* reset counter */
            c = 0;

            /* fill in incoming arcs */
            arc = Node->inarcs;
            while ( arc != NULL )
            {
               SCIP_EXPR* exprQaIn = NULL;
               SCIP_EXPR* exprZaIn = NULL;
               SCIP_EXPR* exprMUaIn = NULL;

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaIn, arc->flowvar, NULL, NULL) );
               assert( arc->positiveFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaIn, arc->positiveFlowBinvar, NULL, NULL) );

               exprQZprod[0] = exprQaIn;
               exprQZprod[1] = exprZaIn;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZprodResult[c], 2, exprQZprod, 1.0, NULL, NULL) ); /* Pairwise q_a*z_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaIn) );

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaIn, arc->flowvar, NULL, NULL) );
               assert( arc->positiveFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaIn, arc->positiveFlowBinvar, NULL, NULL) );
               assert( arc->mixingRatioMass != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaIn, arc->mixingRatioMass, NULL, NULL) );

               exprQZMuprod[0] = exprQZprodResult[c];
               exprQZMuprod[1] = exprMUaIn;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResult[c], 2, exprQZMuprod, 1.0, NULL, NULL) ); /* Pairwise q_a*z_a*mu_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaIn) );

               arc = arc->next_inarc;
               c++;
            }

            assert( c == numInarcs );
            c = 0;

            /* need to add the inflow value as offset - see details above */
            if ( numOutarcs == 0 )
            {
               copyInflow = inflow;
               copyInflowMU = inflowMU;
            }

            /* now the sums over the above created producs are build */

            /* Sum over inarcs q_a*z_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMin, numInarcs, exprQZprodResult, NULL, copyInflow, NULL, NULL) );

            /* Sum over inarcs q_a*z_a*mu_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMinMu, numInarcs, exprQZMuprodResult, NULL, copyInflowMU * copyInflow, NULL, NULL) );

            /*free arrays*/
            for (int i = 0; i < numInarcs; i++)
            {
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZprodResult[i]) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZMuprodResult[i]) );
            }

            SCIPfreeBufferArray(scip, &exprQZMuprodResult);
            SCIPfreeBufferArray(scip, &exprQZprodResult);
         }

         if ( numOutarcs > 0 )
         {
            /* Analog calculation now for the out arcs. */
            /* The loop below builds the pairwise producs q_a*z_a and q_a*z_a*mu_a for the outarcs. */
            /* NOTE: MASS% mixing Ratio is used i.e. mixingRatioMass. */
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZprodResultO, numOutarcs) );
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResultO, numOutarcs) );

            c = 0;
            arc = Node->outarcs;

            /* fill in outgoing arcs */
            while ( arc != NULL )
            {
               SCIP_EXPR* exprQaOut = NULL;
               SCIP_EXPR* exprZaOut = NULL;
               SCIP_EXPR* exprMUaOut = NULL;

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaOut, arc->flowvar, NULL, NULL) );
               assert( arc->negativeFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaOut, arc->negativeFlowBinvar, NULL, NULL) );
               exprQZprod[0] = exprQaOut;
               exprQZprod[1] = exprZaOut;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZprodResultO[c], 2, exprQZprod, 1.0, NULL, NULL) ); /* pairwise q_a * z_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaOut) );

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaOut, arc->flowvar, NULL, NULL) );
               assert( arc->negativeFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaOut, arc->negativeFlowBinvar, NULL, NULL) );
               assert( arc->mixingRatioMass != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaOut, arc->mixingRatioMass, NULL, NULL) );
               exprQZMuprod[0] = exprQZprodResultO[c];
               exprQZMuprod[1] = exprMUaOut;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResultO[c], 2, exprQZMuprod, 1.0, NULL, NULL) ); /* pairwise q_a * z_a * mu_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZaOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaOut) );

               arc = arc->next_outarc;
               c++;
            }

            assert( c == numOutarcs );

            /* In the case that the node has only inflow on outarcs i.e. the flow is negative we multiply it by -1. */
            SCIP_Real minusOne[numOutarcs];

            if ( numInarcs == 0 )
            {
               copyInflow = inflow;
               copyInflowMU = inflowMU;

               for (int i = 0; i < numOutarcs; i++)
                  minusOne[i] = -1.0;
            }
            else
            {
               for (int i = 0; i < numOutarcs; i++)
                  minusOne[i] = 1.0;
            }

            /* now the sums over the above created producs are built */

            /* Sum over inarcs q_a*z_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMout, numOutarcs, exprQZprodResultO, minusOne, copyInflow, NULL, NULL) );

            /* Sum over inarcs q_a*z_a*mu_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMoutMu, numOutarcs, exprQZMuprodResultO, minusOne, copyInflowMU * copyInflow, NULL, NULL) );

            /*free arrays*/
            for (int i = 0; i < numOutarcs; i++)
            {
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZprodResultO[i]) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZMuprodResultO[i]) );
            }
            SCIPfreeBufferArray(scip, &exprQZMuprodResultO);
            SCIPfreeBufferArray(scip, &exprQZprodResultO);
         }

         /* In the following the mixing equation (Equation (7) of the mixing paper) is built. */
         if ( numInarcs > 0 && numOutarcs > 0 )
         {
            /* Now we first build the numerator for the mixing eqaution and then the denominator down below. */
            /* Note that the equation will be written in non fractional form (no div by 0). */

            /* Numerator sum of the mixing equation - EQ (7) in the paper. */
            exprZarr[0] = exprSUMinMu;
            exprZarr[1] = exprSUMoutMu;
            SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 2, exprZarr, coefPlusOneMinusOne, inflowMU *inflow, NULL, NULL) );

            /* Mu_Node*NominatorSum - Denominator sum muliplied by mu */
            exprN[0] = exprSUMin;
            exprN[1] = exprSUMout;
            SCIP_CALL( SCIPcreateExprSum(scip, &exprNsum, 2, exprN, coefPlusOneMinusOne, inflow, NULL, NULL) );

            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );
            exprProdN[0] = exprMuNode;
            exprProdN[1] = exprNsum;
            SCIP_CALL( SCIPcreateExprProduct(scip, &exprMuSumProd, 2, exprProdN, 1.0, NULL, NULL) );

            exprFsumm[0] = exprMuSumProd;
            exprFsumm[1] = exprZsum;
         }
         else if ( numInarcs > 0 && numOutarcs == 0)
         {
            /* Denominator is only sum over inarcs = exprSUMinMu */

            /* Mu_Node * NominatorSum */
            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );
            exprProdN[0] = exprMuNode;
            exprProdN[1] = exprSUMin;
            SCIP_CALL( SCIPcreateExprProduct(scip, &exprMuSumProd, 2, exprProdN, 1.0, NULL, NULL) );

            exprFsumm[0] = exprMuSumProd;
            exprFsumm[1] = exprSUMinMu;
         }
         else if ( numInarcs == 0 && numOutarcs > 0 )
         {
            /* Denominator is only sum over outarcs = exprSUMoutMu */

            /* Mu_Node * NominatorSum. Note that the inflow here is negative which is why it is multiplied with -1 above !*/
            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );
            exprProdN[0] = exprMuNode;
            exprProdN[1] = exprSUMout;
            SCIP_CALL( SCIPcreateExprProduct(scip, &exprMuSumProd, 2, exprProdN, 1.0, NULL, NULL) );

            exprFsumm[0] = exprMuSumProd;
            exprFsumm[1] = exprSUMoutMu;
         }

         /* Final equation */
         /* Mu_Node*NominatorSum - DenominatorSum = 0 */
         SCIP_CALL( SCIPcreateExprSum(scip, &exprFinMu, 2, exprFsumm, coefPlusOneMinusOne, 0.0, NULL, NULL) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-NodeMixingRatioCalculation#%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprFinMu, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, mixingCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );

         /* free all expr */
         if ( numInarcs > 0 && numOutarcs > 0 )
         {
            /*release expr*/
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMinMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMoutMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMin) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMout) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprZsum) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuSumProd) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprFinMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprNsum) );
         }
         else if ( numInarcs > 0 && numOutarcs == 0 )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMinMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMin) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuSumProd) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprFinMu) );
         }
         else if ( numInarcs == 0 && numOutarcs > 0 )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMout) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMoutMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuNode) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprMuSumProd) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprFinMu) );
         }

         {
            /* Now we iterate over all arcs and first calculate and set the molar percentage mixing ratio for each arc. */
            /* Afterwards we set nodeMixingRatio= mixingRatioMass (Arc). */
            /* Note that each constraint is generated twice since we do not know the flow direction in advance. */
            /* Coupling with the flowBinVars ensures that only the correct is active. */
            GAS_Arc* arcs[2];
            SCIP_CONS* linkConsIn;
            SCIP_CONS* linkConsMuMass;

            coefPlusOneMinusOne[0] = 1;
            coefPlusOneMinusOne[1] = -1;

            arcs[0] = Node->outarcs;
            arcs[1] = Node->inarcs;

            c = 0; /* counter */

            for (int i = 0; i < 2; i++)
            {
               arc = arcs[i];

               while ( arc != NULL )
               {
                  SCIP_EXPR* exprSum[2];
                  SCIP_EXPR* exprFIN;
                  SCIP_EXPR* exprZ;
                  SCIP_EXPR* exprMuArc;
                  SCIP_EXPR* exprMuMassArc;
                  SCIP_EXPR* exprDiff;
                  SCIP_EXPR* exprDiffZ;
                  SCIP_EXPR* exprMixRatioNode;

                  if ( ! probdata->linearMixing )
                  {
                     /* Use binary variable in order to determine the arcs with outflow. */
                     if ( i == 0 )
                     {
                        /* flow on an outarc really flows out */
                        assert( arc->positiveFlowBinvar != NULL );
                        SCIP_CALL( SCIPcreateExprVar(scip, &exprZ, arc->positiveFlowBinvar, NULL, NULL) );
                     }
                     else
                     {
                        /* flow on an inarc flows out */
                        assert( arc->negativeFlowBinvar != NULL );
                        SCIP_CALL( SCIPcreateExprVar(scip, &exprZ, arc->negativeFlowBinvar, NULL, NULL) );
                     }

                     /* calculating the mole% mr for the arcs. It still gets double linked and the binary variables are used. */
                     {
                        SCIP_EXPR* expr[2];
                        SCIP_Real coef[1];
                        SCIP_EXPR* C2inv;
                        SCIP_EXPR* exprC2;
                        SCIP_EXPR* C1fracC2;
                        SCIP_EXPR* link;
                        SCIP_Real M2 = 1.0 / probdata->molarMass2;

                        coef[0] = (1.0 / probdata->molarMass1) - (1.0 / probdata->molarMass2); /* this coef is negative */

                        assert( Node->nodeMixingRatio != NULL );
                        SCIP_CALL( SCIPcreateExprVar(scip, &exprMixRatioNode, Node->nodeMixingRatio, NULL, NULL) );

                        expr[0] = exprMixRatioNode;
                        SCIP_CALL( SCIPcreateExprSum(scip, &exprC2, 1, expr, coef, M2, NULL, NULL) );
                        SCIP_CALL( SCIPcreateExprPow(scip, &C2inv, exprC2, -1.0, NULL, NULL) );

                        SCIP_CALL( SCIPreleaseExpr(scip, &exprMixRatioNode) );
                        SCIP_CALL( SCIPcreateExprVar(scip, &exprMixRatioNode, Node->nodeMixingRatio, NULL, NULL) );

                        expr[0] = exprMixRatioNode;
                        expr[1] = C2inv;
                        SCIP_CALL( SCIPcreateExprProduct(scip, &C1fracC2, 2, expr, 1.0 / probdata->molarMass1, NULL, NULL) );

                        assert( arc->mixingRatio != NULL );
                        SCIP_CALL( SCIPcreateExprVar(scip, &exprMuArc, arc->mixingRatio, NULL, NULL) );

                        expr[0] = exprMuArc;
                        expr[1] = C1fracC2;
                        SCIP_CALL( SCIPcreateExprSum(scip, &link, 2, expr, coefPlusOneMinusOne, 0.0, NULL, NULL) );

                        expr[0] = link;
                        expr[1] = exprZ;
                        SCIP_CALL( SCIPcreateExprProduct(scip, &exprFIN, 2, expr, 1.0, NULL, NULL) );

                        /* Here we calculate the MOL% arc mixing ratio. */
                        (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-7Mol-linkNodeMixRatioToArcMixRatioAndCalcMolePercent#%s", arc->id);
                        SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &linkConsIn, consname, exprFIN, 0.0, 0.0) );
                        SCIP_CALL( SCIPaddCons(scip, linkConsIn) );
                        SCIP_CALL( SCIPreleaseCons(scip, &linkConsIn) );

                        /* free all expressions used to form root expression */
                        SCIP_CALL( SCIPreleaseExpr(scip, &exprMixRatioNode) );
                        SCIP_CALL( SCIPreleaseExpr(scip, &exprC2) );
                        SCIP_CALL( SCIPreleaseExpr(scip, &C2inv) );
                        SCIP_CALL( SCIPreleaseExpr(scip, &C1fracC2) );
                        SCIP_CALL( SCIPreleaseExpr(scip, &exprMuArc) );
                        SCIP_CALL( SCIPreleaseExpr(scip, &link) );
                        SCIP_CALL( SCIPreleaseExpr(scip, &exprFIN) );
                        SCIP_CALL( SCIPreleaseExpr(scip, &exprZ) );
                     }
                  }

                  /* Use binary variable in order to determine the arcs with outflow. */
                  if ( i == 0 )
                  {
                     /* flow on an outarc really flows out */
                     assert( arc->positiveFlowBinvar != NULL );
                     SCIP_CALL( SCIPcreateExprVar(scip, &exprZ, arc->positiveFlowBinvar, NULL, NULL) );
                  }
                  else
                  {
                     /* flow on an inarc flows out */
                     assert( arc->negativeFlowBinvar != NULL );
                     SCIP_CALL( SCIPcreateExprVar(scip, &exprZ, arc->negativeFlowBinvar, NULL, NULL) );
                  }

                  assert( Node->nodeMixingRatio != NULL );
                  SCIP_CALL( SCIPcreateExprVar(scip, &exprMixRatioNode, Node->nodeMixingRatio, NULL, NULL) );

                  /* Here we link nodeMixRatio to arcMixingRatioMass in the following way z_a * (nodeMu - arcMu) = 0. */
                  assert( arc->mixingRatioMass != NULL );
                  SCIP_CALL( SCIPcreateExprVar(scip, &exprMuMassArc, arc->mixingRatioMass, NULL, NULL) );
                  exprSum[0] = exprMixRatioNode;
                  exprSum[1] = exprMuMassArc;
                  SCIP_CALL( SCIPcreateExprSum(scip, &exprDiff, 2, exprSum, coefPlusOneMinusOne, 0.0, NULL, NULL) );
                  exprSum[0] = exprDiff;
                  exprSum[1] = exprZ;
                  SCIP_CALL( SCIPcreateExprProduct(scip, &exprDiffZ, 2, exprSum, 1.0, NULL, NULL) );

                  (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "TYPE-7Mass-linkNodeMixRatioToArcMixRatio#%s_%s", arc->id, Node->id );
                  SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &linkConsMuMass, consname, exprDiffZ, 0.0, 0.0) );
                  SCIP_CALL( SCIPaddCons(scip, linkConsMuMass) );
                  SCIP_CALL( SCIPreleaseCons(scip, &linkConsMuMass) );

                  /* free all expressions used to form root expression */
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprZ) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprMixRatioNode) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprMuMassArc) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprDiff) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprDiffZ) );

                  /* increment counter */
                  c++;

                  if ( i == 0 )
                     arc = arc->next_outarc;
                  else if ( i == 1 )
                     arc = arc->next_inarc;
               }
            }

            c = 0; /* reset counter */
         }
      }
   }
#endif /* SCIP_VERSION >= 800 */

   return SCIP_OKAY;
}

/** generate mixing ratio caclulation "constraint" for each arc and node */
static
SCIP_RETCODE generateFlowConsMixing(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
#if SCIP_VERSION >= 800
   GAS_Network* Network = probdata->network;
   SCIP_CONS* mixingCons;
   char consname[SCIP_MAXSTRLEN];
   int j;
   int c = 0;
   int sum = 0;
   SCIP_Bool foundCS = FALSE;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   for (j = 0; j < Network->numnodes; j++)
   {
      GAS_Node* Node = &(Network->nodes_ptr[j]);
      GAS_Arc* arc;
      int numInarcs = 0;
      int numOutarcs = 0;

      assert( Node != NULL );

      /* Determine the number of in and outarcs */
      arc = Node->inarcs;
      while ( arc != NULL )
      {
         numInarcs++;
         arc = arc->next_inarc;
      }

      arc = Node->outarcs;
      while ( arc != NULL )
      {
         numOutarcs++;
         arc = arc->next_outarc;
      }

      sum = numInarcs + numOutarcs;

      if ( sum == 1 && Node->type == ENTRY )
      {

         if ( Node->inarcs != NULL )
            arc = Node->inarcs;
         else
            arc = Node->outarcs;

         assert( Node->nodeMixingRatio != NULL );
         assert( arc->mixingRatioMass != NULL );

         if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
         {
            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 1, i.e., fix mu = 1\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 1.0) );

            /* Second the arc_mixingRatioMass can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatioMass, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatioMass, 1.0) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatio, 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatio, 1.0) );
            }
         }
         else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
         {
            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 0, i.e., fix mu = 0\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 0.0) );

            /* Second the arc_mixingRatioMass can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatioMass, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatioMass, 0.0) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatio, 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatio, 0.0) );
            }
         }
         /* otherwise the gas is a mixture and we need to calculate mu */
         else
         {
            SCIP_Real mu;
            SCIP_Real muMol;
            SCIP_Real M1 = 1.0 / probdata->molarMass1;
            SCIP_Real M2 = 1.0 / probdata->molarMass2;

            /* calculate the mass% mu */
            mu = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);
            /* calculate mol% mu */
            muMol = (mu * M1) / (mu * (M1 - M2) + M2);
            SCIPdebugMsg(scip, "Mixed Gas, molarMass = %f, mu = %f, muMol = %f\n", Node->molarMass, mu, muMol);

            /* First the node mixing ratio can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, mu) );

            /* Second the arc_mixingRatioMass can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatioMass, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatioMass, mu) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               SCIP_CALL( SCIPchgVarLb(scip, arc->mixingRatio, muMol) );
               SCIP_CALL( SCIPchgVarUb(scip, arc->mixingRatio, muMol) );
            }
         }
      }

      /* If an entry only has a compressor as outarc then the above won't detect it because of the bypass arc. We have to check this here. */
      /* only check outarcs because there is no in-arc compressor */
      arc = Node->outarcs;

      while ( arc != NULL )
      {
         if ( arc->type == CS )
            foundCS = TRUE;
         arc = arc->next_outarc;
      }

      /* only continue if a CS is found */
      if ( sum == 2 && Node->type == ENTRY && foundCS )
      {
         GAS_Arc* arc1 = NULL;
         GAS_Arc* arc2 = NULL;
         SCIP_Bool secondRun = FALSE;

         arc = Node->outarcs;
         /* set arc1 and arc2 */

         while ( arc != NULL )
         {
            if ( ! secondRun )
            {
               arc1 = arc;
               secondRun = TRUE;
               arc = arc->next_outarc;
            }
            else if ( secondRun )
            {
               arc2 = arc;
               arc = arc->next_outarc;
            }
         }

         assert( Node->nodeMixingRatio != NULL );
         assert( arc->mixingRatioMass != NULL );

         /* Generate constraints to fix mu on Node, CS and Bypass */
         if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
         {
            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 1, i.e., fix mu = 1\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 1.0) );

            /* Second the arc_mixingRatioMass can be fixed */
            assert( arc1->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatioMass, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatioMass, 1.0) );
            assert( arc2->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatioMass, 1.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatioMass, 1.0) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               assert( arc1->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatio, 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatio, 1.0) );
               assert( arc2->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatio, 1.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatio, 1.0) );
            }
         }
         else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2 ) )
         {

            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 0, i.e., fix mu = 0\n");
            /* First the node mixing ratio can be fixed */
            SCIPdebugMsg(scip, "Constituent 1, i.e., fix mu = 1\n");
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, 0.0) );

            /* Second the arc_mixingRatioMass can be fixed */
            assert( arc1->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatioMass, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatioMass, 0.0) );
            assert( arc2->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatioMass, 0.0) );
            SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatioMass, 0.0) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               assert( arc1->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatio, 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatio, 0.0) );
               assert( arc2->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatio, 0.0) );
               SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatio, 0.0) );
            }
         }
         /* otherwise the gas is a mixture and we need to calculate mu */
         else
         {
            SCIP_Real mu;
            SCIP_Real muMol;
            SCIP_Real M1 = 1.0 / probdata->molarMass1;
            SCIP_Real M2 = 1.0 / probdata->molarMass2;

            /* calculate the mass% mu */
            mu = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);
            /* calculate mol% mu */
            muMol = (mu * M1) / (mu * (M1 - M2) + M2);
            SCIPdebugMsg(scip, "Mixed Gas, molarMass = %f, mu = %f, muMol = %f\n", Node->molarMass, mu, muMol );

            /* First the node mixing ratio can be fixed */
            SCIP_CALL( SCIPchgVarLb(scip, Node->nodeMixingRatio, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, Node->nodeMixingRatio, mu) );

            /* Second the arc_mixingRatioMass can be fixed */
            assert( arc1->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatioMass, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatioMass, mu) );
            assert( arc2->mixingRatioMass != NULL );
            SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatioMass, mu) );
            SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatioMass, mu) );

            /* Third, if molar% is also used then arcMixingRatio can be fixed */
            if ( ! probdata->linearMixing )
            {
               assert( arc1->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc1->mixingRatio, muMol) );
               SCIP_CALL( SCIPchgVarUb(scip, arc1->mixingRatio, muMol) );
               assert( arc2->mixingRatio != NULL );
               SCIP_CALL( SCIPchgVarLb(scip, arc2->mixingRatio, muMol) );
               SCIP_CALL( SCIPchgVarUb(scip, arc2->mixingRatio, muMol) );
            }
         }
      }
      /* if sum of in and outarcs = 2 and node is no entry OR exit, we can fix mu and the nodeMixingRatio */
      else if ( sum == 2 && Node->type == INNODE )
      {
         SCIP_CONS* fixMuCons;
         SCIP_CONS* fixMuMassCons;
         SCIP_CONS* fixMuNodeCons;
         SCIP_VAR* vars[2];
         SCIP_Bool secondRun = FALSE;
         SCIP_Real vals[2];
         GAS_Arc* arc1 = NULL;
         GAS_Arc* arc2 = NULL;
         GAS_Arc* arcs[2];

         arcs[0] = Node->outarcs;
         arcs[1] = Node->inarcs;
         vals[0] = 1;
         vals[1] = -1;

         /* set arc1 and arc2 */
         for (int i = 0; i < 2; i++)
         {
            arc = arcs[i];
            while ( arc != NULL )
            {
               /* Check if outarc */
               if ( i == 0 && !secondRun )
               {
                  arc1 = arc;
                  secondRun = TRUE;
                  arc = arc->next_outarc;
               }
               else if ( i == 0 && secondRun )
               {
                  arc2 = arc;
                  arc = arc->next_outarc;
               }

               /* Check if inarc */
               if ( i == 1 && !secondRun )
               {
                  arc1 = arc;
                  secondRun = TRUE;
                  arc = arc->next_inarc;
               }
               else if ( i == 1 && secondRun )
               {
                  arc2 = arc;
                  arc = arc->next_inarc;
               }
            }
         }

         /* Now that arc1 and arc2 are set, we generate the contraint */
         /* mu_arc1 = mu_arc2 */
         if ( ! probdata->linearMixing )
         {
            assert( arc1 != NULL );
            assert( arc2 != NULL );
            assert( arc1->mixingRatio != NULL );
            assert( arc2->mixingRatio != NULL );

            vars[0] = arc1->mixingRatio;
            vars[1] = arc2->mixingRatio;

            /* Create the first constraint for Mu and add it to SCIP */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1A-2Arcs-fixMuCons#%s", Node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuCons, consname, 2, vars, vals, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, fixMuCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &fixMuCons) );
         }

         assert( arc1->mixingRatioMass != NULL );
         assert( arc2->mixingRatioMass != NULL );

         vars[0] = arc1->mixingRatioMass;
         vars[1] = arc2->mixingRatioMass;

         /* Create the second constraint for MuMass and add it to SCIP */
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1A-2Arcs-fixMuMassCons#%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuMassCons, consname, 2, vars, vals, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, fixMuMassCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &fixMuMassCons) );

         /* fix mu node aka nodeMixingRatio */
         assert( Node->nodeMixingRatio != NULL );
         vars[0] = Node->nodeMixingRatio;
         vars[1] = arc1->mixingRatioMass; /* does not matter if arc1 or arc2 is used */
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1N-2Arcs-fixNodeMu#%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuNodeCons, consname, 2, vars, vals, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, fixMuNodeCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &fixMuNodeCons) );
      }
      /* node has 2 or more neighbours and mu cannot be fixed */
      else if ( sum >= 2 )
      {
         SCIP_EXPR* exprQZMuprod[2];
         SCIP_EXPR** exprQZMuprodResult = NULL;
         SCIP_EXPR* exprSUMinMu = NULL;
         SCIP_EXPR** exprQZMuprodResultO = NULL;
         SCIP_EXPR* exprSUMoutMu = NULL;
         SCIP_EXPR* exprZsum = NULL;
         SCIP_EXPR* exprZarr[2];
         SCIP_Real coefPlusOneMinusOne[2];

         /* used in the loops below as offset i.e if the node is an entry the inflow is added to the flow sums */
         SCIP_Real inflow = 0.0;
         SCIP_Real inflowMU = 0.0;
         SCIP_Real copyInflow = 0.0;
         SCIP_Real copyInflowMU = 0.0;

         coefPlusOneMinusOne[0] = 1;
         coefPlusOneMinusOne[1] = -1;

         /* Adding the inflow at an as entry offset is complicated and happens as follows: if either ONLY inarcs or
          * outarcs exist, the offset is added in one of the next two if loops.*/
         /* Otherwise the offset is added in the 3rd if loop (i.e. if numInarcs && numoutarcs >0). */
         /* The complication is cauesd by the fact that the flow value has to be added exactly once. */
         /* The flow is always added in the denominator but only in the numerator if mu=1. */
         if ( Node->type == ENTRY )
         {
            if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
            {
               /* heavy gas => mu=1 */
               inflow = Node->flow;
               inflowMU = 1;
            }
            else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
            {
               /* light gas => mu=0 */
               inflow = Node->flow;
               inflowMU = 0;
            }
            /* gas is a mixture and we first need to calculate mu */
            else
            {
               /* calculate the mass% mu */
               inflowMU = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);
               inflow = Node->flow;
               SCIPdebugMsg(scip, "inflow on node with >=2 arcs and mu= %f, molarMass = %f, node =%s\n", inflowMU, Node->molarMass, Node->id );
            }
         }

         if ( numInarcs > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResult, numInarcs) );

            /* The loop below builds the pairwise produc q_a*mu_a for the inarcs. */
            /* NOTE: MASS% mixing Ratio is used i.e. mixingRatioMass */

            /* reset counter */
            c = 0;

            /* fill in incoming arcs */
            arc = Node->inarcs;
            while ( arc != NULL )
            {
               SCIP_EXPR* exprQaIn = NULL;
               SCIP_EXPR* exprMUaIn = NULL;

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaIn, arc->flowvar, NULL, NULL) );
               assert( arc->mixingRatioMass != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaIn, arc->mixingRatioMass, NULL, NULL) );

               exprQZMuprod[0] = exprQaIn;
               exprQZMuprod[1] = exprMUaIn;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResult[c], 2, exprQZMuprod, 1.0, NULL, NULL) ); /* Pairwise q_a*mu_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaIn) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaIn) );

               arc = arc->next_inarc;
               c++;
            }

            assert( c == numInarcs );
            c = 0;

            /* need to add the inflow value as offset - see details above */
            if ( numOutarcs == 0 )
            {
               copyInflow = inflow;
               copyInflowMU = inflowMU;
            }

            /* now the sum over the above created producs are build */
            /* Sum over inarcs q_a*mu_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMinMu, numInarcs, exprQZMuprodResult, NULL, copyInflowMU * copyInflow, NULL, NULL) );

            /*free arrays*/
            for (int i = 0; i < numInarcs; i++)
            {
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZMuprodResult[i]) );
            }

            SCIPfreeBufferArray(scip, &exprQZMuprodResult);
         }

         if ( numOutarcs > 0 )
         {
            /* Analog calculation now for the out arcs. */
            /* The loop below builds the pairwise producs q_a*z_a and q_a*z_a*mu_a for the outarcs. */
            /* NOTE: MASS% mixing Ratio is used i.e. mixingRatioMass. */
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResultO, numOutarcs) );

            c = 0;
            arc = Node->outarcs;

            /* fill in outgoing arcs */
            while ( arc != NULL )
            {
               SCIP_EXPR* exprQaOut = NULL;
               SCIP_EXPR* exprMUaOut = NULL;

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaOut, arc->flowvar, NULL, NULL) );
               assert( arc->mixingRatioMass != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaOut, arc->mixingRatioMass, NULL, NULL) );

               exprQZMuprod[0] = exprQaOut;
               exprQZMuprod[1] = exprMUaOut;

               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResultO[c], 2, exprQZMuprod, 1.0, NULL, NULL) ); /* pairwise q_a * mu_a product */
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQaOut) );
               SCIP_CALL( SCIPreleaseExpr(scip, &exprMUaOut) );

               arc = arc->next_outarc;
               c++;
            }

            assert( c == numOutarcs );

            /* In the case that the node has only inflow on outarcs i.e. the flow is negative we multiply it by -1. */
            SCIP_Real minusOne[numOutarcs];

            if ( numInarcs == 0 )
            {
               copyInflow = inflow;
               copyInflowMU = inflowMU;

               for (int i = 0; i < numOutarcs; i++)
                  minusOne[i] = -1.0;
            }
            else
            {
               for (int i = 0; i < numOutarcs; i++)
                  minusOne[i] = 1.0;
            }

            /* now the sums over the above created producs are built */
            /* Sum over inarcs q_a*z_a*mu_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMoutMu, numOutarcs, exprQZMuprodResultO, minusOne, copyInflowMU * copyInflow, NULL, NULL) );

            /*free arrays*/
            for (int i = 0; i < numOutarcs; i++)
            {
               SCIP_CALL( SCIPreleaseExpr(scip, &exprQZMuprodResultO[i]) );
            }
            SCIPfreeBufferArray(scip, &exprQZMuprodResultO);
         }

         /* In the following the mixing flow conservation equation is built. */
         if ( numInarcs > 0 && numOutarcs > 0 )
         {
            if ( Node->type != EXIT )
            {
               exprZarr[0] = exprSUMinMu;
               exprZarr[1] = exprSUMoutMu;
               SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 2, exprZarr, coefPlusOneMinusOne, inflowMU *inflow, NULL, NULL) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprZsum, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );
            }
            else
            {
               SCIP_EXPR* nodeMu;
               SCIP_EXPR* exprSum3[3];
               SCIP_Real coef3[3];
               coef3[0] = 1.0;
               coef3[1] = -1.0;
               coef3[2] = Node->flow;
               SCIP_CALL( SCIPcreateExprVar(scip, &nodeMu, Node->nodeMixingRatio, NULL, NULL) );

               exprSum3[0] = exprSUMinMu;
               exprSum3[1] = exprSUMoutMu;
               exprSum3[2] = nodeMu;
               SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 3, exprSum3, coef3, 0.0, NULL, NULL) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprZsum, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );

               SCIP_CALL( SCIPreleaseExpr(scip, &nodeMu) );
            }


         }
         else if ( numInarcs > 0 && numOutarcs == 0)
         {

            if ( Node->type != EXIT )
            {
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprSUMinMu, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );
            }
            else
            {
               SCIP_EXPR* nodeMu;
               SCIP_EXPR* exprSum3[3];
               SCIP_Real coef3[3];
               coef3[0] = 1.0;
               coef3[1] = Node->flow;
               SCIP_CALL( SCIPcreateExprVar(scip, &nodeMu, Node->nodeMixingRatio, NULL, NULL) );

               exprSum3[0] = exprSUMinMu;
               exprSum3[1] = nodeMu;
               SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 2, exprSum3, coef3, 0.0, NULL, NULL) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprZsum, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );

               SCIP_CALL( SCIPreleaseExpr(scip, &nodeMu) );
            }
         }
         else if ( numInarcs == 0 && numOutarcs > 0 )
         {
            if ( Node->type != EXIT )
            {
               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprSUMoutMu, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );
            }
            else
            {
               SCIP_EXPR* nodeMu;
               SCIP_EXPR* exprSum3[3];
               SCIP_Real coef3[3];
               coef3[0] = 1.0;
               coef3[1] = Node->flow;
               SCIP_CALL( SCIPcreateExprVar(scip, &nodeMu, Node->nodeMixingRatio, NULL, NULL) );

               exprSum3[0] = exprSUMoutMu;
               exprSum3[1] = nodeMu;
               SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 2, exprSum3, coef3, 0.0, NULL, NULL) );

               (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-MixtureFlowCons#%s", Node->id);
               SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprZsum, 0.0, 0.0) );
               SCIP_CALL( SCIPaddCons(scip, mixingCons) );
               SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );

               SCIP_CALL( SCIPreleaseExpr(scip, &nodeMu) );
            }
         }

         /* free all expr */
         if ( numInarcs > 0 && numOutarcs > 0 )
         {
            /*release expr*/
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMinMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMoutMu) );
            SCIP_CALL( SCIPreleaseExpr(scip, &exprZsum) );
         }
         else if ( numInarcs > 0 && numOutarcs == 0 )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMinMu) );
            if( Node->type == EXIT )
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZsum) );
         }
         else if ( numInarcs == 0 && numOutarcs > 0 )
         {
            SCIP_CALL( SCIPreleaseExpr(scip, &exprSUMoutMu) );
            if( Node->type == EXIT )
               SCIP_CALL( SCIPreleaseExpr(scip, &exprZsum) );
         }

         {
            /* Now we iterate over all arcs */
            /* We set  binvar(nodeMixingRatio - mixingRatioMass (Arc))=0. */
            /* Note that each constraint is generated twice since we do not know the flow direction in advance. */
            /* Coupling with the flowBinVars ensures that only the correct is active. */
            GAS_Arc* arcs[2];
            SCIP_CONS* linkConsMuMass;

            coefPlusOneMinusOne[0] = 1;
            coefPlusOneMinusOne[1] = -1;

            arcs[0] = Node->outarcs;
            arcs[1] = Node->inarcs;

            c = 0; /* counter */

            for (int i = 0; i < 2; i++)
            {
               arc = arcs[i];

               while ( arc != NULL )
               {
                  SCIP_EXPR* exprSum[2];
                  SCIP_EXPR* exprZ;
                  SCIP_EXPR* exprMuMassArc;
                  SCIP_EXPR* exprDiff;
                  SCIP_EXPR* exprDiffZ;
                  SCIP_EXPR* exprMixRatioNode;


                  /* Use binary variable in order to determine the arcs with outflow. */
                  if ( i == 0 )
                  {
                     /* flow on an outarc really flows out */
                     assert( arc->positiveFlowBinvar != NULL );
                     SCIP_CALL( SCIPcreateExprVar(scip, &exprZ, arc->positiveFlowBinvar, NULL, NULL) );
                  }
                  else
                  {
                     /* flow on an inarc flows out */
                     assert( arc->negativeFlowBinvar != NULL );
                     SCIP_CALL( SCIPcreateExprVar(scip, &exprZ, arc->negativeFlowBinvar, NULL, NULL) );
                  }

                  assert( Node->nodeMixingRatio != NULL );
                  SCIP_CALL( SCIPcreateExprVar(scip, &exprMixRatioNode, Node->nodeMixingRatio, NULL, NULL) );

                  /* Here we link nodeMixRatio to arcMixingRatioMass in the following way z_a * (nodeMu - arcMu) = 0. */
                  assert( arc->mixingRatioMass != NULL );
                  SCIP_CALL( SCIPcreateExprVar(scip, &exprMuMassArc, arc->mixingRatioMass, NULL, NULL) );
                  exprSum[0] = exprMixRatioNode;
                  exprSum[1] = exprMuMassArc;
                  SCIP_CALL( SCIPcreateExprSum(scip, &exprDiff, 2, exprSum, coefPlusOneMinusOne, 0.0, NULL, NULL) );
                  exprSum[0] = exprDiff;
                  exprSum[1] = exprZ;
                  SCIP_CALL( SCIPcreateExprProduct(scip, &exprDiffZ, 2, exprSum, 1.0, NULL, NULL) );

                  (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "TYPE-7Mass-linkNodeMixRatioToArcMixRatio%s#%s", arc->id, Node->id );
                  SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &linkConsMuMass, consname, exprDiffZ, 0.0, 0.0) );
                  SCIP_CALL( SCIPaddCons(scip, linkConsMuMass) );
                  SCIP_CALL( SCIPreleaseCons(scip, &linkConsMuMass) );

                  /* free all expressions used to form root expression */
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprZ) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprMixRatioNode) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprMuMassArc) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprDiff) );
                  SCIP_CALL( SCIPreleaseExpr(scip, &exprDiffZ) );

                  /* increment counter */
                  c++;

                  if ( i == 0 )
                     arc = arc->next_outarc;
                  else if ( i == 1 )
                     arc = arc->next_inarc;
               }
            }

            c = 0; /* reset counter */
         }
      }
   }
#endif /* SCIP_VERSION >= 800 */

   return SCIP_OKAY;
}

/** The following function is a variant which directly calculates the mol % mixing ratio and does not use the auxilary variables mu_mass
 *
 *  This makes the calculation of the mixing ratio much more involved since we create fraction of flow and  the molar mass of the mixture.
 */
static
SCIP_RETCODE generateMixingRatioMol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
#if SCIP_VERSION >= 800
   char consname[SCIP_MAXSTRLEN];
   GAS_Network* Network = probdata->network;
   SCIP_CONS* mixingCons;
   int j;
   int c = 0;
   int sum = 0;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   for (j = 0; j < Network->numnodes; j++)
   {
      GAS_Node* Node = &(Network->nodes_ptr[j]);
      GAS_Arc* arc;
      int numInarcs = 0;
      int numOutarcs = 0;

      assert( Node != NULL );

      /* Determine the number of in and outarcs */
      arc = Node->inarcs;
      while ( arc != NULL )
      {
         numInarcs++;
         arc = arc->next_inarc;
      }

      arc = Node->outarcs;
      while ( arc != NULL )
      {
         numOutarcs++;
         arc = arc->next_outarc;
      }

      sum = numInarcs + numOutarcs;

      if ( sum == 1 && Node->type == ENTRY )
      {
         SCIP_CONS* fixNodeMuCons;
         SCIP_CONS* fixMuCons;
         SCIP_VAR* vars[2];
         SCIP_VAR* var[1];
         SCIP_Real vals[2];
         SCIP_Real val[1];

         vals[0] = 1.0;
         vals[1] = -1.0;
         val[0] = 1.0;

         var[0] = Node->nodeMixingRatio;

         if ( Node->inarcs != NULL )
            arc = Node->inarcs;
         else
            arc = Node->outarcs;

         assert( Node->nodeMixingRatio != NULL );
         assert( arc->mixingRatio != NULL );

         if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
         {
            SCIPdebugMsg(scip, "Heavy Gas i.e. fix mu = 1\n");
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1N-fixNodeMixRatio%s", Node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixNodeMuCons, consname, 1, var, val, 1.0, 1.0) );
            SCIP_CALL( SCIPaddCons(scip, fixNodeMuCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &fixNodeMuCons) );

            vars[0] = Node->nodeMixingRatio;
            vars[1] = arc->mixingRatio;

            /* Create the constraint (nodeMixingRatio = arcMixingRatio) and add it to SCIP */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1A-fixArcMixRatioToNodeMixRatio%s", arc->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuCons, consname, 2, vars, vals, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, fixMuCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &fixMuCons) );
         }
         else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
         {
            SCIPdebugMsg(scip, "Light Gas i.e. fix mu=0\n");
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1N-fixNodeMixRatio%s", Node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixNodeMuCons, consname, 1, var, val, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, fixNodeMuCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &fixNodeMuCons) );

            vars[0] = Node->nodeMixingRatio;
            vars[1] = arc->mixingRatio;

            /* Create the constraint (nodeMixingRatio = arcMixingRatio) and add it to SCIP */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1A-fixArcMixRatioToNodeMixRatio%s", arc->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuCons, consname, 2, vars, vals, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, fixMuCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &fixMuCons) );
         }
         /* otherwise the gas is a mixture and we need to calculate mu */
         else
         {
            SCIP_Real mu;
            SCIP_Real muMol;
            SCIP_Real M1 = 1.0 / probdata->molarMass1;
            SCIP_Real M2 = 1.0 / probdata->molarMass2;

            /* calculate the mass% mu */
            mu = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);

            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1N-fixNodeMixRatio%s", Node->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixNodeMuCons, consname, 1, var, val, mu, mu) );
            SCIP_CALL( SCIPaddCons(scip, fixNodeMuCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &fixNodeMuCons) );

            /* calculate mol% mu */
            muMol = (mu * M1) / (mu * (M1 - M2) + M2);
            SCIPdebugMsg(scip, "Mixed Gas, molarMass = %f, mu = %f, muMol = %f\n", Node->molarMass, mu, muMol);

            var[0] = arc->mixingRatio;

            /* Create the constraint (arcMixingRatio = muMol ) and add it to SCIP */
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1A-fixArcMixRatioToNodeMixRatio%s", arc->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuCons, consname, 2, var, val, muMol, muMol) );
            SCIP_CALL( SCIPaddCons(scip, fixMuCons) );
            SCIP_CALL( SCIPreleaseCons(scip, &fixMuCons) );
         }
      }
      /* if sum of in and outarcs = 2 and node is no entry OR exit, we can fix mu and the nodeMixingRatio */
      else if ( sum == 2 && Node->type == INNODE )
      {
         SCIP_CONS* fixMuCons;
         SCIP_CONS* fixMuNodeCons;
         SCIP_VAR* vars[2];
         SCIP_Bool secondRun = FALSE;
         SCIP_Real vals[2];
         GAS_Arc* arc1 = NULL;
         GAS_Arc* arc2 = NULL;
         GAS_Arc* arcs[2];

         arcs[0] = Node->outarcs;
         arcs[1] = Node->inarcs;
         vals[0] = 1;
         vals[1] = -1;

         /* set arc1 and arc2 */
         for (int i = 0; i < 2; i++)
         {
            arc = arcs[i];
            while ( arc != NULL )
            {
               /* Check if outarc */
               if ( i == 0 && !secondRun )
               {
                  arc1 = arc;
                  secondRun = TRUE;
                  arc = arc->next_outarc;
               }
               else if ( i == 0 && secondRun )
               {
                  arc2 = arc;
                  arc = arc->next_outarc;
               }

               /* Check if inarc */
               if ( i == 1 && !secondRun )
               {
                  arc1 = arc;
                  secondRun = TRUE;
                  arc = arc->next_inarc;
               }
               else if ( i == 1 && secondRun )
               {
                  arc2 = arc;
                  arc = arc->next_inarc;
               }
            }
         }

         /* Now that arc1 and arc2 are set, we generate the contraint */
         /* mu_arc1 = mu_arc2 */
         assert( arc1 != NULL );
         assert( arc2 != NULL );
         assert( arc1->mixingRatio != NULL );
         assert( arc2->mixingRatio != NULL );

         vars[0] = arc1->mixingRatio;
         vars[1] = arc2->mixingRatio;

         /* Create the first constraint for Mu and add it to SCIP */
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1-2Arcs-fixMuCons%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuCons, consname, 2, vars, vals, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, fixMuCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &fixMuCons) );

         /* here MOL% is used */
         /* fix mu node aka nodeMixingRatio */
         assert( Node->nodeMixingRatio != NULL );
         vars[0] = Node->nodeMixingRatio;
         vars[1] = arc1->mixingRatio; /* does not matter if arc1 or arc2 is used */
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-1N-2Arcs-fixNodeMu%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &fixMuNodeCons, consname, 2, vars, vals, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, fixMuNodeCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &fixMuNodeCons) );
      }
      /* node has 2 or more neighbours and mu cannot be fixed */
      else if ( sum >= 2 )
      {
         SCIP_EXPR* exprQZprod[3];
         SCIP_EXPR** exprQZprodResult = NULL;
         SCIP_EXPR* exprQZMuprod[2];
         SCIP_EXPR** exprQZMuprodResult = NULL;
         SCIP_EXPR* exprSUMin = NULL;
         SCIP_EXPR* exprSUMinMu = NULL;
         SCIP_EXPR** exprQZprodResultO = NULL;
         SCIP_EXPR** exprQZMuprodResultO = NULL;
         SCIP_EXPR* exprSUMout = NULL;
         SCIP_EXPR* exprSUMoutMu = NULL;
         SCIP_EXPR* exprZsum = NULL;
         SCIP_EXPR* exprZarr[2];
         SCIP_Real coefPlusOneMinusOne[2];
         SCIP_EXPR* exprNsum = NULL;
         SCIP_EXPR* exprN[2];
         SCIP_EXPR* exprMuNode = NULL;
         SCIP_EXPR* exprProdN[2];
         SCIP_EXPR* exprMuSumProd = NULL;
         SCIP_EXPR* exprFsumm[2];
         SCIP_EXPR* exprFinMu;
         /* used in the loops below as offset i.e if the node is an entry the inflow is added to the flow sums */
         SCIP_Real inflow = 0.0;
         SCIP_Real inflowMU = 0.0;
         SCIP_Real muMass = 0.0;
         SCIP_Real copyInflow = 0.0;
         SCIP_Real copyInflowMU = 0.0;

         coefPlusOneMinusOne[0] = 1;
         coefPlusOneMinusOne[1] = -1;

         /* Adding the inflow at an as entry offset is complicated and happens as follows: if either ONLY inarcs or
          * outarcs exist, the offset is added in one of the next two if loops.*/
         /* Otherwise the offset is added in the 3rd if loop (i.e. if numInarcs && numoutarcs >0). */
         /* The complication is cauesd by the fact that the flow value has to be added exactly once. */
         /* The flow is always added in the denominator but only in the numerator if mu=1. */
         if ( Node->type == ENTRY )
         {
            if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass1) )
            {
               /* heavy gas => mu=1 */
               inflow = Node->flow / probdata->molarMass1;
               inflowMU = 1;
            }
            else if ( SCIPisEQ(scip, Node->molarMass, probdata->molarMass2) )
            {
               /* light gas => mu=0 */
               inflow = Node->flow / probdata->molarMass2;
               inflowMU = 0;
            }
            /* gas is a mixture and we first need to calculate mu */
            else
            {
               SCIP_Real M1 = 1.0 / probdata->molarMass1;
               SCIP_Real M2 = 1.0 / probdata->molarMass2;
               SCIP_Real MofMu;

               /* calculate the mass% mu */
               muMass = (Node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);

               /* calculate mol% mu */
               inflowMU = (muMass * M1) / (muMass * (M1 - M2) + M2);
               MofMu = probdata->molarMass2 + (probdata->molarMass1 - probdata->molarMass2) * inflowMU;
               inflow = Node->flow / MofMu ; /*need to calculate the number of moles*/
               SCIPdebugMsg(scip, "inflow on node with >=2 arcs and mu= %f, molarMass = %f, node =%s\n", inflowMU, Node->molarMass, Node->id);
            }
         }

         if ( numInarcs > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZprodResult, numInarcs) );
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResult, numInarcs) );

            /* The loop below builds the pairwise producs q_a*z_a and q_a*z_a*mu_a for the inarcs. */
            /* NOTE: MASS% mixing Ratio is used i.e. mixingRatioMass */

            /* reset counter */
            c = 0;

            /* fill in incoming arcs */
            arc = Node->inarcs;
            while ( arc != NULL )
            {
               SCIP_EXPR* exprQaIn = NULL;
               SCIP_EXPR* exprZaIn = NULL;
               SCIP_EXPR* exprMUaIn = NULL;
               SCIP_EXPR* exprMolarMass = NULL;
               SCIP_EXPR* exprInv = NULL;
               SCIP_Real coef = (probdata->molarMass1 - probdata->molarMass2);

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaIn, arc->flowvar, NULL, NULL) );

               assert( arc->positiveFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaIn, arc->positiveFlowBinvar, NULL, NULL) );

               assert( arc->mixingRatio != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaIn, arc->mixingRatio, NULL, NULL) );

               SCIP_CALL( SCIPcreateExprSum(scip, &exprMolarMass, 1, &exprMUaIn, &coef, probdata->molarMass2, NULL,NULL) );
               SCIP_CALL( SCIPcreateExprPow(scip, &exprInv, exprMolarMass, -1.0, NULL, NULL) );

               exprQZprod[0] = exprQaIn;
               exprQZprod[1] = exprZaIn;
               exprQZprod[2] = exprInv;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZprodResult[c], 3, exprQZprod, 1.0, NULL, NULL) ); /* pairwise q_a * z_a / M(mu) product */

               exprQZMuprod[0] = exprQZprodResult[c];
               exprQZMuprod[1] = exprMUaIn;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResult[c], 2, exprQZMuprod, 1.0, NULL, NULL) ); /* pairwise q_a * z_a * mu_a / M(mu) product */

               arc = arc->next_inarc;
               c++;
            }

            assert( c == numInarcs );
            c = 0;

            /* need to add the inflow value as offset - see details above */
            if ( numOutarcs == 0 )
            {
               copyInflow = inflow;
               copyInflowMU = inflowMU;
            }

            /* now the sums over the above created producs are build */
            /*the first sum calculates the total number of moles of inflow*/
            /* Sum over inarcs q_a * z_a / M(mu) prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMin, numInarcs, exprQZprodResult, NULL, copyInflow, NULL, NULL) );

            /* Sum over inarcs q_a * z_a * mu_a / M(mu) prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMinMu, numInarcs, exprQZMuprodResult, NULL, copyInflowMU * copyInflow, NULL, NULL) );

            SCIPfreeBufferArray(scip, &exprQZMuprodResult);
            SCIPfreeBufferArray(scip, &exprQZprodResult);
         }

         if ( numOutarcs > 0 )
         {
            /* Analog calculation now for the out arcs. */
            /* The loop below builds the pairwise producs q_a*z_a and q_a*z_a*mu_a for the outarcs. */
            /* NOTE: MASS% mixing Ratio is used i.e. mixingRatioMass. */
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZprodResultO, numOutarcs) );
            SCIP_CALL( SCIPallocBufferArray(scip, &exprQZMuprodResultO, numOutarcs) );

            c = 0;
            arc = Node->outarcs;

            /* fill in outgoing arcs */
            while ( arc != NULL )
            {
               SCIP_EXPR* exprQaOut = NULL;
               SCIP_EXPR* exprZaOut = NULL;
               SCIP_EXPR* exprMUaOut = NULL;
               SCIP_EXPR* exprMolarMass = NULL;
               SCIP_EXPR* exprInv = NULL;
               SCIP_Real coef = (probdata->molarMass1 - probdata->molarMass2);

               assert( arc->flowvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprQaOut, arc->flowvar, NULL, NULL) );

               assert( arc->negativeFlowBinvar != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprZaOut, arc->negativeFlowBinvar, NULL, NULL) );

               assert( arc->mixingRatio != NULL );
               SCIP_CALL( SCIPcreateExprVar(scip, &exprMUaOut, arc->mixingRatio, NULL, NULL) );

               SCIP_CALL( SCIPcreateExprSum(scip, &exprMolarMass, 1, &exprMUaOut, &coef, probdata->molarMass2, NULL, NULL) );
               SCIP_CALL( SCIPcreateExprPow(scip, &exprInv, exprMolarMass, -1.0, NULL, NULL) );

               exprQZprod[0] = exprQaOut;
               exprQZprod[1] = exprZaOut;
               exprQZprod[2] = exprInv;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZprodResultO[c], 3, exprQZprod, 1.0, NULL, NULL) ); /* pairwise q_a * z_a / M(mu) product */

               exprQZMuprod[0] = exprQZprodResultO[c];
               exprQZMuprod[1] = exprMUaOut;
               SCIP_CALL( SCIPcreateExprProduct(scip, &exprQZMuprodResultO[c], 2, exprQZMuprod, 1.0, NULL, NULL) ); /* pairwise q_a * z_a * mu_a product */

               arc = arc->next_outarc;
               c++;
            }

            assert( c == numOutarcs );

            if ( numInarcs == 0 )
            {
               copyInflow = inflow; /*note that this is the total number of moles entering the node*/
               copyInflowMU = inflowMU; /* number of moles of gas 1*/
            }

            /* now the sums over the above created producs are built */

            /* Sum over inarcs q_a*z_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMout, numOutarcs, exprQZprodResultO, NULL, copyInflow, NULL, NULL) );

            /* Sum over inarcs q_a*z_a*mu_a prod + inflow if needed. see above */
            SCIP_CALL( SCIPcreateExprSum(scip, &exprSUMoutMu, numOutarcs, exprQZMuprodResultO, NULL, copyInflowMU * copyInflow, NULL, NULL) );

            SCIPfreeBufferArray(scip, &exprQZMuprodResultO);
            SCIPfreeBufferArray(scip, &exprQZprodResultO);
         }

         /* In the following the mixing equation (Equation (7) of the mixing paper) is built. */
         if ( numInarcs > 0 && numOutarcs > 0 )
         {
            /* Now we first build the numerator for the mixing eqaution and then the denominator down below. */
            /* Note that the equation will be written in non fractional form (no div by 0). */

            /* Numerator sum of the mixing equation - EQ (7) in the paper. */
            exprZarr[0] = exprSUMinMu;
            exprZarr[1] = exprSUMoutMu;
            SCIP_CALL( SCIPcreateExprSum(scip, &exprZsum, 2, exprZarr, coefPlusOneMinusOne, inflowMU *inflow, NULL, NULL) );

            /* Mu_Node*NominatorSum - Denominator sum muliplied by mu */
            exprN[0] = exprSUMin;
            exprN[1] = exprSUMout;
            SCIP_CALL( SCIPcreateExprSum(scip, &exprNsum, 2, exprN, coefPlusOneMinusOne, inflow, NULL, NULL) );

            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );
            exprProdN[0] = exprMuNode;
            exprProdN[1] = exprNsum;
            SCIP_CALL( SCIPcreateExprProduct(scip, &exprMuSumProd, 2, exprProdN, 1.0, NULL, NULL) );

            exprFsumm[0] = exprMuSumProd;
            exprFsumm[1] = exprZsum;
         }
         else if ( numInarcs > 0 && numOutarcs == 0)
         {
            /* Denominator is only sum over inarcs = exprSUMinMu */

            /* Mu_Node * NominatorSum */
            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );
            exprProdN[0] = exprMuNode;
            exprProdN[1] = exprSUMin;
            SCIP_CALL( SCIPcreateExprProduct(scip, &exprMuSumProd, 2, exprProdN, 1.0, NULL, NULL) );

            exprFsumm[0] = exprMuSumProd;
            exprFsumm[1] = exprSUMinMu;
         }
         else if ( numInarcs == 0 && numOutarcs > 0 )
         {
            /* Denominator is only sum over outarcs = exprSUMoutMu */

            /* Mu_Node * NominatorSum */
            assert( Node->nodeMixingRatio != NULL );
            SCIP_CALL( SCIPcreateExprVar(scip, &exprMuNode, Node->nodeMixingRatio, NULL, NULL) );
            exprProdN[0] = exprMuNode;
            exprProdN[1] = exprSUMout;
            SCIP_CALL( SCIPcreateExprProduct(scip, &exprMuSumProd, 2, exprProdN, 1.0, NULL, NULL) );

            exprFsumm[0] = exprMuSumProd;
            exprFsumm[1] = exprSUMoutMu;
         }

         /* Final equation */
         /* Mu_Node*NominatorSum - DenominatorSum = 0 */
         SCIP_CALL( SCIPcreateExprSum(scip, &exprFinMu, 2, exprFsumm, coefPlusOneMinusOne, 0.0, NULL, NULL) );

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "TYPE-6-NodeMixingRatioCalculation_%s", Node->id);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &mixingCons, consname, exprFinMu, 0.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, mixingCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &mixingCons) );

         /* recursively free all expressions used to form root expression */
         SCIP_CALL( GASexprFreeRecurisve(scip, &exprFinMu) );

         {
            /* Now we iterate over all arcs and first calculate and set the molar percentage mixing ratio for each arc. */
            /* Afterwards we set nodeMixingRatio= mixingRatioMass (Arc). */
            /* Note that each constraint is generated twice since we do not know the flow direction in advance. */
            /* Coupling with the flowBinVars ensures that only the correct is active. */
            GAS_Arc* arcs[2];
            SCIP_CONS* linkConsMuMass;

            coefPlusOneMinusOne[0] = 1;
            coefPlusOneMinusOne[1] = -1;

            arcs[0] = Node->outarcs;
            arcs[1] = Node->inarcs;

            c = 0; /* counter */

            for (int i = 0; i < 2; i++)
            {
               arc = arcs[i];

               while ( arc != NULL )
               {
                  SCIP_EXPR* exprSum[2];
                  SCIP_EXPR* exprZ;
                  SCIP_EXPR* exprMuMassArc;
                  SCIP_EXPR* exprDiff;
                  SCIP_EXPR* exprDiffZ;
                  SCIP_EXPR* exprMixRatioNode;

                  /* Use binary variable in order to determine the arcs with outflow. */
                  if ( i == 0 )
                  {
                     /* flow on an outarc really flows out */
                     assert( arc->positiveFlowBinvar != NULL );
                     SCIP_CALL( SCIPcreateExprVar(scip, &exprZ, arc->positiveFlowBinvar, NULL, NULL) );
                  }
                  else
                  {
                     /* flow on an inarc flows out */
                     assert( arc->negativeFlowBinvar != NULL );
                     SCIP_CALL( SCIPcreateExprVar(scip, &exprZ, arc->negativeFlowBinvar, NULL, NULL) );
                  }

                  assert( Node->nodeMixingRatio != NULL );
                  SCIP_CALL( SCIPcreateExprVar(scip, &exprMixRatioNode, Node->nodeMixingRatio, NULL, NULL) );

                  /* Here we link nodeMixRatio to arcMixingRatio in the following way z_a * (nodeMu - arcMu) = 0. */
                  assert( arc->mixingRatio != NULL );
                  SCIP_CALL( SCIPcreateExprVar(scip, &exprMuMassArc, arc->mixingRatio, NULL, NULL) );
                  exprSum[0] = exprMixRatioNode;
                  exprSum[1] = exprMuMassArc;
                  SCIP_CALL( SCIPcreateExprSum(scip, &exprDiff, 2, exprSum, coefPlusOneMinusOne, 0.0, NULL, NULL) );
                  exprSum[0] = exprDiff;
                  exprSum[1] = exprZ;
                  SCIP_CALL( SCIPcreateExprProduct(scip, &exprDiffZ, 2, exprSum, 1.0, NULL, NULL) );

                  (void)SCIPsnprintf( consname, SCIP_MAXSTRLEN, "TYPE-7Mass-linkNodeMixRatioToArcMixRatio%s", arc->id );
                  SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &linkConsMuMass, consname, exprDiffZ, 0.0, 0.0) );
                  SCIP_CALL( SCIPaddCons(scip, linkConsMuMass) );
                  SCIP_CALL( SCIPreleaseCons(scip, &linkConsMuMass) );

                  /* recursively free all expressions used to form root expression */
                  SCIP_CALL( GASexprFreeRecurisve(scip, &exprDiffZ) );

                  /* increment counter */
                  c++;

                  if ( i == 0 )
                     arc = arc->next_outarc;
                  else if ( i == 1 )
                     arc = arc->next_inarc;
               }
            }

            c = 0; /* reset counter */
         }
      }
   }
#endif /* SCIP_VERSION >= 800 */

   return SCIP_OKAY;
}



/** generate constraint for shortpipes; there is no pressure loss, thus p_w = p_v */
static
SCIP_RETCODE generateShortPipeConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   arcposition         /**< position in arc array */
   )
{
   SCIP_CONS* ShortPipe_Cons;
   char consname[SCIP_MAXSTRLEN];
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   /* Name of the constraint */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "ShortPipeCons#%s#", probdata->network->arcs_ptr[arcposition].id);

   /* Array, that points to variables pressurevars */
   vars[0] = probdata->network->arcs_ptr[arcposition].targetnode->pressurevar;
   vars[1] = probdata->network->arcs_ptr[arcposition].sourcenode->pressurevar;

   /* Array, with coefficients (1 ,-1) */
   vals[0] =  1 ;
   vals[1] = -1;

   /* Create constraint and add it to scip */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &ShortPipe_Cons, consname, 2, vars, vals, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, ShortPipe_Cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &ShortPipe_Cons) );

   return SCIP_OKAY;
}

/** generate valve constraints
 *
 * The following function generates four constraints.
 *
 * valve_pressure_cons1: (pressureMax_v - pressureMin_u) * valve_binvar + pressurevar_v - pressurevar_u <= pressureMax_v - pressureMin_u
 * valve_pressure_cons2: (pressureMax_u - pressureMin_v) * valve_binvar + pressurevar_u - pressurevar_v <= pressureMax_u - pressureMin_v
 * valve_flow_cons1: flowMin(in arc) * valve_binvar <= flowvar
 * valve_flow_cons2: flowMax(in arc) * valve_binvar => flowvar
 */
static
SCIP_RETCODE generateValveConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   arcposition,        /**< position of arc */
   int                   valveposition       /**< position of valve */
   )
{
   SCIP_VAR* Sourcevariable;
   SCIP_VAR* Targetvariable;
   SCIP_Real TargetMaxMinusSourceMin;
   SCIP_Real SourceMaxMinusTargetMin;
   GAS_Node* Sourcenode;
   GAS_Node* Targetnode;
   GAS_Valve* valve;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3];
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* valve_pressure_cons;
   SCIP_CONS* valve_flow_cons;

   valve = probdata->network->arcs_ptr[arcposition].detailed_info;

   Sourcenode = probdata->network->arcs_ptr[arcposition].sourcenode;
   Targetnode = probdata->network->arcs_ptr[arcposition].targetnode;
   TargetMaxMinusSourceMin = MIN(Targetnode->pressureMax - Sourcenode->pressureMin, valve->pressureDifferentialMax);
   SourceMaxMinusTargetMin = MIN(Sourcenode->pressureMax - Targetnode->pressureMin, valve->pressureDifferentialMax);
   Sourcevariable = probdata->network->arcs_ptr[arcposition].sourcenode->pressurevar;
   Targetvariable = probdata->network->arcs_ptr[arcposition].targetnode->pressurevar;

   /*
    *  create constraints that link the pressure variables if the valve is open
    */
   /* Name of the constraint */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "VALVE_pressureCons1#%s#", probdata->network->arcs_ptr[arcposition].id);

   /*Array, that points to variables valve_binvar and pressurevars */
   vars[0] = probdata->VALVE_binvars[valveposition];
   vars[1] = Targetvariable;
   vars[2] = Sourcevariable;
   /* Array, with coefficients (TargetMaxMinusSourceMin, 1 ,-1) */
   vals[0] =  TargetMaxMinusSourceMin;
   vals[1] =  1;
   vals[2] = -1;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &valve_pressure_cons, consname, 3, vars, vals, -SCIPinfinity(scip) ,TargetMaxMinusSourceMin) );
   SCIP_CALL( SCIPaddCons(scip, valve_pressure_cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &valve_pressure_cons) );

   /* Array, with coefficients (SourceMaxMinusTargetMin,-1 ,1) */
   vals[0] =  SourceMaxMinusTargetMin;
   vals[1] = -1;
   vals[2] =  1;

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "VALVE_pressureCons2#%s#", probdata->network->arcs_ptr[arcposition].id);
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &valve_pressure_cons, consname, 3, vars , vals, -SCIPinfinity(scip) , SourceMaxMinusTargetMin) );
   SCIP_CALL( SCIPaddCons(scip, valve_pressure_cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &valve_pressure_cons) );

   /*
    *  create the link between binary and flow variable:
    *  q == 0.0 if binary variable is 0
    */
   /* Next we create the constraint valve_binvar * flowmin(in arc) <= flowvar <= valve_binvar * flowmax(in arc) using cons_varbound.h */
   if ( probdata->network->arcs_ptr[arcposition].flowMin != 0.0 )
   {
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "VALVE_flowCons1#%s#", probdata->network->arcs_ptr[arcposition].id);

      /* create the variable bound constraint and add it to SCIP */
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &valve_flow_cons, consname, probdata->flowvars[arcposition], probdata->VALVE_binvars[valveposition], -1 *(probdata->network->arcs_ptr[arcposition].flowMin), 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, valve_flow_cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &valve_flow_cons) );
   }

   if ( probdata->network->arcs_ptr[arcposition].flowMax != 0.0 )
   {
      /* create the variable bound constraint and add it to SCIP */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "VALVE_flowCons2#%s#", probdata->network->arcs_ptr[arcposition].id);

      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &valve_flow_cons, consname, probdata->flowvars[arcposition], probdata->VALVE_binvars[valveposition], -1 * (probdata->network->arcs_ptr[arcposition].flowMax), - SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCons(scip, valve_flow_cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &valve_flow_cons) );
   }

   return SCIP_OKAY;
}


/** generate control valve constraints
 *
 *  The following function generates four constraints.
 *   controlvalve_flow_cons1: flowMin(in arc) * controlvalve_binvar <= flowvar
 *   controlvalve_flow_cons2: flowMax(in arc) * controlvalve_binvar => flowvar
 *
 *   controlvalve_pressure_cons1: (pressureMax_v - pressureMin_u + pressureDifferentialMin) * controlvalve_binvar
 *    + pressurevar_v - pressurevar_u <= pressureMax_v - pressureMin_u
 *
 *   controlvalve_pressure_cons2: (pressureMax_u - pressureMin_v - pressureDifferentialMax) * controlvalve_binvar
 *    + pressurevar_u - pressurevar_v <= pressureMax_u - pressureMin_v
 *
 *  This function interprets all bounds given in the XML-file as refering to the control valve itself, excluding
 *  resistors.
 */
static
SCIP_RETCODE generateControlValveConstraintsExternalResistors(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GAS_Arc*              arc                 /**< controlValve arc */
   )
{
   GAS_Controlvalve* cv;
   GAS_Node*         sourcenode;
   GAS_Node*         targetnode;
   SCIP_VAR*         cvbinvar;
   SCIP_VAR*         sourcepressure;
   SCIP_VAR*         targetpressure;
   SCIP_VAR*         flowvar;
   SCIP_VAR*         vars[3];
   SCIP_Real         vals[3];
   SCIP_Real         lhs;
   SCIP_Real         rhs;
   SCIP_Real         pInMin;
   SCIP_Real         pOutMax;
   SCIP_Real         coeff;
   SCIP_CONS*        cons;
   char              consname[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( arc != NULL );

   cv                      = arc->detailed_info;
   sourcenode              = arc->sourcenode;
   targetnode              = arc->targetnode;
   cvbinvar                = cv->binvar;
   sourcepressure          = sourcenode->pressurevar;
   targetpressure          = targetnode->pressurevar;
   flowvar                 = arc->flowvar;

   pInMin  = cv->pressureInMin;
   pOutMax = cv->pressureOutMax;

   /* relax pressure bounds if option relaxCVbounds is used */
   if ( probdata->relaxCVbounds )
   {
      pInMin  = MAX(1.01325, pInMin - probdata->relaxCVvalue);
      pOutMax = pOutMax + probdata->relaxCVvalue;
   }

   /*
    *  there has to be no flow if cv is inactive
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_flowCons#%s#", arc->id);

   lhs   = - SCIPinfinity(scip);
   rhs   = 0.0;
   coeff = - arc->flowMax;

   /* add: flowvar - coeff * cvbinvar <= 0 */
   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, flowvar, cvbinvar, coeff, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* link control valve binvar to posflowbinvar */
   (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_link#%s", arc->id);
   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, arc->positiveFlowBinvar, cv->binvar, -1.0, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /*
    *  if cv is active then pressureInMin and pressureOutMax have to be satisfied
    */
   if ( SCIPisLT(scip, sourcenode->pressureMin, pInMin + cv->pressureLossIn) )
   {
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_minPinCons#%s#", arc->id);

      lhs   = sourcenode->pressureMin;
      rhs   = SCIPinfinity(scip);
      coeff = sourcenode->pressureMin - (pInMin + cv->pressureLossIn);
      /* coeff = - (cv->pressureInMin - sourcenode->pressureMin); */

      /* p_u + (p_u_min - (cv_p_min + cv_ploss_in)) * cvbinvar >= p_u_min */
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, sourcepressure, cvbinvar, coeff, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   if ( SCIPisGT(scip, targetnode->pressureMax, pOutMax - cv->pressureLossOut) )
   {
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_maxPoutCons#%s#", arc->id);

      lhs   = - SCIPinfinity(scip);
      rhs   = targetnode->pressureMax;
      coeff = targetnode->pressureMax - (pOutMax - cv->pressureLossOut);
      /* coeff = rhs - cv->pressureOutMax; */

      /* p_v + (p_v_max - (cv_p_max - cv_ploss_out)) * cvbinvar <= p_v_max */
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, targetpressure, cvbinvar, coeff, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /*  first constraint for pressure difference:
    *  p_u - p_v + (p_u_min - p_v_max - cvdiffmin) * cvbinvar >= (p_u_min - p_v_max)  <=>
    *  p_v - p_u + (p_v_max - p_u_min + cvdiffmin) * cvbinvar <= (p_v_max - p_u_min)
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_minPressureDiff#%s#", arc->id);

   lhs = sourcenode->pressureMin - targetnode->pressureMax;
   rhs = SCIPinfinity(scip);

   vars[0] = sourcepressure;
   vars[1] = targetpressure;
   vars[2] = cvbinvar;

   vals[0] = 1.0;
   vals[1] = - 1.0;
   vals[2] = lhs - (cv->pressureDifferentialMin + cv->pressureLossIn + cv->pressureLossOut);
   /* vals[2] = lhs - cv->pressureDifferentialMin; */

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 3, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /*  second constraint for pressure difference:
    *  p_u - p_v + (p_u_max - p_v_min - cvdiffmin) * cvbinvar <= (p_u_max - p_v_min)
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_maxPressureDiff#%s#", arc->id);

   lhs = - SCIPinfinity(scip);
   rhs = sourcenode->pressureMax - targetnode->pressureMin;

   vals[2] = rhs - (cv->pressureDifferentialMax + cv->pressureLossIn + cv->pressureLossOut);
   /* vals[2] = rhs - cv->pressureDifferentialMax; */

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 3, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /*
    *  constraints only necessary if there is a bypass
    */
   if ( cv->internalBypassRequired )
   {
      GAS_Valve* bypass;
      SCIP_VAR*  binvarBypass;

      bypass       = cv->bypass->detailed_info;
      binvarBypass = bypass->valve_binvar;

      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_binvarSum#%s#", arc->id);

      /* lhs = 0.0; */
      lhs = -SCIPinfinity(scip);
      rhs = 1.0;

      vars[0] = cvbinvar;
      vars[1] = binvarBypass;

      vals[0] = 1.0;
      vals[1] = 1.0;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 2, vars, vals, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      if ( SCIPisPositive(scip, arc->flowMin) )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_minFlowCons#%s#", arc->id);

         lhs = arc->flowMin;
         rhs = arc->flowMax;

         vars[0] = flowvar;
         vars[1] = cv->bypass->flowvar;

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 2, vars, vals, lhs, rhs) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   return SCIP_OKAY;
}


/** generate control valve constraints
 *
 *  The following function generates four constraints.
 *   controlvalve_flow_cons1: flowMin(in arc) * controlvalve_binvar <= flowvar
 *   controlvalve_flow_cons2: flowMax(in arc) * controlvalve_binvar => flowvar
 *
 *   controlvalve_pressure_cons1: (pressureMax_v - pressureMin_u + pressureDifferentialMin) * controlvalve_binvar
 *    + pressurevar_v - pressurevar_u <= pressureMax_v - pressureMin_u
 *
 *   controlvalve_pressure_cons2: (pressureMax_u - pressureMin_v - pressureDifferentialMax) * controlvalve_binvar
 *    + pressurevar_u - pressurevar_v <= pressureMax_u - pressureMin_v
 *
 *  This function interprets all bounds given in the XML-file as refering to the control valve station, including
 *  resistors.
 */
static
SCIP_RETCODE generateControlValveConstraintsInternalResistors(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GAS_Arc*              arc                 /**< controlValve arc */
   )
{
   GAS_Controlvalve* cv;
   GAS_Node*         sourcenode;
   GAS_Node*         targetnode;
   SCIP_VAR*         cvbinvar;
   SCIP_VAR*         sourcepressure;
   SCIP_VAR*         targetpressure;
   SCIP_VAR*         flowvar;
   SCIP_VAR*         vars[3];
   SCIP_Real         vals[3];
   SCIP_Real         lhs;
   SCIP_Real         rhs;
   SCIP_Real         pInMin;
   SCIP_Real         pOutMax;
   SCIP_Real         coeff;
   SCIP_CONS*        cons;
   char              consname[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( arc != NULL );

   cv                      = arc->detailed_info;
   sourcenode              = arc->sourcenode;
   targetnode              = arc->targetnode;
   cvbinvar                = cv->binvar;
   sourcepressure          = sourcenode->pressurevar;
   targetpressure          = targetnode->pressurevar;
   flowvar                 = arc->flowvar;

   pInMin  = cv->pressureInMin;
   pOutMax = cv->pressureOutMax;

   /* relax pressure bounds if option relaxCVbounds is used */
   if ( probdata->relaxCVbounds )
   {
      pInMin  = MAX(1.01325, pInMin - probdata->relaxCVvalue);
      pOutMax = pOutMax + probdata->relaxCVvalue;
   }

   /*
    *  there has to be no flow if cv is inactive
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_flowCons#%s#", arc->id);

   lhs   = - SCIPinfinity(scip);
   rhs   = 0.0;
   coeff = - arc->flowMax;

   /* add: flowvar - coeff * cvbinvar <= 0 */
   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, flowvar, cvbinvar, coeff, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* link control valve binvar to posflowbinvar */
   (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_link#%s", arc->id);
   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, arc->positiveFlowBinvar, cv->binvar, -1.0, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /*
    *  if cv is active then pressureInMin and pressureOutMax have to be satisfied
    */
   if ( SCIPisLT(scip, sourcenode->pressureMin, pInMin) )
   {
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_minPinCons#%s#", arc->id);

      lhs   = sourcenode->pressureMin;
      rhs   = SCIPinfinity(scip);
      coeff = sourcenode->pressureMin - pInMin;

      /* p_u + (p_u_min - cv_p_min) * cvbinvar >= p_u_min */
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, sourcepressure, cvbinvar, coeff, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   if ( SCIPisGT(scip, targetnode->pressureMax, pOutMax) )
   {
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_maxPoutCons#%s#", arc->id);

      lhs   = - SCIPinfinity(scip);
      rhs   = targetnode->pressureMax;
      coeff = targetnode->pressureMax - pOutMax;

      /* p_v + (p_v_max - cv_p_max) * cvbinvar <= p_v_max */
      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &cons, consname, targetpressure, cvbinvar, coeff, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /*  first constraint for pressure difference:
    *  p_u - p_v + (p_u_min - p_v_max - (cvdiffmin + pLossIn + pLossOut)) * cvbinvar >= (p_u_min - p_v_max)  <=>
    *  p_v - p_u + (p_v_max - p_u_min + (cvdiffmin + pLossIn + pLossOut)) * cvbinvar <= (p_v_max - p_u_min)
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_minPressureDiff#%s#", arc->id);

   lhs = sourcenode->pressureMin - targetnode->pressureMax;
   rhs = SCIPinfinity(scip);

   vars[0] = sourcepressure;
   vars[1] = targetpressure;
   vars[2] = cvbinvar;

   vals[0] = 1.0;
   vals[1] = - 1.0;
   vals[2] = lhs - (cv->pressureDifferentialMin + cv->pressureLossIn + cv->pressureLossOut);

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 3, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /*  second constraint for pressure difference:
    *  p_u - p_v + (p_u_max - p_v_min - (cvdiffmax + pLossIn + pLossOut)) * cvbinvar <= (p_u_max - p_v_min)
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_maxPressureDiff#%s#", arc->id);

   lhs = - SCIPinfinity(scip);
   rhs = sourcenode->pressureMax - targetnode->pressureMin;

   vals[2] = rhs - (cv->pressureDifferentialMax + cv->pressureLossIn + cv->pressureLossOut);

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 3, vars, vals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /*
    *  constraints only necessary if there is a bypass
    */
   if ( cv->internalBypassRequired )
   {
      GAS_Valve* bypass;
      SCIP_VAR*  binvarBypass;

      bypass       = cv->bypass->detailed_info;
      binvarBypass = bypass->valve_binvar;

      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_binvarSum#%s#", arc->id);

      lhs = -SCIPinfinity(scip);
      rhs = 1.0;

      vars[0] = cvbinvar;
      vars[1] = binvarBypass;

      vals[0] = 1.0;
      vals[1] = 1.0;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 2, vars, vals, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      if ( SCIPisPositive(scip, arc->flowMin) )
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CV_minFlowCons#%s#", arc->id);

         lhs = arc->flowMin;
         rhs = arc->flowMax;

         vars[0] = flowvar;
         vars[1] = cv->bypass->flowvar;

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 2, vars, vals, lhs, rhs) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   return SCIP_OKAY;
}


/** generate additional inequalities */
static
SCIP_RETCODE generateAdditionalFacets(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GAS_Arc*              arc                 /**< a compressor arc */
   )
{
   GAS_CS* compressor;
   GAS_CSBoxCons* boxcons;
   SCIP_VAR* sourcevariable;
   SCIP_VAR* targetvariable;
   SCIP_VAR* flowvariable;
   SCIP_VAR* vars[3];
   SCIP_CONS* facet_cons;
   char consname[SCIP_MAXSTRLEN];
   int i;
   int j;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   compressor = arc->detailed_info;

   for (j = 0; j < compressor->numconfigurations; ++j)
   {
      boxcons = &(compressor->boxcons[j]);
      sourcevariable = boxcons->BCM_pin;
      targetvariable = boxcons->BCM_pout;
      flowvariable = boxcons->BCM_flow;

      if ( strcmp(boxcons->a, "massFlow") == 0 )
      {
         assert( strcmp(boxcons->a_unit, "kg_per_s") == 0);
         vars[0] = flowvariable;
      }
      else if ( strcmp(boxcons->a, "pressureIn") == 0 )
      {
         assert( strcmp(boxcons->a_unit, "bar") == 0);
         vars[0] = sourcevariable;
      }
      else if ( strcmp(boxcons->a, "pressureOut") == 0 )
      {
         assert( strcmp(boxcons->a_unit, "bar") == 0);
         vars[0] = targetvariable;
      }
      if ( strcmp(boxcons->b, "massFlow") == 0 )
      {
         assert( strcmp(boxcons->b_unit, "kg_per_s") == 0);
         vars[1] = flowvariable;
      }
      else if ( strcmp(boxcons->b, "pressureIn") == 0 )
      {
         assert( strcmp(boxcons->b_unit, "bar") == 0);
         vars[1] = sourcevariable;
      }
      else if ( strcmp(boxcons->b, "pressureOut") == 0 )
      {
         assert( strcmp(boxcons->b_unit, "bar") == 0);
         vars[1] = targetvariable;
      }
      if ( strcmp(boxcons->c, "massFlow") == 0 )
      {
         assert( strcmp(boxcons->c_unit, "kg_per_s") == 0);
         vars[2] = flowvariable;
      }
      else if ( strcmp(boxcons->c, "pressureIn") == 0 )
      {
         assert( strcmp(boxcons->c_unit, "bar") == 0);
         vars[2] = sourcevariable;
      }
      else if ( strcmp(boxcons->c, "pressureOut") == 0 )
      {
         assert( strcmp(boxcons->c_unit, "bar") == 0);
         vars[2] = targetvariable;
      }

      SCIPdebugMessage("facets:<%i>\n", boxcons->numfacets);
      for (i = 0; i < boxcons->numfacets; i++)
      {
         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_facet#%s_%s#%i", arc->id, boxcons->id, i);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &facet_cons, consname, 3, vars, boxcons->facetcoeff[i], -SCIPinfinity(scip), boxcons->rhs[i]) );
         SCIP_CALL( SCIPaddCons(scip, facet_cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &facet_cons) );
      }
   }

   return SCIP_OKAY;
}

/** generate constraints for a nonlinear resistor with fixed flow direction */
static
SCIP_RETCODE generateNonLinResistorForStation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GAS_Arc*              arc,                /**< either compressor arc or controlvalve arc */
   SCIP_Bool             inRes               /**< boolean wheter or not in resistor should be created  */
   )
{
   GAS_CS*          cs;
   /* GAS_Controlvalve cv; */
   SCIP_VAR*        sourcepressure;
   SCIP_VAR*        targetpressure;
   SCIP_VAR*        pressureDiff;
   SCIP_VAR*        flow;
   SCIP_VAR*        vars3[3];
   SCIP_VAR*        vars4[4];
   SCIP_Real        ca = SCIP_INVALID;
   SCIP_Real        vals3[3];
   SCIP_Real        vals4[4];
   SCIP_CONS*       constraint;
   char             consname[SCIP_MAXSTRLEN];
   char             prefix[SCIP_MAXSTRLEN];

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( arc != NULL );

   flow = arc->flowvar;

   if ( arc->type != CS )
   {
      SCIPerrorMessage("Arc type %d not implemented for the function generateNonLinResistorForStation().\n", arc->type);
      return SCIP_ERROR;
   }

   cs = arc->detailed_info;
   if ( inRes )
   {
      (void) SCIPsnprintf(prefix, SCIP_MAXSTRLEN, "CS_ResInCons");
      sourcepressure = arc->sourcenode->pressurevar;
      targetpressure = cs->NLRin_pressure;
      pressureDiff   = cs->NLRin_pressureDiff;

      ca = ConstantNonlinearResistor(probdata->network, 0.5 * (arc->sourcenode->pressureMin + arc->sourcenode->pressureMax),
         cs->dragFactorIn, cs->diameterIn, probdata->papay);
   }
   else
   {
      (void) SCIPsnprintf(prefix, SCIP_MAXSTRLEN, "CS_ResOutCons");
      sourcepressure = cs->NLRout_pressure;
      targetpressure = arc->targetnode->pressurevar;
      pressureDiff   = cs->NLRout_pressureDiff;

      ca = ConstantNonlinearResistor(probdata->network, 0.5 * (arc->targetnode->pressureMin + arc->targetnode->pressureMax),
         cs->dragFactorOut, cs->diameterOut, probdata->papay);
   }
   assert( ca != SCIP_INVALID );

   /*
    *  create the coupling of the pressure variables with the pressure difference
    */

   /* Create the constraint: Delta = p_u - p_v*/
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s1#%s#", prefix, arc->id);

   /* variable array  */
   vars3[0] = pressureDiff;
   vars3[1] = sourcepressure;
   vars3[2] = targetpressure;

   /* coefficients array */
   vals3[0] =  1;
   vals3[1] = -1;
   vals3[2] =  1;

   /* Create linear constraint and add it to scip */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 3, vars3, vals3, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, constraint) );
   SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

   /*
    *  create the coupling between the pressureDifference and the binary variables:
    *  pressureDiff <= maxDiff * (1 - \sum config_binvars)
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s2#%s#", prefix, arc->id);

   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &constraint, consname, pressureDiff, cs->compressor_binvar,
         - SCIPcomputeVarUbLocal(scip, pressureDiff), - SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddCons(scip, constraint) );
   SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

    /*
     *  create the nonlinear resistor constraint
     *  p_u^2 - p_v^2 + pDiff^2 - 2 * ca * q^2 = 0.0
     */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "%s3#%s#", prefix, arc->id);

   vars4[0] = sourcepressure;
   vars4[1] = targetpressure;
   vars4[2] = pressureDiff;
   vars4[3] = flow;

   vals4[0] = 1.0;
   vals4[1] = - 1.0;
   vals4[2] = 1.0;
   vals4[3] = - 2.0 * ca;

#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
   SCIP_CALL( SCIPcreateConsBasicQuadraticNonlinear(scip, &constraint, consname, 0, NULL, NULL, 4, vars4, vars4, vals4, 0.0, 0.0) );
#else
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &constraint, consname, 0, NULL, NULL, 4, vars4, vars4, vals4, 0.0, 0.0) );
#endif
   SCIP_CALL( SCIPaddCons(scip, constraint) );
   SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

   return SCIP_OKAY;
}

/** generates constraints for the box constrained compressor model, see documentation for details */
static
SCIP_RETCODE generateBoxConstraintModel(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GAS_Arc*              arc                 /**< the compressor arc */
   )
{
   GAS_CS*        cs;
   GAS_CSBoxCons* boxcons;
   GAS_Node*      Sourcenode;
   GAS_Node*      Targetnode;
   SCIP_VAR*      sourcepressure;
   SCIP_VAR*      targetpressure;
   SCIP_VAR*      flow;
   SCIP_VAR*      vars2[2];
   SCIP_VAR*      vars3[3];
   SCIP_VAR*      longvars[(((GAS_CS*)arc->detailed_info)->numconfigurations + 1)];
   SCIP_Real      vals2[2];
   SCIP_Real      vals3[3];
   SCIP_Real      longvals[(((GAS_CS*)arc->detailed_info)->numconfigurations + 1)];
   SCIP_Real      lhs;
   SCIP_Real      rhs;
   SCIP_CONS*     constraint;
   SCIP_CONS*     compressor_binvar_cons;
   char           consname[SCIP_MAXSTRLEN];
   int            k;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->network != NULL );

   cs             = arc->detailed_info;
   Sourcenode     = arc->sourcenode;
   Targetnode     = arc->targetnode;
   sourcepressure = Sourcenode->pressurevar;
   targetpressure = Targetnode->pressurevar;
   flow           = arc->flowvar;

   /* link compressor binvar to posflowbinvar */
   (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_link#%s", arc->id);
   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &compressor_binvar_cons, consname, arc->positiveFlowBinvar, cs->compressor_binvar, -1.0, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, compressor_binvar_cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &compressor_binvar_cons) );

   /*
    *  coupling of compressor_binvar with config_binvars
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_BinVarSum#%s#", arc->id);

   /* fill longvars and longvals */
   longvars[0] = cs->compressor_binvar;
   longvals[0] = 1.0;
   for (k = 1; k <= cs->numconfigurations; ++k)
   {
      longvars[k] = cs->boxcons[k-1].config_binvar;
      longvals[k] = - 1.0;
   }

   lhs = 0.0;
   rhs = 0.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, cs->numconfigurations + 1, longvars, longvals, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, constraint) );
   SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

   /*
    *  flow has to be zero if compressor is turned off
    */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_flowCons#%s#", arc->id);

   /* fill longvars and longvals */
   vars2[0] = flow;
   vals2[0] = 1.0;
   vars2[1] = cs->compressor_binvar;
   vals2[1] = - arc->flowMax;

   lhs = - SCIPinfinity(scip);
   rhs = 0.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 2, vars2, vals2, lhs, rhs) );
   SCIP_CALL( SCIPaddCons(scip, constraint) );
   SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

   /*
    *  constraints only needed if there is a bypass
    */
   if ( cs->internalBypassRequired == 1 )
   {
      /* compressor can only be operating in one configuration or be in bypass mode or off */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_or_Bypass#%s#", arc->id);

      vars2[0] = probdata->VALVE_binvars[cs->bypassPosition];
      vals2[0] = 1.0;
      vars2[1] = cs->compressor_binvar;
      vals2[1] = 1.0;

      lhs = 0.0;
      rhs = 1.0;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 2, vars2, vals2, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      /* If flow on the compressorStation arc should be positive. This can either be on the bypass or the compressor itself. */
      if ( arc->flowMin > 0.0 )
      {
         SCIP_VAR*  flowValve;
         int valveposition = arc->arcposition + 1;

         assert( probdata->network->arcs_ptr[valveposition].type == VALVE );

         flowValve = probdata->network->arcs_ptr[valveposition].flowvar;

         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_MinFlowCons#%s#", arc->id);
         vars2[0] = arc->flowvar;
         vars2[1] = flowValve;
         vals2[0] = 1.0;
         vals2[1] = 1.0;

         lhs = arc->flowMin;
         rhs = arc->flowMax;

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 2, vars2, vals2, lhs, rhs) );
         SCIP_CALL( SCIPaddCons(scip, constraint) );
         SCIP_CALL( SCIPreleaseCons(scip, &constraint) );
      }
   }

   /*
    *  generate constraints for nonlinear in or out resistors
    */
   if ( ! SCIPisZero(scip, cs->dragFactorIn) )
   {
      SCIP_CALL( generateNonLinResistorForStation(scip, probdata, arc, TRUE) );
      sourcepressure = cs->NLRin_pressure;
   }

   if ( ! SCIPisZero(scip, cs->dragFactorOut) )
   {
      SCIP_CALL( generateNonLinResistorForStation(scip, probdata, arc, FALSE) );
      targetpressure = cs->NLRout_pressure;
   }

   /*
    *  generate the box constraint model for each configuration
    */
   for (k = 0; k < cs->numconfigurations; ++k)
   {
      boxcons = &(cs->boxcons[k]);

      /* absolute pressure increase has to be in the given bounds:
       * boxcons->pressureIncAbsMin <=  BCM_pout - BCM_pin <= boxcons->pressureIncAbsMax */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "BCM_AbsPressIncCons#%s_%s#", arc->id, boxcons->id);

      vars2[0] = boxcons->BCM_pout;
      vars2[1] = boxcons->BCM_pin;
      vals2[0] = 1.0;
      vals2[1] = - 1.0;

      lhs = boxcons->pressureIncAbsMin;
      rhs = boxcons->pressureIncAbsMax;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 2, vars2, vals2, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      /* minimal relative pressure increase
       * boxcons->pressureIncRelMin * BCM_pin <= BCM_pout */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "BCM_minRelPressInc#%s_%s#", arc->id, boxcons->id);
      vals2[0] = - 1.0;
      vals2[1] = boxcons->pressureIncRelMin;

      lhs = - SCIPinfinity(scip);
      rhs = 0.0;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 2, vars2, vals2, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      /* maximal relative pressure increase
       * BCM_pout - boxcons->pressureIncRelMax * BCM_pin <= 0.0 */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "BCM_maxRelPressInc#%s_%s#", arc->id, boxcons->id);

      vals2[0] = 1.0;
      vals2[1] = - boxcons->pressureIncRelMax;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 2, vars2, vals2, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      /* if the compressor is running with this configuration then the flow
       * over the compressor arc has to be equal to the flow of this configuration:
       * PART 1: q_a - q_c <= (1 - s_c ) (q_a^max - q_c^min) */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "BCM_flowCons1#%s_%s#", arc->id, boxcons->id);

      vars3[0] = flow;
      vars3[1] = boxcons->BCM_flow;
      vars3[2] = boxcons->config_binvar;
      vals3[0] = 1;
      vals3[1] = -1;
      vals3[2] = arc->flowMax - boxcons->massFlowMin;

      lhs = - SCIPinfinity(scip);
      rhs = vals3[2];

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 3, vars3, vals3, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      /* PART 2: (1 - s_c ) (q_a^min - q_c^max) <= q_a - q_c  */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "BCM_flowCons2#%s_%s#", arc->id, boxcons->id);

      /* note that the lower bound on q_a depends on the presence or absence of a bypass and the sign of flowMin */
      vals3[2] = SCIPcomputeVarLbLocal(scip, flow) - boxcons->massFlowMax;

      lhs = vals3[2];
      rhs = SCIPinfinity(scip);

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 3, vars3, vals3, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      /* if the compressor is running with this configuration then the pressure
       * at the end of the in-resistor and the pressure at the compressor entrance
       * have to be equal.
       * the pressure at the in-resistor end is either the pressure at the sourcenode - pressureLossIn
       * or the pressure CS_pressureNLRin, sourcepressure is set to the corresponding pressure.
       * PART 1: sourcepressure - BCM_pin <= b_cf * pressureLossIn + (1 - b_cf) * (sourcepressure_max - pin_min) */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "BCM_linkPinCons1#%s_%s#", arc->id, boxcons->id);
      vars3[0] = sourcepressure;
      vars3[1] = boxcons->BCM_pin;
      vars3[2] = boxcons->config_binvar;
      vals3[0] = 1.0;
      vals3[1] = - 1.0;
      vals3[2] = Sourcenode->pressureMax - (boxcons->pressureInMin + cs->pressureLossIn);

      lhs = - SCIPinfinity(scip);
      rhs = Sourcenode->pressureMax - boxcons->pressureInMin;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 3, vars3, vals3, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      /* PART 2:  b_cf * pressureLossIn + (1 - b_cf) * (sourcepressure_min - pin_max) <= sourcepressure - BCM_pin*/
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "BCM_linkPinCons2#%s_%s#", arc->id, boxcons->id);

      vals3[2] = Sourcenode->pressureMin - (boxcons->pressureInMax + cs->pressureLossIn);

      lhs = Sourcenode->pressureMin - boxcons->pressureInMax;
      rhs = SCIPinfinity(scip);

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 3, vars3, vals3, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      /* Analog for the out-pressures
       * PART 1: BCM_pout - targetpressure <= b_cf * pLossOut + (1-b_cf) * (pout_max - targetpressure_min) */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "BCM_linkPoutCons1#%s_%s#", arc->id, boxcons->id);
      vars3[0] = boxcons->BCM_pout;
      vars3[1] = targetpressure;
      vars3[2] = boxcons->config_binvar;
      vals3[0] = 1.0;
      vals3[1] = - 1.0;
      vals3[2] = boxcons->pressureOutMax - (Targetnode->pressureMin  + cs->pressureLossOut);

      lhs = - SCIPinfinity(scip);
      rhs = boxcons->pressureOutMax - Targetnode->pressureMin;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 3, vars3, vals3, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );

      /* PART 2: b_cf * pLossOut + (1-b_cf) * (pout_min - targetpressure_max) <= BCM_pout - targetpressure */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "BCM_linkPoutCons2#%s_%s#", arc->id, boxcons->id);

      vals3[2] = boxcons->pressureOutMin -( Targetnode->pressureMax  + cs->pressureLossOut);

      lhs = boxcons->pressureOutMin - Targetnode->pressureMax;
      rhs = SCIPinfinity(scip);

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &constraint, consname, 3, vars3, vals3, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(scip, constraint) );
      SCIP_CALL( SCIPreleaseCons(scip, &constraint) );
   }

   return SCIP_OKAY;
}

/** generates the constraints for an ideal compressor
 *
 *  See model (6.2) in the dissertation of Oliver Habeck.
 */
static
SCIP_RETCODE generateIdealCSConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GAS_Arc*              arc                 /**< compressor arc */
   )
{
   char          consname[SCIP_MAXSTRLEN];
   SCIP_CONS*    compressor_pressure_cons1;
   SCIP_CONS*    compressor_pressure_cons2;
   SCIP_CONS*    compressor_binvar_cons;
   GAS_CS*       compressor;
   SCIP_VAR*     Sourcevariable;
   SCIP_VAR*     Targetvariable;
   GAS_Node*     Sourcenode;
   GAS_Node*     Targetnode;
   SCIP_VAR*     Flowvariable;
   SCIP_VAR*     vars2[2];
   SCIP_Real     vals2[2];
   SCIP_VAR*     vars3[3];
   SCIP_Real     vals3[3];
   SCIP_CONS*    flow_cons;
   SCIP_Real     deltaPressureMax;
   SCIP_Real     deltaPressureMin;
   SCIP_CONS*    pressure_diff_cons1;
   SCIP_CONS*    pressure_diff_cons2;
   SCIP_Real     pIncMin;
   SCIP_Real     pIncMax;

   compressor     = arc->detailed_info;
   Sourcenode     = arc->sourcenode;
   Targetnode     = arc->targetnode;
   Sourcevariable = Sourcenode->pressurevar;
   Targetvariable = Targetnode->pressurevar;
   Flowvariable   = arc->flowvar;

   /*
    *  possible pressure increase!
    */
   pIncMin = probdata->idealCSminIncrease;
   pIncMax = probdata->idealCSmaxIncrease;

   if ( pIncMin > pIncMax )
   {
      SCIPerrorMessage("Wrong parameters for pressure increase in an ideal compressor station set:\n");
      SCIPerrorMessage("The minimal pressure increase has to be less or equal than the maximal increase.\n");
      return SCIP_PARAMETERWRONGVAL;
   }

   /* link compressor binvar to posflowbinvar */
   (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_link#%s", arc->id);
   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &compressor_binvar_cons, consname, arc->positiveFlowBinvar, compressor->compressor_binvar, -1.0, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, compressor_binvar_cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &compressor_binvar_cons) );

   /*
    *  flow is zero if compressor is turned off
    */
   /* flow on q_a has to be zero if compressor is off:  q_a - b_a * flowMax \leq 0.0 */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_flowCons#%s#", arc->id);
   SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &flow_cons, consname, Flowvariable, compressor->compressor_binvar, - arc->flowMax, - SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddCons(scip, flow_cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &flow_cons) );

   /*
    *  either compressor is on or bypass is open.
    *  If the minFlow on the compressor arc is positive, then either flow on cs or valve has to be positive
    */
   /* check if internalBypassRequired is 1 or 0 */
   if ( compressor->internalBypassRequired == 1 )
   {
      SCIP_CONS* BinVarSumCons;

      /* compressor is either in bypass, on or off */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_BinVarSum#%s#", arc->id);

      vars2[0] = compressor->compressor_binvar;
      vars2[1] = probdata->VALVE_binvars[compressor->bypassPosition];
      vals2[0] = 1.0;
      vals2[1] = 1.0;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &BinVarSumCons, consname, 2, vars2, vals2, 0.0, 1.0) );
      SCIP_CALL( SCIPaddCons(scip, BinVarSumCons) );
      SCIP_CALL( SCIPreleaseCons(scip, &BinVarSumCons) );

      /* if should be positive flow on the compressorStation
       * this can either be on the bypass or the compressor itself
       */
      if ( arc->flowMin > 0.0 )
      {
         SCIP_CONS* MinFlowCons;
         SCIP_VAR*  flowValve;
         int valveposition = arc->arcposition + 1;

         assert( probdata->network->arcs_ptr[valveposition].type == VALVE );

         flowValve = probdata->network->arcs_ptr[valveposition].flowvar;

         (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_MinFlowCons#%s#", arc->id);
         vars2[0] = arc->flowvar;
         vars2[1] = flowValve;

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &MinFlowCons, consname, 2, vars2, vals2, arc->flowMin, SCIPinfinity(scip)) );
         SCIP_CALL( SCIPaddCons(scip, MinFlowCons) );
         SCIP_CALL( SCIPreleaseCons(scip, &MinFlowCons) );
      }
   }

   /*
    *  pressure at the entrance has to be big enough
    */
   if ( compressor->pressureInMin > Sourcenode->pressureMin )
   {
      /* create linear constraint compressor_pressure_cons1: c^in_min * b_a - p^in <= 0  <=> p^in >= c^in_min * b_a */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_minPinCons#%s#", arc->id);

      vars2[0] = compressor->compressor_binvar;
      vars2[1] = Sourcevariable;
      vals2[0] = compressor->pressureInMin;
      vals2[1] = -1;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &compressor_pressure_cons1, consname, 2, vars2, vals2, -SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCons(scip, compressor_pressure_cons1) );
      SCIP_CALL( SCIPreleaseCons(scip, &compressor_pressure_cons1) );
   }

   /*
    *  cannot produce pressure above pressureOutMax
    */
   if ( compressor->pressureOutMax < Targetnode->pressureMax )
   {
      /* create linear constraint compressor_pressure_cons2: p^out + (p^out_max - c^out_max) * b_a <= p^out_max */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "CS_maxPoutCons#%s#", arc->id);

      vars2[0] = Targetvariable;
      vars2[1] = compressor->compressor_binvar;
      vals2[0] = 1;
      vals2[1] = Targetnode->pressureMax - compressor->pressureOutMax;

      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &compressor_pressure_cons2, consname, 2, vars2, vals2, -SCIPinfinity(scip), Targetnode->pressureMax) );
      SCIP_CALL( SCIPaddCons(scip, compressor_pressure_cons2) );
      SCIP_CALL( SCIPreleaseCons(scip, &compressor_pressure_cons2) );
   }

   /*
    *  pressure can be increased between pIncMin and pIncMax;
    */
   if ( probdata->csusemaxincrease )
   {
      /* constraint concerning pressure */
      deltaPressureMax = (Targetnode->pressureMax) - (Sourcenode->pressureMin);
      deltaPressureMin = (Targetnode->pressureMin) - (Sourcenode->pressureMax);

      /* constraint name */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "idealCS_Cons1#%s#", arc->id);

      /* array, that points to variables */
      vars3[0] = Targetvariable;
      vars3[1] = Sourcevariable;
      vars3[2] = compressor->compressor_binvar;

      /* array, with coefficients */
      vals3[0] =  1;
      vals3[1] =  - 1;
      vals3[2] = deltaPressureMax - pIncMax;

      /* create linear constraint: p^out - p^in + (\Delta_max - p^inc_max) * b_a <= \Delta_max */
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &pressure_diff_cons1, consname, 3, vars3, vals3, -SCIPinfinity(scip), deltaPressureMax) );
      SCIP_CALL( SCIPaddCons(scip, pressure_diff_cons1) );
      SCIP_CALL( SCIPreleaseCons(scip, &pressure_diff_cons1) );

      /* constraint name */
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "idealCS_Cons2#%s#", arc->id);

      /* array, with coefficients */
      vals3[0] =  1;
      vals3[1] =  - 1;
      vals3[2] = deltaPressureMin - pIncMin;

      /* create linear constraint: p^out - p^in + (\Delta_min - p^inc_min) * b_a >= \Delta_min */
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &pressure_diff_cons2, consname, 3, vars3, vals3, deltaPressureMin, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, pressure_diff_cons2) );
      SCIP_CALL( SCIPreleaseCons(scip, &pressure_diff_cons2) );
   }

   return SCIP_OKAY;
}

/** calls the correct function for generating the compressor station constraints depending on the model which is used */
static
SCIP_RETCODE generateCompressorStationConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GAS_Arc*              arc                 /**< compressor arc */
   )
{
   assert( scip != NULL );
   assert( probdata != NULL );
   assert( arc != NULL );

   if ( probdata->boxConstraintModel )
   {
      /* use box constraint model*/
      SCIP_CALL( generateBoxConstraintModel(scip, probdata, arc) );
      if ( probdata->additionalFacets )
      {
         SCIP_CALL( generateAdditionalFacets(scip, probdata, arc) );
      }
   }
   else
   {
      /* use idealized compressor model */
      SCIP_CALL( generateIdealCSConstraints(scip, probdata, arc) );
   }

   return SCIP_OKAY;
}

/** The following function generates the constraints for a nonlinear resistor,
 *  including the coupling of the pressure variables with the pressureSquare
 *  variable and the flow and pressure difference variables with their abspower
 *  variabiales
 */
static
SCIP_RETCODE generateNonlinearResistorConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   arcposition         /**< poairion of arc */
   )
{
   SCIP_VAR*     Sourcepressure;
   SCIP_VAR*     Targetpressure;
   SCIP_VAR*     Flowvariable;
   SCIP_VAR*     SourcePIvar;
   SCIP_VAR*     TargetPIvar;
   GAS_Node*     Sourcenode;
   GAS_Node*     Targetnode;
   SCIP_Real     SourceMaxMinusTargetMin;
   SCIP_Real     SourceMinMinusTargetMax;
   SCIP_VAR*     vars3[3];
   SCIP_VAR*     vars4[4];
   SCIP_Real     vals3[3];
   SCIP_Real     vals4[4];
   SCIP_Real     meanpressure;
   SCIP_CONS*    AbsFlow_cons;
   SCIP_CONS*    AbsDelta_cons;
   SCIP_CONS*    pressureDiffCons;
   SCIP_CONS*    NLResistorCons;
   char          consname[SCIP_MAXSTRLEN];
   GAS_Arc*      gas_arc;
   GAS_Resistor* resistor;
   int           resistorposition;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   gas_arc                 = &probdata->network->arcs_ptr[arcposition];
   resistor                = gas_arc->detailed_info;
   resistorposition        = resistor->NLR_position;
   Sourcenode              = probdata->network->arcs_ptr[arcposition].sourcenode;
   Targetnode              = probdata->network->arcs_ptr[arcposition].targetnode;
   Sourcepressure          = Sourcenode->pressurevar;
   Targetpressure          = Targetnode->pressurevar;
   Flowvariable            = gas_arc->flowvar;
   SourceMaxMinusTargetMin = (Sourcenode -> pressureMax) - (Targetnode -> pressureMin);
   SourceMinMinusTargetMax = (Sourcenode -> pressureMin)  - (Targetnode -> pressureMax);
   meanpressure            = ComputeMeanPressure(Sourcenode->pressureMin, Sourcenode->pressureMax, Targetnode->pressureMin, Targetnode->pressureMax);
   SourcePIvar             = Sourcenode->PIvar;
   TargetPIvar             = Targetnode->PIvar;

   /*
    *  create the coupling of the flowvariable with its variable for |q|*q
    */

   /* Create the constraint: |flowvar|*flowvar = AbsFlowVar  */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "NLR_AbsFlowCons#%s#", gas_arc->id);

#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
   SCIP_CALL( SCIPcreateConsBasicSignpowerNonlinear(scip, &AbsFlow_cons, consname, Flowvariable, probdata->NLR_AbsFlowVars[resistorposition], 2.0, 0.0, -1.0, 0.0, 0.0) );
#else
   SCIP_CALL( SCIPcreateConsBasicAbspower(scip, &AbsFlow_cons, consname, Flowvariable, probdata->NLR_AbsFlowVars[resistorposition], 2.0, 0.0, -1.0, 0.0, 0.0) );
#endif
   SCIP_CALL( SCIPaddCons(scip, AbsFlow_cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &AbsFlow_cons) );

   /*
    *  create the coupling of the pressure variables with the Delta variable for the difference
    */

   /* Create the constraint: Delta = p_u - p_v*/
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "NLR_PressureDiffCons#%s#", gas_arc->id);

   /* variable array  */
   vars3[0] = probdata->NLR_DeltaVars[resistorposition];
   vars3[1] = Sourcepressure;
   vars3[2] = Targetpressure;

   /* coefficients array */
   vals3[0] =  1;
   vals3[1] = -1;
   vals3[2] =  1;

   /* Create linear constraint and add it to scip */
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &pressureDiffCons, consname, 3, vars3, vals3, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, pressureDiffCons) );
   SCIP_CALL( SCIPreleaseCons(scip, &pressureDiffCons) );

   /*
    *  create the coupling of the Delta variable with the abspower variable
    */

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "NLR_AbsDeltaCons#%s#", gas_arc->id);

#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
   SCIP_CALL( SCIPcreateConsBasicSignpowerNonlinear(scip, &AbsDelta_cons, consname, probdata->NLR_DeltaVars[resistorposition],
         probdata->NLR_AbsDeltaVars[resistorposition], 2.0, 0.0, -1.0, 0.0, 0.0) );
#else
   SCIP_CALL( SCIPcreateConsBasicAbspower(scip, &AbsDelta_cons, consname, probdata->NLR_DeltaVars[resistorposition],
         probdata->NLR_AbsDeltaVars[resistorposition], 2.0, 0.0, -1.0, 0.0, 0.0) );
#endif
   SCIP_CALL( SCIPaddCons(scip, AbsDelta_cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &AbsDelta_cons) );

   /*
    *  create the coupling of the Delta variable with the flow direction variable
    */
   if ( ! probdata->noFlowBinvars )
   {
      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "NLR_positivePressureDiffCons#%s#", gas_arc->id);

      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &pressureDiffCons, consname, probdata->NLR_DeltaVars[resistorposition],
            gas_arc->positiveFlowBinvar, - SourceMaxMinusTargetMin, - SCIPinfinity(scip), 0.0) );
      SCIP_CALL( SCIPaddCons(scip, pressureDiffCons) );
      SCIP_CALL( SCIPreleaseCons(scip, &pressureDiffCons) );

      (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "NLR_negativePressureDiffCons#%s#", gas_arc->id);

      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &pressureDiffCons, consname, probdata->NLR_DeltaVars[resistorposition],
            gas_arc->negativeFlowBinvar, - SourceMinMinusTargetMax, 0.0, SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, pressureDiffCons) );
      SCIP_CALL( SCIPreleaseCons(scip, &pressureDiffCons) );
   }

   /*
    *  create the nonlinear resistor constraint with the squared pressure variables and abspower variables
    */

   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "NLR_resistorCons#%s#", gas_arc->id);

   /* variable array */
   vars4[0] = SourcePIvar;
   vars4[1] = TargetPIvar;
   vars4[2] = probdata->NLR_AbsDeltaVars[resistorposition];
   vars4[3] = probdata->NLR_AbsFlowVars[resistorposition];

   /* coefficients array */
   vals4[0] =  1;
   vals4[1] = -1;
   vals4[2] =  1;
   vals4[3] =  -2.0 * ConstantNonlinearResistor(probdata->network, meanpressure, resistor->dragFactor, resistor->diameter, probdata->papay);

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &NLResistorCons, consname, 4, vars4, vals4, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, NLResistorCons) );
   SCIP_CALL( SCIPreleaseCons(scip, &NLResistorCons) );

   return SCIP_OKAY;
}

/** generate linear resitor constraints */
static
SCIP_RETCODE generateLinearResistorConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   arcposition,        /**< position of arc */
   int                   resistorposition    /**< position of resistor */
   )
{
   SCIP_VAR* Sourcevariable;
   SCIP_VAR* Targetvariable;
   SCIP_VAR* Flowvariable;
   SCIP_Real epsilon;
   SCIP_Real vals3[3];
   SCIP_VAR* vars3[3];
   SCIP_Real vals4[4];
   SCIP_VAR* vars4[4];
   SCIP_Real vals5[5];
   SCIP_VAR* vars5[5];
   SCIP_CONS* LinResCons;
   char consname[SCIP_MAXSTRLEN];
   GAS_Arc* gas_arc;
   GAS_Resistor* resistor;

   assert( probdata != NULL );
   assert( probdata->network != NULL );

   gas_arc = &probdata->network->arcs_ptr[arcposition];
   resistor = gas_arc->detailed_info;
   Sourcevariable = probdata->network->arcs_ptr[arcposition].sourcenode->pressurevar;
   Targetvariable = probdata->network->arcs_ptr[arcposition].targetnode->pressurevar;
   Flowvariable = gas_arc->flowvar;
   epsilon = probdata->network->normDensity / 3600.0;

   /* constraint LinResCons1: p_v - p_u = - pL * LR_negFlowDir + pL * LR_posFlowDir + pL * LR_smoothingFlow */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LR_pressureLossEQ#%s#", gas_arc->id);

   vals5[0] = 1;
   vals5[1] = -1;
   vals5[2] = -resistor->pressureLoss;
   vals5[3] = resistor->pressureLoss;
   vals5[4] = resistor->pressureLoss;

   vars5[0] = Targetvariable;
   vars5[1] = Sourcevariable;
   vars5[2] = probdata->LR_negFlowDir[resistorposition];
   vars5[3] = probdata->LR_posFlowDir[resistorposition];
   vars5[4] = probdata->LR_smoothingFlow[resistorposition];

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &LinResCons, consname, 5, vars5, vals5, 0.0, 0.0) );
   SCIP_CALL( SCIPaddCons(scip, LinResCons) );
   SCIP_CALL( SCIPreleaseCons(scip, &LinResCons) );

   /* constraint LinResCons2: flowMin * LR_negFlowDir + epsilon * LR_posFlowDir + epsilon * LR_smoothingFlow - flowvar <= 0 */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LR_flowCons1#%s#", gas_arc->id);

   vals4[0] = gas_arc->flowMin;
   vals4[1] = epsilon;
   vals4[2] = epsilon;
   vals4[3] = - 1.0;

   vars4[0] = probdata->LR_negFlowDir[resistorposition];
   vars4[1] = probdata->LR_posFlowDir[resistorposition];
   vars4[2] = probdata->LR_smoothingFlow[resistorposition];
   vars4[3] = Flowvariable;

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &LinResCons, consname, 4, vars4, vals4, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddCons(scip, LinResCons) );
   SCIP_CALL( SCIPreleaseCons(scip, &LinResCons) );

   /*  constraint LinResCons3: 0.0 <= - epsilon * LR_negFlowDir + flowMax * LR_posFlowDir + epsilon * LR_smoothingFlow - flowvar */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LR_flowCons2#%s#", gas_arc->id);

   vals4[0] = - epsilon;
   vals4[1] = gas_arc->flowMax;
   vals4[2] = epsilon;
   vals4[3] = - 1;

   /* vars3 stays the same as for LinResCons3*/

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &LinResCons, consname, 4, vars4, vals4, 0.0, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, LinResCons) );
   SCIP_CALL( SCIPreleaseCons(scip, &LinResCons) );

   /* constraint LinResCons4: LR_negFlowDir + LR_posFlowDir - LR_smoothingFlow  <= 1 */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LR_flowDirCons1#%s#", gas_arc->id);

   vals3[0] =  1.0;
   vals3[1] =  1.0;
   vals3[2] = -1.0;

   vars3[0] = probdata->LR_negFlowDir[resistorposition];
   vars3[1] = probdata->LR_posFlowDir[resistorposition];
   vars3[2] = probdata->LR_smoothingFlow[resistorposition];

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &LinResCons, consname, 3, vars3, vals3, -SCIPinfinity(scip), 1.0) );
   SCIP_CALL( SCIPaddCons(scip, LinResCons) );
   SCIP_CALL( SCIPreleaseCons(scip, &LinResCons) );

   /*  constraint LinResCons5:  LR_negFlowDir + LR_posFlowDir +LR_smoothingFlow <= 1.0 */
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LR_flowDirCons2#%s#", gas_arc->id);

   vals3[2] = 1.0;

   /* vars3 stays the same as for LinResCons5*/

   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &LinResCons, consname, 3, vars3, vals3, -SCIPinfinity(scip), 1.0) );
   SCIP_CALL( SCIPaddCons(scip, LinResCons) );
   SCIP_CALL( SCIPreleaseCons(scip, &LinResCons) );

   return SCIP_OKAY;
}

/** generate objective constraints */
static
SCIP_RETCODE generateObjectiveConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;                  /* problem data */
   SCIP_CONS* cons;
   char consname[SCIP_MAXSTRLEN];
   int i;
   GAS_Node* N;
   SCIP_VAR* vars2[2];
   SCIP_Real vals2[2];

   assert( scip != NULL );

   probdata = SCIPgetProbData(scip);

   assert( probdata != NULL );
   assert( probdata->objective != NULL );

   /*
    *  generate constraints for minimizing the (maximal) violation of pressure bounds
    */
   if ( probdata->relaxLowerBounds || probdata->relaxUpperBounds )
   {
      vars2[1] = probdata->objective;
      vals2[0] = 1.0;

      for (i = 0; i < probdata->network->numnodes; ++i)
      {
         SCIP_Real bound;

         N = &(probdata->network->nodes_ptr[i]);

         if ( probdata->minSlackPerBound )
            vars2[1] = N->pSlackVar;

         vars2[0] = N->pressurevar;

         if ( probdata->relaxLowerBounds )
         {
            bound = MAX(N->netPressureMin, N->scnPressureMin);
            if ( SCIPisLT(scip, N->pressureMin, bound) )
            {
               vals2[1] = 1.0;

               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "relaxLBcons#%s#", N->id);
               SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 2, vars2, vals2, bound, SCIPinfinity(scip)) );
               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }
         }

         if ( probdata->relaxUpperBounds )
         {
            bound = MIN3(N->netPressureMax, N->scnPressureMax, N->arcPressureMax);
            if ( SCIPisGT(scip, N->pressureMax, bound) )
            {
               vals2[1] = - 1.0;

               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "relaxUBcons#%s#", N->id);
               SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 2, vars2, vals2, -SCIPinfinity(scip), bound) );
               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }
         }
      }
   }

   if ( probdata->relaxCVbounds && probdata->network->numcontrolvalves > 0 )
   {
      GAS_Arc* arc;
      GAS_Controlvalve* cv;
      SCIP_VAR* vars[3];
      SCIP_Real vals[3];

      for (i = 0; i < probdata->network->numarcs; ++i)
      {
         arc = &(probdata->network->arcs_ptr[i]);

         if ( arc->type != CONTROLVALVE )
            continue;

         cv = arc->detailed_info;

         vars[0] = arc->targetnode->pressurevar;
         vars[1] = cv->binvar;
         if ( probdata->minSlack )
            vars[2] = probdata->objective;
         else if ( probdata->minSlackPerBound )
            vars[2] = arc->targetnode->pSlackVar;

         vals[0] = 1.0;
         vals[1] = arc->targetnode->pressureMax - (cv->pressureOutMax - cv->pressureLossOut);
         vals[2] = - 1.0;

         if ( SCIPisPositive(scip, vals[1]) )
         {
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "relaxPoutMaxCons#%s#", arc->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 3, vars, vals, -SCIPinfinity(scip), arc->targetnode->pressureMax) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }

         vars[0] = arc->sourcenode->pressurevar;
         if ( probdata->minSlackPerBound )
            vars[2] = arc->sourcenode->pSlackVar;

         vals[1] = arc->sourcenode->pressureMin - (cv->pressureInMin + cv->pressureLossIn);
         vals[2] = 1.0;

         if ( SCIPisNegative(scip, vals[1]) )
         {
            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "relaxPinMinCons#%s#", arc->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 3, vars, vals, arc->sourcenode->pressureMin, SCIPinfinity(scip)) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }
   }

   return SCIP_OKAY;
}

/** generate component cuts
 *
 *  Produce connected components, when cutting at active elements. Then the cut accross each component \f$S\f$ with non-zero
 *  flow has to provide enough capacity:
 *  \f[
 *    \sum_{a \in \delta^+(S)} u_a - \sum_{a \in \delta^-(S)} \ell_a \geq \sum_{v \in S} b_v,
 *  \f]
 *  where \f$\ell\f$ and \f$u\f$ are lower and upper bounds on the flow, respectively. The supply/demand is given by \f$b\f$.
 *
 *  These cuts can be strengthened with the right hand side and often yield set-covering constraints (which are
 *  detected in presolving).
 */
static
SCIP_RETCODE generateComponentCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network             /**< network */
   )
{
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   GAS_Node** queue;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Bool* inqueue;
   int* component;
   int nvisited = 0;
   int ngen = 0;
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, network->numarcs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, network->numarcs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &queue, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inqueue, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &component, network->numnodes) );

   /* init data */
   for (i = 0; i < network->numnodes; ++i)
   {
      inqueue[i] = FALSE;
      component[i] = -1;
   }

   /* loop through nodes */
   for (i = 0; i < network->numnodes; ++i)
   {
      if ( component[i] < 0 )
      {
         SCIP_Real balance = 0.0;
         int lqueue = 0;
         int uqueue = 1;

         queue[0] = &network->nodes_ptr[i];
         assert( 0 <= queue[0]->nodeposition && queue[0]->nodeposition < network->numnodes );
         inqueue[queue[0]->nodeposition] = TRUE;
         assert( queue[0]->nodeposition == i );

         /* BFS loop through graph to compute component */
         while ( uqueue - lqueue > 0 )
         {
            GAS_Arc* arc;
            GAS_Node* v;
            GAS_Node* w;

            v = queue[lqueue++];
            assert( 0 <= v->nodeposition && v->nodeposition < network->numnodes );
            inqueue[v->nodeposition] = FALSE;
            component[v->nodeposition] = i;
            ++nvisited;

            /* update balance of component */
            balance += v->flow;

            /* check neighbors */
            arc = v->inarcs;
            while ( arc != NULL )
            {
               assert( 0 <= arc->arcposition && arc->arcposition < network->numarcs );
               w = arc->sourcenode;
               assert( 0 <= w->nodeposition && w->nodeposition < network->numnodes );
               if ( component[w->nodeposition] < 0 && ! inqueue[w->nodeposition] )
               {
                  /* only use non-active acrs */
                  if ( arc->type != VALVE && arc->type != CONTROLVALVE && arc->type != CS )
                  {
                     queue[uqueue++] = w;
                     inqueue[w->nodeposition] = TRUE;
                  }
               }
               arc = arc->next_inarc;
            }

            arc = v->outarcs;
            while ( arc != NULL )
            {
               assert( 0 <= arc->arcposition && arc->arcposition < network->numarcs );
               w = arc->targetnode;
               assert( 0 <= w->nodeposition && w->nodeposition < network->numnodes );
               if ( component[w->nodeposition] < 0 && ! inqueue[w->nodeposition] )
               {
                  /* only use non-active acrs */
                  if ( arc->type != VALVE && arc->type != CONTROLVALVE && arc->type != CS )
                  {
                     queue[uqueue++] = w;
                     inqueue[w->nodeposition] = TRUE;
                  }
               }
               arc = arc->next_outarc;
            }
         }

         /* if balance is nonzero, generate cut */
         if ( ! SCIPisZero(scip, balance) )
         {
            SCIP_Bool flipped = FALSE;
            SCIP_Real rhs;
            int cnt = 0;
            int k;

            if ( balance < 0.0 )
            {
               flipped = TRUE;
               rhs = - balance;
            }
            else
               rhs = balance;

            /* loop through nodes in component */
            for (k = 0; k < network->numnodes; ++k)
            {
               /* only consider nodes in current component */
               if ( component[k] == i )
               {
                  GAS_Node* node;
                  GAS_Valve* valve;
                  GAS_Controlvalve* cv;
                  GAS_CS* cs;
                  GAS_Arc* arc;

                  node = &network->nodes_ptr[k];

                  /* loop over outgoing arcs */
                  arc = node->outarcs;
                  while ( arc != NULL )
                  {
                     assert( SCIPisLE(scip, arc->flowMin, 0.0) );
                     assert( SCIPisGE(scip, arc->flowMax, 0.0) );

                     /* if arc is accross the cut */
                     if ( component[arc->targetnode->nodeposition] != i )
                     {
                        if ( arc->type == VALVE )
                        {
                           valve = (GAS_Valve*) arc->detailed_info;
                           vars[cnt] = valve->valve_binvar;
                        }
                        else if ( arc->type == CONTROLVALVE )
                        {
                           cv = (GAS_Controlvalve*) arc->detailed_info;
                           vars[cnt] = cv->binvar;
                        }
                        else if ( arc->type == CS )
                        {
                           cs = (GAS_CS*) arc->detailed_info;
                           vars[cnt] = cs->compressor_binvar;
                        }
                        else
                        {
                           SCIPerrorMessage("Unkown active element accross the cut.\n");
                           return SCIP_ERROR;
                        }

                        if ( flipped )
                           vals[cnt] = -arc->flowMin;
                        else
                           vals[cnt] = arc->flowMax;

                        /* possibly strengthen cut */
                        if ( vals[cnt] > rhs )
                           vals[cnt] = rhs;
                        ++cnt;
                     }
                     arc = arc->next_outarc;
                  }

                  /* loop over incoming arcs */
                  arc = node->inarcs;
                  while ( arc != NULL )
                  {
                     assert( SCIPisLE(scip, arc->flowMin, 0.0) );
                     assert( SCIPisGE(scip, arc->flowMax, 0.0) );

                     /* if arc is accross the cut */
                     if ( component[arc->sourcenode->nodeposition] != i )
                     {
                        if ( arc->type == VALVE )
                        {
                           valve = (GAS_Valve*) arc->detailed_info;
                           vars[cnt] = valve->valve_binvar;
                        }
                        else if ( arc->type == CONTROLVALVE )
                        {
                           cv = (GAS_Controlvalve*) arc->detailed_info;
                           vars[cnt] = cv->binvar;
                        }
                        else if ( arc->type == CS )
                        {
                           cs = (GAS_CS*) arc->detailed_info;
                           vars[cnt] = cs->compressor_binvar;
                        }
                        else
                        {
                           SCIPerrorMessage("Unkown active element accross the cut.\n");
                           return SCIP_ERROR;
                        }

                        if ( flipped )
                           vals[cnt] = arc->flowMax;
                        else
                           vals[cnt] = -arc->flowMin;

                        /* possibly strengthen cut */
                        if ( vals[cnt] > rhs )
                           vals[cnt] = rhs;
                        ++cnt;
                     }
                     arc = arc->next_inarc;
                  }
               }
            }

            if ( cnt > 0 )
            {
               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "componentcut#%d", i);
               SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, cnt, vars, vals, rhs, SCIPinfinity(scip)) );
               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
               ++ngen;
            }
         }
      }
      if ( nvisited == network->numnodes )
         break;
   }

   SCIPfreeBufferArray(scip, &component);
   SCIPfreeBufferArray(scip, &inqueue);
   SCIPfreeBufferArray(scip, &queue);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   SCIPinfoMessage(scip, NULL, "\nAdded %d components cut constraints.\n", ngen);

   return SCIP_OKAY;
}


/** auxiliary routine to check for only pipe arc excepting compressor and bypass arcs */
static
GAS_Arc* findOnlyPipeArc(
   GAS_Node*             w,                  /**< node to check for */
   GAS_Arc*              cs_arc,             /**< compressor arc that should be ignored */
   GAS_Arc*              bypass_arc          /**< bypass arc that should be ignored */
   )
{
   GAS_Arc* arc;
   GAS_Arc* pipe_arc;

   /* check incoming arcs */
   arc = w->inarcs;
   pipe_arc = NULL;
   while ( arc != NULL )
   {
      /* skip compressor or bypass arcs */
      if ( arc == cs_arc )
      {
         arc = arc->next_inarc;
         continue;
      }
      if ( arc == bypass_arc )
      {
         arc = arc->next_inarc;
         continue;
      }

      if ( arc->type != PIPE )
         break;
      if ( pipe_arc != NULL )
         break;
      pipe_arc = arc;
      arc = arc->next_inarc;
   }

   if ( arc != NULL )
      return NULL;

   if ( pipe_arc == NULL )
   {
      arc = w->outarcs;
      while ( arc != NULL )
      {
         /* skip compressor or bypass arcs */
         if ( arc == cs_arc )
         {
            arc = arc->next_outarc;
            continue;
         }
         if ( arc == bypass_arc )
         {
            arc = arc->next_outarc;
            continue;
         }

         if ( arc->type != PIPE )
            break;
         if ( pipe_arc != NULL )
            break;
         pipe_arc = arc;
         arc = arc->next_outarc;
      }
   }

   return pipe_arc;
}


/** generate symmetry breaking constraints for parallel compressors
 *
 *  If idential active constraints are parallel, we can add symmetry breaking constraints.
 *
 *  These inequalities will likely be detected automatically in the future.
 *
 *  They currently only seem to help for the symmetrized version of GasLib-135-sym.
 */
static
SCIP_RETCODE generateParallelCompressors(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GAS_Network*          network             /**< network */
   )
{
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   GAS_Arc** csarcs;
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   int i;

   SCIP_CALL( SCIPallocBufferArray(scip, &csarcs, network->numarcs) );
   vals[0] = 1.0;
   vals[1] = -1.0;

   /* loop through nodes */
   for (i = 0; i < network->numnodes; ++i)
   {
      GAS_Arc* arc;
      GAS_Node* v;
      int ncs = 0;

      v = &network->nodes_ptr[i];

      /* check neighbors for compressor stations */
      arc = v->inarcs;
      while ( arc != NULL )
      {
         if ( arc->type == CS )
            csarcs[ncs++] = arc;
         arc = arc->next_inarc;
      }

      /* check outarcs if no compressors have been found */
      if ( ncs == 0 )
      {
         arc = v->outarcs;
         while ( arc != NULL )
         {
            if ( arc->type == CS )
               csarcs[ncs++] = arc;
            arc = arc->next_outarc;
         }
      }

      /* check for pairs of compressor stations */
      if ( ncs >= 2 )
      {
         GAS_Arc* cs1_arc;
         GAS_Arc* cs2_arc;
         GAS_CS* cs1;
         GAS_CS* cs2;
         GAS_Arc* pipe1_arc;
         GAS_Arc* pipe2_arc;
         GAS_Pipe* pipe1;
         GAS_Pipe* pipe2;
         GAS_Node* w1;
         GAS_Node* w2;
         int k;
         int l;

         for (k = 0; k < ncs; ++k)
         {
            cs1_arc = csarcs[k];
            cs1 = (GAS_CS*) cs1_arc->detailed_info;
            if ( cs1->numPistonCS > 0 )
               continue;
            if ( cs1->numPistonCS > 0 )
               continue;

            /* check incoming arcs */
            w1 = cs1_arc->sourcenode;
            if ( w1 == v )
               w1 = cs1_arc->targetnode;
            assert( cs1_arc->sourcenode == v || cs1_arc->targetnode == v );

            pipe1_arc = findOnlyPipeArc(w1, cs1_arc, cs1->bypass);
            if ( pipe1_arc == NULL )
               continue;
            pipe1 = (GAS_Pipe*) pipe1_arc->detailed_info;

            for (l = k+1; l < ncs; ++l)
            {
               /* check whether compressor stations are idential */
               cs2_arc = csarcs[l];
               cs2 = (GAS_CS*) cs2_arc->detailed_info;
               if ( cs2->numPistonCS > 0 )
                  continue;

               if ( cs1_arc->flowMin != cs2_arc->flowMin ) /*lint !e777*/
                  continue;
               if ( cs1_arc->flowMax != cs2_arc->flowMax ) /*lint !e777*/
                  continue;

               if ( cs1->pressureLossIn != cs2->pressureLossIn ) /*lint !e777*/
                  continue;
               if ( cs1->pressureLossOut != cs2->pressureLossOut ) /*lint !e777*/
                  continue;
               if ( cs1->dragFactorIn != cs2->dragFactorIn ) /*lint !e777*/
                  continue;
               if ( cs1->dragFactorOut != cs2->dragFactorOut ) /*lint !e777*/
                  continue;

               /* the in diameter only is improtant if dragFactorIn is nonzero */
               if ( ! SCIPisZero(scip, cs1->dragFactorIn) )
               {
                  if ( cs1->diameterIn != cs2->diameterIn ) /*lint !e777*/
                     continue;
               }

               /* the out diameter only is improtant if dragFactorIn is nonzero */
               if ( ! SCIPisZero(scip, cs1->dragFactorOut) )
               {
                  if ( cs1->dragFactorOut != cs2->dragFactorOut ) /*lint !e777*/
                     continue;
               }

               if ( cs1->pressureInMin != cs2->pressureInMin ) /*lint !e777*/
                  continue;
               if ( cs1->pressureOutMax != cs2->pressureOutMax ) /*lint !e777*/
                  continue;
               if ( cs1->numTurboCS != cs2->numTurboCS )
                  continue;
               if ( cs1->numconfigurations != cs2->numconfigurations )
                  continue;
               if ( cs1->internalBypassRequired != cs2->internalBypassRequired )
                  continue;

               if ( probdata->boxConstraintModel )
               {
                  /* do additional checks about the equality of the box constraints ... */
               }

               /* check incoming arcs */
               w2 = cs2_arc->sourcenode;
               if ( w2 == v )
                  w2 = cs2_arc->targetnode;
               assert( cs2_arc->sourcenode == v || cs2_arc->targetnode == v );

               pipe2_arc = findOnlyPipeArc(w2, cs2_arc, cs2->bypass);
               if ( pipe2_arc == NULL )
                  continue;
               pipe2 = (GAS_Pipe*) pipe2_arc->detailed_info;

               /* check whehter the other node agrees */
               if ( pipe1_arc->sourcenode == w1 && pipe2_arc->sourcenode == w2 )
               {
                  if ( pipe1_arc->targetnode != pipe2_arc->targetnode )
                     continue;
               }
               if ( pipe1_arc->sourcenode == w1 && pipe2_arc->targetnode == w2 )
               {
                  if ( pipe1_arc->targetnode != pipe2_arc->sourcenode )
                     continue;
               }
               if ( pipe1_arc->targetnode == w1 && pipe2_arc->sourcenode == w2 )
               {
                  if ( pipe1_arc->sourcenode != pipe2_arc->targetnode )
                     continue;
               }
               if ( pipe1_arc->targetnode == w1 && pipe2_arc->targetnode == w2 )
               {
                  if ( pipe1_arc->sourcenode != pipe2_arc->sourcenode )
                     continue;
               }

               /* compare pipes */
               if ( pipe1_arc->flowMin != pipe2_arc->flowMin ) /*lint !e777*/
                  continue;
               if ( pipe1_arc->flowMax != pipe2_arc->flowMax ) /*lint !e777*/
                  continue;

               if ( pipe1->length != pipe2->length ) /*lint !e777*/
                  continue;
               if ( pipe1->diameter != pipe2->diameter ) /*lint !e777*/
                  continue;
               if ( pipe1->roughness != pipe2->roughness ) /*lint !e777*/
                  continue;
               if ( pipe1->pressureMax != pipe2->pressureMax ) /*lint !e777*/
                  continue;
               if ( pipe1->slope != pipe2->slope ) /*lint !e777*/
                  continue;

               vars[0] = cs1->compressor_binvar;
               vars[1] = cs2->compressor_binvar;

               SCIPinfoMessage(scip, NULL, "Adding symmetry breaking inequality for parallel compressors <%s> and <%s>.\n", csarcs[k]->id, csarcs[l]->id);
               (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "cssymbreak#%s#%s", cs1_arc->id, cs2_arc->id);
               SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 2, vars, vals, 0.0, SCIPinfinity(scip)) );
               SCIP_CALL( SCIPaddCons(scip, cons) );
               SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &csarcs);

   return SCIP_OKAY;
}


/** run DFS to compute articulation nodes for passive subnetworks */
static
void passiveArticulationPoints(
   GAS_Network*          network,            /**< network */
   GAS_Node*             v,                  /**< current node */
   int                   d,                  /**< depth */
   int*                  parent,             /**< partent nodes */
   int*                  depth,              /**< depth of nodes */
   int*                  low,                /**< lowest depth of decendents */
   SCIP_Bool*            passive,            /**< pointer to store whether the subnetwork is passive */
   SCIP_Bool*            acyclic,            /**< pointer to store whether the subnetwork is acyclic */
   SCIP_Bool*            allfixed,           /**< whether the flows on all arcs in the subnetwork are fixed */
   int*                  articulationnodes,  /**< array to store articulation nodes */
   GAS_Arc**             articulationarcs,   /**< array to store arcs that create components of articulation nodes */
   int*                  narticulationnodes  /**< pointer to store the number of articulation nodes */
   )
{
   GAS_Node* w;
   GAS_Arc* arc;
   SCIP_Bool subpassive;
   SCIP_Bool subacyclic;
   SCIP_Bool suballfixed;
   int nchildren = 0;
   SCIP_Bool marked = FALSE;

   assert( v != NULL );
   assert( 0 <= v->nodeposition && v->nodeposition < network->numnodes );

   depth[v->nodeposition] = d;
   low[v->nodeposition] = d;

   *passive = TRUE;
   *acyclic = TRUE;
   *allfixed = TRUE;

   /* check in-neighbors */
   arc = v->inarcs;
   while ( arc != NULL )
   {
      w = arc->sourcenode;
      assert( 0 <= w->nodeposition && w->nodeposition < network->numnodes );

      /* if any arc in the component is not active, the whole component is not passive */
      if ( arc->type != PIPE && arc->type != SHORTPIPE )
         *passive = FALSE;

      /* check whether any arc is not fixed */
      if ( REALABS(arc->flowMax - arc->flowMin) > 1e-9 )
         *allfixed = FALSE;

      /* if node has not been visited */
      if ( depth[w->nodeposition] < 0 )
      {
         parent[w->nodeposition] = v->nodeposition;
         passiveArticulationPoints(network, w, d + 1, parent, depth, low, &subpassive, &subacyclic, &suballfixed, articulationnodes, articulationarcs, narticulationnodes);
         *passive = *passive && subpassive;
         *acyclic = *acyclic && subacyclic;
         *allfixed = *allfixed && suballfixed;
         ++nchildren;

         assert( low[w->nodeposition] >= 0 );

         /* we have an articulation point if the low point is above the current node */
         if ( low[w->nodeposition] >= d && ! marked )
         {
            /* check whether we found a passive and cyclic subnetwork (acyclic subnetworks will be dealt with automatically) in which not all flows are fixed */
            if ( subpassive && ! subacyclic && ! suballfixed )
            {
               /* if the node is not the root and an articulation point or it is the root and has children */
               if ( parent[v->nodeposition] >= 0 || ( parent[v->nodeposition] < 0 && nchildren > 1 ) )
               {
                  articulationnodes[*narticulationnodes] = v->nodeposition;
                  articulationarcs[*narticulationnodes] = arc;
                  ++(*narticulationnodes);
                  assert( *narticulationnodes <= network->numnodes );
                  marked = TRUE;
               }
            }
         }

         /* update low point of current node to min(low[v], low[w]) */
         if ( low[w->nodeposition] < low[v->nodeposition] )
            low[v->nodeposition] = low[w->nodeposition];
      }
      else if ( w->nodeposition != parent[v->nodeposition] )
      {
         *acyclic = FALSE;

         /* update low point of current node to min(low[v], depth[w]) */
         if ( depth[w->nodeposition] < low[v->nodeposition] )
            low[v->nodeposition] = depth[w->nodeposition];
      }
      arc = arc->next_inarc;
   }

   /* check out-neighbors */
   arc = v->outarcs;
   while ( arc != NULL )
   {
      w = arc->targetnode;
      assert( 0 <= w->nodeposition && w->nodeposition < network->numnodes );

      /* if any arc in the component is not active, the whole component is not passive */
      if ( arc->type != PIPE && arc->type != SHORTPIPE )
         *passive = FALSE;

      /* check whether any are is not fixed */
      if ( REALABS(arc->flowMax - arc->flowMin) > 1e-9 )
         *allfixed = FALSE;

      /* if node has not been visited */
      if ( depth[w->nodeposition] < 0 )
      {
         parent[w->nodeposition] = v->nodeposition;
         passiveArticulationPoints(network, w, d + 1, parent, depth, low, &subpassive, &subacyclic, &suballfixed, articulationnodes, articulationarcs, narticulationnodes);
         *passive = *passive && subpassive;
         *acyclic = *acyclic && subacyclic;
         *allfixed = *allfixed && suballfixed;
         ++nchildren;

         assert( low[w->nodeposition] >= 0 );

         /* we have an articulation point if the low point is above the current node */
         if ( low[w->nodeposition] >= d && ! marked )
         {
            /* check whether we found a passive and cyclic subnetwork (acyclic subnetworks will be dealt with automatically) in which not all flows are fixed */
            if ( subpassive && ! subacyclic && ! suballfixed )
            {
               /* if the node is not the root and an articulation point or it is the root and has children */
               if ( parent[v->nodeposition] >= 0 || ( parent[v->nodeposition] < 0 && nchildren > 1 ) )
               {
                  articulationnodes[*narticulationnodes] = v->nodeposition;
                  articulationarcs[*narticulationnodes] = arc;
                  ++(*narticulationnodes);
                  assert( *narticulationnodes <= network->numnodes );
                  marked = TRUE;
               }
            }
         }

         /* update low point of current node to min(low[v], low[w]) */
         if ( low[w->nodeposition] < low[v->nodeposition] )
            low[v->nodeposition] = low[w->nodeposition];
      }
      else if ( parent[v->nodeposition] != w->nodeposition )
      {
         *acyclic = FALSE;

         /* update low point of current node to min(low[v], depth[w]) */
         if ( depth[w->nodeposition] < low[v->nodeposition] )
            low[v->nodeposition] = depth[w->nodeposition];
      }
      arc = arc->next_outarc;
   }
}


/** create passive subnetwork */
static
SCIP_RETCODE createSubnetwork(
   SCIP*                 scip,               /**< SCIP data structure */
   GAS_Network*          network,            /**< network */
   GAS_Network**         subnetworkptr,      /**< subnetwork */
   GAS_Node*             root,               /**< articulation point */
   GAS_Arc*              rootarc,            /**< arc creating passive subnetwork */
   int*                  arc_subnettonet     /**< map from subnet arcs to original arcs */
   )
{
   GAS_Network* subnetwork;
   GAS_Node** nodequeue;
   GAS_Node* firstnode;
   SCIP_Real balance;
   int* node_nettosubnet;
   int* node_subnettonet;
   int* arc_nettosubnet;
   int lqueue;
   int uqueue;
   int nsubnodes = 0;
   int nsubarcs = 0;
   int nshortpipes = 0;
   int i;
   int j;

   assert( rootarc->type == PIPE || rootarc->type == SHORTPIPE );

   SCIP_CALL( SCIPallocBufferArray(scip, &nodequeue, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &node_nettosubnet, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &node_subnettonet, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &arc_nettosubnet, network->numarcs) );

   for (i = 0; i < network->numnodes; ++i)
   {
      node_nettosubnet[i] = -1;
      node_subnettonet[i] = -1;
   }

   for (j = 0; j < network->numarcs; ++j)
   {
      arc_nettosubnet[j] = -1;
      arc_subnettonet[j] = -1;
   }

   /* treat root */
   node_nettosubnet[root->nodeposition] = nsubnodes;
   node_subnettonet[nsubnodes] = root->nodeposition;
   ++nsubnodes;

   /* init queue */
   if ( rootarc->targetnode != root )
   {
      assert( rootarc->sourcenode == root );
      firstnode = rootarc->targetnode;
   }
   else
   {
      assert( rootarc->targetnode == root );
      firstnode = rootarc->sourcenode;
   }

   assert( rootarc->type == PIPE || rootarc->type == SHORTPIPE );
   assert( arc_nettosubnet[rootarc->arcposition] < 0 );
   if ( rootarc->type == SHORTPIPE )
      ++nshortpipes;
   arc_nettosubnet[rootarc->arcposition] = nsubarcs;
   arc_subnettonet[nsubarcs] = rootarc->arcposition;
   ++nsubarcs;

   assert( node_nettosubnet[firstnode->nodeposition] < 0 );
   node_nettosubnet[firstnode->nodeposition] = nsubnodes;
   node_subnettonet[nsubnodes] = firstnode->nodeposition;
   ++nsubnodes;
   balance = firstnode->flow;

   /* init queue */
   nodequeue[0] = firstnode;
   lqueue = 0;
   uqueue = 1;

   /* fill in mapping functions */
   while ( uqueue > lqueue )
   {
      GAS_Arc* arc;
      GAS_Node* v;
      GAS_Node* w;

      v = nodequeue[lqueue++];
      assert( lqueue <= network->numnodes );

      arc = v->inarcs;
      while ( arc != NULL )
      {
         assert( arc->type == PIPE || arc->type == SHORTPIPE );
         if ( arc_nettosubnet[arc->arcposition] < 0 )
         {
            if ( arc->type == SHORTPIPE )
               ++nshortpipes;
            arc_nettosubnet[arc->arcposition] = nsubarcs;
            arc_subnettonet[nsubarcs] = arc->arcposition;
            ++nsubarcs;
         }

         w = arc->sourcenode;
         if ( node_nettosubnet[w->nodeposition] < 0 )
         {
            node_nettosubnet[w->nodeposition] = nsubnodes;
            node_subnettonet[nsubnodes] = w->nodeposition;
            ++nsubnodes;

            nodequeue[uqueue++] = w;
            assert( uqueue <= network->numnodes );

            if ( w != root )
               balance += w->flow;
         }
         arc = arc->next_inarc;
      }

      arc = v->outarcs;
      while ( arc != NULL )
      {
         assert( arc->type == PIPE || arc->type == SHORTPIPE );
         if ( arc_nettosubnet[arc->arcposition] < 0 )
         {
            if ( arc->type == SHORTPIPE )
               ++nshortpipes;
            arc_nettosubnet[arc->arcposition] = nsubarcs;
            arc_subnettonet[nsubarcs] = arc->arcposition;
            ++nsubarcs;
         }

         w = arc->targetnode;
         if ( node_nettosubnet[w->nodeposition] < 0 )
         {
            node_nettosubnet[w->nodeposition] = nsubnodes;
            node_subnettonet[nsubnodes] = w->nodeposition;
            ++nsubnodes;

            nodequeue[uqueue++] = w;
            assert( uqueue <= network->numnodes );

            if ( w != root )
               balance += w->flow;
         }
         arc = arc->next_outarc;
      }
   }

   /* create subnetwork */
   SCIP_CALL( SCIPallocBlockMemory(scip, &subnetwork) );
   *subnetworkptr = subnetwork;

   subnetwork->numnodes                  = nsubnodes;
   subnetwork->numarcs                   = nsubarcs;
   subnetwork->maxnumincidentarcs        = 0;
   subnetwork->numflowdirvars            = 0;
   subnetwork->numpipes                  = nsubarcs - nshortpipes;
   subnetwork->numshortpipes             = nshortpipes;
   subnetwork->numvalves                 = 0;
   subnetwork->numcontrolvalves          = 0;
   subnetwork->numcompressor             = 0;
   subnetwork->numcsconfigurations       = 0;
   subnetwork->numcsvarsfornlrs          = 0;
   subnetwork->numnonlinresistor         = 0;
   subnetwork->numlinresistor            = 0;
   subnetwork->numpivars                 = 0;
   subnetwork->normDensity               = network->normDensity;
   subnetwork->pseudocriticalPressure    = network->pseudocriticalPressure;
   subnetwork->pseudocriticalTemperature = network->pseudocriticalTemperature;
   subnetwork->gasTemperature            = network->gasTemperature;
   subnetwork->molarMass1                = network->molarMass1;
   subnetwork->molarMass2                = network->molarMass2;
   /* subnetwork->specificHeatCap1 = network->specificHeatCap1; */
   /* subnetwork->specificHeatCap2 = network->specificHeatCap1; */
   subnetwork->nodes_ptr                 = NULL;
   subnetwork->arcs_ptr                  = NULL;
   subnetwork->spanningTree              = NULL;
   subnetwork->numBasisCycles            = 0;
   subnetwork->firstCycle                = NULL;
   subnetwork->lastCycle                 = NULL;
   subnetwork->firstCombinedCycle        = NULL;

   /* allocates memory for subnodes and subarcs */
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &subnetwork->nodes_ptr, nsubnodes) );
   SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &subnetwork->arcs_ptr, nsubarcs) );

   /* fill in sunetwork nodes */
   for (i = 0; i < nsubnodes; ++i)
   {
      GAS_Node* subnode;
      GAS_Node* node;
      GAS_Arc* arc;

      assert( 0 <= node_subnettonet[i] && node_subnettonet[i] < network->numnodes );
      node = &(network->nodes_ptr[node_subnettonet[i]]);

      subnode = &(subnetwork->nodes_ptr[i]);
      subnode->outarcs = NULL;
      subnode->inarcs = NULL;

      /* determine inarcs: find first arc that is present in the subnetwork */
      arc = node->inarcs;
      while ( arc != NULL )
      {
         if ( arc_nettosubnet[arc->arcposition] >= 0 )
         {
            subnode->inarcs = &(subnetwork->arcs_ptr[arc_nettosubnet[arc->arcposition]]);
            break;
         }
         arc = arc->next_inarc;
      }

      /* determine outarcs: find first arc that is present in the subnetwork */
      arc = node->outarcs;
      while ( arc != NULL )
      {
         if ( arc_nettosubnet[arc->arcposition] >= 0 )
         {
            subnode->outarcs = &(subnetwork->arcs_ptr[arc_nettosubnet[arc->arcposition]]);
            break;
         }
         arc = arc->next_outarc;
      }

      /* set/copy values */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(subnode->id), node->id, strlen(node->id)+1) );
      subnode->nodeposition    = i;
      subnode->geoWGS84Long    = node->geoWGS84Long;
      subnode->geoWGS84Lat     = node->geoWGS84Lat;
      subnode->x               = node->x;
      subnode->y               = node->y;
      subnode->flow            = node->flow;
      subnode->scnFlowMax      = node->scnFlowMax;
      subnode->scnFlowMin      = node->scnFlowMin;
      subnode->pressureMin     = node->pressureMin;
      subnode->pressureMax     = node->pressureMax;
      subnode->netPressureMin  = node->netPressureMin;
      subnode->netPressureMax  = node->netPressureMax;
      subnode->scnPressureMin  = node->scnPressureMin;
      subnode->scnPressureMax  = node->scnPressureMax;
      subnode->arcPressureMax  = node->arcPressureMax;
      subnode->numneighbors    = 0;
      subnode->numincidentarcs = 0;
      subnode->pressurevar     = NULL;
      subnode->needPIvar       = TRUE;
      subnode->PIvar           = NULL;
      subnode->pSlackVar       = NULL;
      subnode->type            = node->type;
      subnode->nodeMixingRatio = node->nodeMixingRatio;

      ++subnetwork->numpivars;
   }

   /* adjust balance of root */
   {
      GAS_Node* subroot;

      assert( 0 <= node_nettosubnet[root->nodeposition] && node_nettosubnet[root->nodeposition] < network->numnodes );
      subroot = &(subnetwork->nodes_ptr[node_nettosubnet[root->nodeposition]]);
      subroot->flow = -balance;
   }

   /* fill in subnetwork arcs */
   for (j = 0; j < nsubarcs; ++j)
   {
      GAS_Arc* subarc;
      GAS_Arc* arc;
      GAS_Arc* arc2;
      GAS_Pipe* pipe;
      GAS_Pipe* subpipe;

      assert( 0 <= arc_subnettonet[j] && arc_subnettonet[j] < network->numarcs );
      arc = &(network->arcs_ptr[arc_subnettonet[j]]);
      subarc = &(subnetwork->arcs_ptr[j]);

      /* set/copy values */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(subarc->id), arc->id, strlen(arc->id)+1) );
      subarc = &(subnetwork->arcs_ptr[j]);
      subarc->sourcenode = &(subnetwork->nodes_ptr[node_nettosubnet[arc->sourcenode->nodeposition]]);
      subarc->targetnode = &(subnetwork->nodes_ptr[node_nettosubnet[arc->targetnode->nodeposition]]);
      ++subarc->sourcenode->numincidentarcs;
      ++subarc->targetnode->numincidentarcs;

      subarc->arcposition= j;
      subarc->type = arc->type;
      subarc->flowMax = arc->flowMax;
      subarc->flowMin = arc->flowMin;
      subarc->flowvar = NULL;
      subarc->positiveFlowBinvar = NULL;
      subarc->negativeFlowBinvar = NULL;
      subarc->next_inarc = NULL;
      subarc->next_outarc = NULL;
      subarc->detailed_info = NULL;

      assert( arc->type == PIPE || arc->type == SHORTPIPE );

      if ( arc->type == PIPE )
      {
         pipe = (GAS_Pipe*) arc->detailed_info;
         SCIP_CALL( SCIPallocBlockMemory(scip, &subpipe) );
         subarc->detailed_info = subpipe;

         subpipe->length      = pipe->length;
         subpipe->diameter    = pipe->diameter;
         subpipe->roughness   = pipe->roughness;
         subpipe->pressureMax = pipe->pressureMax;
         subpipe->slope       = pipe->slope;
         subpipe->PIdiffVar   = NULL;
      }
      /* nothing to do for short pipes */

      /* find first arc that is present in the subnetwork */
      arc2 = arc->next_inarc;
      while ( arc2 != NULL )
      {
         if ( arc_nettosubnet[arc2->arcposition] >= 0 )
         {
            subarc->next_inarc = &(subnetwork->arcs_ptr[arc_nettosubnet[arc2->arcposition]]);;
            break;
         }
         arc2 = arc2->next_inarc;
      }

      /* find first arc that is present in the subnetwork */
      arc2 = arc->next_outarc;
      while ( arc2 != NULL )
      {
         if ( arc_nettosubnet[arc2->arcposition] >= 0 )
         {
            subarc->next_outarc = &(subnetwork->arcs_ptr[arc_nettosubnet[arc2->arcposition]]);;
            break;
         }
         arc2 = arc2->next_outarc;
      }
   }

   /* determine maximal number of incident arcs */
   for (i = 0; i < nsubnodes; ++i)
   {
      GAS_Node* node;

      node = &(subnetwork->nodes_ptr[i]);
      if ( node->numincidentarcs > subnetwork->maxnumincidentarcs )
         subnetwork->maxnumincidentarcs = node->numincidentarcs;
   }

   SCIPfreeBufferArray(scip, &arc_nettosubnet);
   SCIPfreeBufferArray(scip, &node_subnettonet);
   SCIPfreeBufferArray(scip, &node_nettosubnet);
   SCIPfreeBufferArray(scip, &nodequeue);

   SCIPinfoMessage(scip, NULL, "Found subnetwork: root = %s (balance = %f), #nodes = %d, #pipes = %d, #short pipes = %d.\n", root->id, balance, nsubnodes, nsubarcs, nshortpipes);

   return SCIP_OKAY;
}


/** try to precompute the flows in passive components
 *
 *  We search for passive components in the network using articulation points. For such components one can precompute
 *  the flow, which can then be fixed.
 *
 *  Note that this procedure only seems to be valid for the algebraic model, where the pipe parameters are constant and
 *  do not depend on the flow or pressure.
 */
static
SCIP_RETCODE preprocessPassiveComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   GAS_Network*          network             /**< network */
   )
{
   GAS_Node* v;
   int* depth;
   int* low;
   int* parent;
   int* articulationnodes;
   GAS_Arc** articulationarcs;
   int narticulationnodes = 0;
   int* arc_subnettonet;
   SCIP_Bool passive;
   SCIP_Bool acyclic;
   SCIP_Bool allfixed;
   int i;

   assert( probdata->algebraic );

   /* do not run recursively - this is the case if no network name is given */
   if ( probdata->netname == NULL )
      return SCIP_OKAY;

   SCIPinfoMessage(scip, NULL, "\nSearching for passive subnetworks ...\n");

   SCIP_CALL( SCIPallocBufferArray(scip, &depth, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &parent, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &low, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &articulationnodes, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &articulationarcs, network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &arc_subnettonet, network->numarcs) );

   for (i = 0; i < network->numnodes; ++i)
   {
      depth[i] = -1;
      parent[i] = -1;
      low[i] = -1;
   }

   /* search for articulation points with passive and cyclic subnetworks */
   v = &network->nodes_ptr[0];
   passiveArticulationPoints(network, v, 0, parent, depth, low, &passive, &acyclic, &allfixed, articulationnodes, articulationarcs, &narticulationnodes);

   /* loop through articulation points */
   if ( narticulationnodes > 0 )
   {
      for (i = 0; i < narticulationnodes; ++i)
      {
         SCIP* subscip;
         GAS_Network* subnetwork;
         SCIP_PROBDATA* subprobdata;
         SCIP_SOL* bestsol;
         int j;

         v = &network->nodes_ptr[articulationnodes[i]];

         /* initialize subscip */
         SCIP_CALL( SCIPcreate(&subscip) );

         /* include default SCIP plugins */
         SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

         /* create subnetwork */
         SCIP_CALL( createSubnetwork(subscip, network, &subnetwork, v, articulationarcs[i], arc_subnettonet) );

         /* create probdata */
         SCIP_CALL( SCIPallocBlockMemory(subscip, &subprobdata) );
         subprobdata->network             = subnetwork;
         subprobdata->netname             = NULL;
         subprobdata->scnname             = NULL;
         subprobdata->objective           = NULL;
         subprobdata->flowvars            = NULL;
         subprobdata->netflowvars         = NULL;
         subprobdata->positiveFlowBinvars = NULL;
         subprobdata->negativeFlowBinvars = NULL;
         subprobdata->pressurevars        = NULL;
         subprobdata->PIvars              = NULL;
         subprobdata->pSlackVars          = NULL;
         subprobdata->PIdiffvars          = NULL;
         subprobdata->VALVE_binvars       = NULL;
         subprobdata->CV_binvars          = NULL;
         subprobdata->NLR_DeltaVars       = NULL;
         subprobdata->NLR_AbsDeltaVars    = NULL;
         subprobdata->NLR_AbsFlowVars     = NULL;
         subprobdata->LR_smoothingFlow    = NULL;
         subprobdata->LR_posFlowDir       = NULL;
         subprobdata->LR_negFlowDir       = NULL;
         subprobdata->CS_binvars          = NULL;
         subprobdata->CS_resistorvars     = NULL;
         subprobdata->BCM_flowCS          = NULL;
         subprobdata->BCM_pressureIn      = NULL;
         subprobdata->BCM_pressureOut     = NULL;
         subprobdata->Lambda              = NULL;
         subprobdata->mixingRatio         = NULL;
         subprobdata->mixingRatioMass     = NULL;
         subprobdata->nodeMixingRatio     = NULL;
         subprobdata-> nodeMixingRatioMol = NULL;

         /* model options */
         subprobdata->noFlowBinvars       = TRUE; /* probdata->noFlowBinvars;*/
         subprobdata->binaryFlowCons      = FALSE; /* probdata->binaryFlowCons; */
         subprobdata->noCycles            = probdata->noCycles;
         subprobdata->allCycles           = probdata->allCycles;
         subprobdata->algebraic           = probdata->algebraic;
         subprobdata->binVarEQone         = probdata->binVarEQone;
         subprobdata->papay               = probdata->papay;
         subprobdata->boxConstraintModel  = probdata->boxConstraintModel;
         subprobdata->additionalFacets    = probdata->additionalFacets;
         subprobdata->charDiagram         = probdata->charDiagram;
         subprobdata->relaxLowerBounds    = probdata->relaxLowerBounds;
         subprobdata->relaxUpperBounds    = probdata->relaxUpperBounds;
         subprobdata->relaxCVbounds       = probdata->relaxCVbounds;
         subprobdata->relaxLBvalue        = probdata->relaxLBvalue;
         subprobdata->relaxUBvalue        = probdata->relaxUBvalue;
         subprobdata->relaxCVvalue        = probdata->relaxCVvalue;
         subprobdata->resistorCVInternal  = probdata->resistorCVInternal;
         subprobdata->idealCSminIncrease  = probdata->idealCSminIncrease;
         subprobdata->idealCSmaxIncrease  = probdata->idealCSmaxIncrease;
         subprobdata->csusemaxincrease    = probdata->csusemaxincrease;
         subprobdata->addComponentCuts    = FALSE;
         subprobdata->addParallelCSCuts   = FALSE;
         subprobdata->reduceFlowConservation = FALSE;
         subprobdata->mixing              = FALSE; /* mixing can be turned of since only flow is fixed */
         subprobdata->linearMixing        = FALSE;
         subprobdata->nodeMu              = probdata->nodeMu;
         subprobdata->mol                 = probdata->mol;
         subprobdata->flowConsMixing      = probdata->flowConsMixing;
         subprobdata->approx              = FALSE;
         subprobdata->aggregateParallelFlowdirVars = probdata->aggregateParallelFlowdirVars;
         subprobdata->preprocessPassiveComponents = FALSE;

         /* objective functions */
         subprobdata->minCompSum          = FALSE;
         subprobdata->minPressSum         = FALSE;
         subprobdata->noObjective         = TRUE;
         subprobdata->powerLoss           = FALSE;
         subprobdata->minSlack            = FALSE;
         subprobdata->minSlackPerBound    = FALSE;
         subprobdata->maxFlow             = FALSE;
         subprobdata->minLambda           = FALSE;
         subprobdata->addComponentCuts    = FALSE;
         subprobdata->addParallelCSCuts   = FALSE;

         /* create problem */
         SCIP_CALL( GAScreateProbNetwork(subscip, "subnet", subprobdata, subnetwork) );

#ifdef SCIP_DEBUG
         SCIP_CALL( GASwriteNetwork(subscip, subprobdata, "subnet.net", FALSE) );
#endif

         /* generate gas problem */
         SCIP_CALL( GASgenerateModel(subscip) );

#ifdef SCIP_DEBUG
         SCIP_CALL( SCIPwriteOrigProblem(subscip, "subnetwork.cip", "cip", FALSE) );
#endif

         /* turn of output */
         SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

         /* solve problem */
         SCIP_CALL( SCIPsolve(subscip) );

         assert( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL );

         bestsol = SCIPgetBestSol(subscip);
         assert( bestsol != NULL );

         /* fix flows */
         for (j = 0; j < subnetwork->numarcs; ++j)
         {
            GAS_Arc* arc;
            GAS_Arc* subarc;
            SCIP_Real flowval;

            subarc = &(subnetwork->arcs_ptr[j]);
            assert( subarc->flowvar != NULL );
            flowval = SCIPgetSolVal(subscip, bestsol, subarc->flowvar);
            arc = &(network->arcs_ptr[arc_subnettonet[subarc->arcposition]]);
            if ( SCIPisGT(scip, flowval, SCIPvarGetLbLocal(arc->flowvar)) || SCIPisLT(scip, flowval, SCIPvarGetUbLocal(arc->flowvar)) )
            {
#if 0
	       SCIPinfoMessage(scip, NULL, "Fixing flow value of arc <%s> to %g.\n", arc->id, flowval);
#endif
               SCIP_CALL( SCIPchgVarLb(scip, arc->flowvar, flowval) );
               SCIP_CALL( SCIPchgVarUb(scip, arc->flowvar, flowval) );
            }
            else
	    {
#if 0
               SCIPinfoMessage(scip, NULL, "Already fixed flow value of arc <%s> = %g.\n", arc->id, flowval);
#endif
	    }
            /* also change bounds in network itself to avoid conflicting recursive calls */
            assert( SCIPisGE(scip, flowval, arc->flowMin) && SCIPisLE(scip, flowval, arc->flowMax) );
            arc->flowMin = flowval;
            arc->flowMax = flowval;
         }

         SCIP_CALL( SCIPfree(&subscip) );
      }
   }
   else
      SCIPinfoMessage(scip, NULL, "Found no passive articulation points.\n");

   SCIPfreeBufferArray(scip, &arc_subnettonet);
   SCIPfreeBufferArray(scip, &articulationarcs);
   SCIPfreeBufferArray(scip, &articulationnodes);
   SCIPfreeBufferArray(scip, &low);
   SCIPfreeBufferArray(scip, &parent);
   SCIPfreeBufferArray(scip, &depth);

   return SCIP_OKAY;
}


/** generate aggregation equation for flow direction variables of parallel arcs */
static
SCIP_RETCODE generateAggregationParallelArcs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   GAS_Network* network;
   char consname[SCIP_MAXSTRLEN];
   SCIP_VAR* vars[2];
   SCIP_Real vals[2];
   SCIP_CONS* cons;
   int ngen = 0;
   int i;

   /* exit if there are no binary flow variables */
   if ( probdata->noFlowBinvars )
      return SCIP_OKAY;

   network = probdata->network;
   assert( network != NULL );

   vals[0] = 1.0;
   vals[1] = -1.0;

   /* loop through nodes */
   for (i = 0; i < network->numnodes; ++i)
   {
      GAS_Arc* arc1;
      GAS_Node* v;

      v = &network->nodes_ptr[i];
      arc1 = v->outarcs;
      while ( arc1 != NULL )
      {
         GAS_Arc* arc2;
         GAS_Node* w;

         /* there are no flow direction variables for controlvalves and compressor stations */
         if ( arc1->type == CS || arc1->type == CONTROLVALVE )
         {
            arc1 = arc1->next_outarc;
            continue;
         }

         assert( arc1->positiveFlowBinvar != NULL );

         w = arc1->targetnode;

         /* search for parallel arcs to w */
         arc2 = arc1;
         arc2 = arc2->next_outarc;

         /* TODO: This code assumes that the out-neighbors are sorted - ensure this! */
         while ( arc2 != NULL && arc2->targetnode == w )
         {
            /* there are no flow direction variables for controlvalves and compressor stations */
            if ( arc2->type == CS || arc2->type == CONTROLVALVE )
            {
               arc2 = arc2->next_outarc;
               continue;;
            }

            assert( arc2->positiveFlowBinvar != NULL );

            vars[0] = arc1->positiveFlowBinvar;
            vars[1] = arc2->positiveFlowBinvar;

            (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "posflowbinvarequal#%s#%s", arc1->id, arc2->id);
            SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 2, vars, vals, 0.0, 0.0) );
#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            ++ngen;

            arc2 = arc2->next_outarc;
         }

         /* continue with search at next arc after target node w */
         arc1 = arc2;
      }
   }
   SCIPinfoMessage(scip, NULL, "\nAdded %d flow direction aggregation constraints.\n", ngen);

   return SCIP_OKAY;
}


/** create stationary gas transport model */
SCIP_RETCODE GASgenerateModel(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   GAS_Resistor* resistor;
   int valveposition = 0;
   int linresistorposition = 0;
   int j;

   /* get problem data */
   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );
   assert( probdata->network != NULL );

   /* generate variables */
   SCIP_CALL( generateVariables(scip) );

   /* link objective variable to the corresponding constraints */
   if ( probdata->minSlack || probdata->minSlackPerBound )
   {
      SCIP_CALL( generateObjectiveConstraints(scip) );
   }

   /* create the flow conservation constraints */
   SCIP_CALL( generateFlowConstraint(scip, probdata) );

   if ( ! probdata->noFlowBinvars )
   {
      if ( probdata->binaryFlowCons )
      {
         /* create binary flow conservation constraints */
         SCIP_CALL( generateBinaryFlowConservationConstraints(scip, probdata) );
      }

      /* create the coupling of the flowDirectionBinvars with the flow */
      SCIP_CALL( generateFlowDirectionConstraints(scip, probdata) );
   }

   /* generate aggregations between flow direction variables for parallel arcs */
   if ( probdata->aggregateParallelFlowdirVars )
   {
      SCIP_CALL( generateAggregationParallelArcs(scip, probdata) );
   }

   /* create cuts that forbid flow in a cycle */
   if ( probdata->network->firstCycle != NULL )
   {
      SCIP_CALL( generateNoCircularFlowConstraint(scip, probdata) );
   }

   /* create the constraints for the PIvars */
   if ( probdata->network->numpivars != 0 )
   {
      SCIP_CALL( generatePressureSquaredConstraints(scip, probdata) );
   }

   /* create mixing constraints */
   if ( probdata->mixing || probdata->linearMixing || probdata->approx )
   {
      SCIP_CALL( adjustMixigRatioLowerBound(scip, probdata) );

      if ( probdata->nodeMu )
      {
         SCIP_CALL( generateMixingRatioNode(scip, probdata) );
      }
      else if ( probdata->mol )
      {
         SCIP_CALL( generateMixingRatioMol(scip, probdata) );
      }
      else if ( probdata->flowConsMixingNode )
      {
         SCIP_CALL ( generateFlowConsMixingNode(scip, probdata) );
      }
      else if ( probdata->flowConsMixing )
      {
         SCIP_CALL ( generateFlowConsMixing(scip, probdata) );
      }
      else
      {
         SCIP_CALL( generateMixingRatio(scip, probdata) );
      }
   }

   /* generate constraints on the arcs */
   for (j = 0; j < probdata->network->numarcs; j++)
   {
      if ( probdata->network->arcs_ptr[j].type == PIPE )
      {
         if ( ! probdata->algebraic )
         {
            SCIP_CALL( generateODEModel(scip, probdata, j) );
         }
         else
         {
            /* in both functions it will be checked if the normal algebraic model or the mixing should be used */
            if ( probdata->nodeMu )
            {
               SCIP_CALL( generateAlgebraicModelNodeMu(scip, j) );
            }
            else if ( probdata->flowConsMixingNode )
            {
               SCIP_CALL( generateAlgebraicModelNodeMu(scip, j) );
            }
            else
            {
               SCIP_CALL( generateAlgebraicModel(scip, j) );
            }
         }
      }
      else if ( probdata->network->arcs_ptr[j].type == SHORTPIPE )
      {
         SCIP_CALL( generateShortPipeConstraints(scip, probdata, j) );
      }
      else if ( probdata->network->arcs_ptr[j].type == VALVE )
      {
         SCIP_CALL( generateValveConstraints(scip, probdata, j, valveposition) );
         ++valveposition;
      }
      else if ( probdata->network->arcs_ptr[j].type == CONTROLVALVE )
      {
         if ( probdata->resistorCVInternal )
         {
            SCIP_CALL( generateControlValveConstraintsInternalResistors(scip, probdata, &(probdata->network->arcs_ptr[j])) );
         }
         else
         {
            SCIP_CALL( generateControlValveConstraintsExternalResistors(scip, probdata, &(probdata->network->arcs_ptr[j])) );
         }
      }
      else if ( probdata->network->arcs_ptr[j].type == CS )
      {
         SCIP_CALL( generateCompressorStationConstraints(scip, probdata, &(probdata->network->arcs_ptr[j])) );
      }
      else if ( probdata->network->arcs_ptr[j].type == RESISTOR )
      {
         resistor = probdata->network->arcs_ptr[j].detailed_info;
         if ( resistor->type == NONLINEAR )
         {
            SCIP_CALL( generateNonlinearResistorConstraints(scip, probdata, j) );
         }
         else if ( resistor->type == LINEAR )
         {
            SCIP_CALL( generateLinearResistorConstraints(scip, probdata, j, linresistorposition) );
            ++linresistorposition;
         }
         else
         {
            SCIPerrorMessage("Unknown resistor type. No constraint created.\n");
         }
      }
   }

   if ( probdata->addComponentCuts )
   {
      SCIP_CALL( generateComponentCuts(scip, probdata->network) );
   }

   if ( probdata->addParallelCSCuts )
   {
      SCIP_CALL( generateParallelCompressors(scip, probdata, probdata->network) );
   }

   if ( probdata->preprocessPassiveComponents && probdata->algebraic && ! probdata->maxFlow )
   {
      SCIP_CALL( preprocessPassiveComponents(scip, probdata, probdata->network) );
   }

   return SCIP_OKAY;
}
