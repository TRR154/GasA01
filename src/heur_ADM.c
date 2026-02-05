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

/**@file   heur_ADM.c
 * @ingroup PRIMALHEURISTICS
 * @brief  LNS heuristic that is adapted to gas transportation problems
 * @author Pascal Börner
 */

 //#define SCIP_MORE_DEBUG
//#define VISUALIZE
//#define NORMAL_TERMINATION
//#define STOP_IF_PRIMAL_SOL_FOUND
//#define STALLNODE_LIMIT
//#define STOP_IF_SUB_SCIP_FOUND_SOL
#define STOP_IF_SUB_SCIP_SMALL_GAP
//#define PRINT_ARC_IDS_WITH_CHANGED_FLOW_DIRECTION
#define NLP_TL
//#define NLP_GAP
#define NLP_STALLNODE
#define NLP_SOL
//#define NDEBUG
//#define CHG_OBJ

#include "heur_ADM.h"
#include "LSFiles_gas.h"
#include <scip/scipdefplugins.h>
#include "generateModel_gas.h"
#include "characteristics_gas.h"
#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/cons_varbound.h"
#include "cons_pipe_ode.h"
#include "probdata_gas.h"
#include <string.h>
#include "heur_gas.h"

/* default values for standard parameters that every primal heuristic has in SCIP */
#define HEUR_NAME             "ADM"
#define HEUR_DESC             "ADM adapated to gas transport"
#define HEUR_DISPCHAR         '~'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0 // not needed
#define HEUR_MAXDEPTH         0 // only root node | 1: in every node
#define HEUR_TIMING           SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */
#define DEFAULT_NODELIMIT     SCIP_LONGINT_MAX       /**< maximal number of nodes to use in sub-SCIP */
#define DEFAULT_TRANSFERSOLS  FALSE     /**< Transfer solutions of main SCIP to sub-SCIP? */
#define DEFAULT_ALREADYRUN    FALSE     /**< Was the heur already run? if yes do not run after restart of scip */


/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          nodelimit;          /**< maximal number of nodes to use in sub-SCIP */
   SCIP_Bool             transfersols;       /**< transfer solutions of main SCIP to sub-SCIP */
   SCIP_Bool             alreadyRun;         /**< was the heuristic already run? */
};

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeADM)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( scip != NULL );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free heuristic data */
   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/*
1. get copy or original problem
2. Modify problem: fix all unfixed mu variables to some fix mu and delete the mu calculation constraints.
3. Solve the modifies problem (no mixing)
4. Determine a flow topological ordering of the nodes and calculate the mixing ratio according to flow of solution.
5. Fix mu in the original problem to the calculated mu and reiterate.
*/

/** Given the flow directions this functions compotes a topological ordering of the nodes */
static
void topSort(
   GAS_Network*          network,            /**< network */
   GAS_Node*             v,                  /**< current node */
   int*                  visited,            /**< visited flag (use -1 for unvisited) */
   GAS_Node**            nodeTopo,           /**< array to store topological ordering */
   int*                  index,              /**< current index in nodeTopo (starts at numnodes - 1) */
   SCIP_SOL*             sol,                /**< solution from which the flow directions are extracted */
   SCIP*                 subscip,            /**< origin of the above solution */
   SCIP_HASHMAP*         varmapfw            /**< hasmap to get the corresponding variable in the subscip */
   )
{
   GAS_Node* w;
   GAS_Arc* arc;
   SCIP_VAR* fbv;

   assert( v != NULL );
   assert( 0 <= v->nodeposition && v->nodeposition < network->numnodes );

   visited[v->nodeposition] = 1; /* mark as visited */

   /* visit all unvisited neighbors  */
   arc = v->outarcs;
   while ( arc != NULL )
   {
      fbv = SCIPhashmapGetImage( varmapfw, arc->positiveFlowBinvar );
      if ( SCIPgetSolVal(subscip, sol, fbv) > 0 )
      {
         w = arc->targetnode;
         assert(0 <= w->nodeposition && w->nodeposition < network->numnodes);

         if ( visited[w->nodeposition] < 0 )
            topSort(network, w, visited, nodeTopo, index, sol, subscip, varmapfw);
      }

      arc = arc->next_outarc;
   }

   arc = v->inarcs;
   while ( arc != NULL )
   {
      fbv = SCIPhashmapGetImage( varmapfw, arc->negativeFlowBinvar );
      if ( SCIPgetSolVal(subscip, sol, fbv) > 0 )
      {
         w = arc->sourcenode;
         assert( 0 <= w->nodeposition && w->nodeposition < network->numnodes );

         if ( visited[w->nodeposition] < 0 )
            topSort(network, w, visited, nodeTopo, index, sol, subscip, varmapfw);
      }

      arc = arc->next_inarc;
   }

   /* post-visit: store node in topological order */
   nodeTopo[*index] = v;
   (*index)--;
}

/** Delets all constraints which are necessary for the calculation of the node mixing ratio */
static
SCIP_RETCODE deleteMixingCons(
   SCIP*                 subscip,            /**< subscip */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   /* delete constraints for calculation of the mixing ratio:
      TYPE-6-NodeMixingRatioCalculation,
      TYPE-7Mass-linkNodeMixRatioToArcMixRatio and constraints that link node to arc mixing ratio
   */
   SCIP_CONS** conss;
   int nconss;
   int i;

   conss = SCIPgetConss(subscip);
   nconss = SCIPgetNConss(subscip);

   for (i = nconss - 1; i >= 0; --i)
   {
      SCIP_CONS* cons;

      cons = conss[i];
      assert( cons != NULL );

      if ( strcmp("nonlinear", SCIPconshdlrGetName(SCIPconsGetHdlr(cons))) == 0 )
      {
         if ( probdata->nodeMu )
         {
            if ( SCIPstrAtStart(SCIPconsGetName(cons), "TYPE-6-NodeMixingRatioCalculation", 33) )
            {
               SCIP_CALL( SCIPdelCons(subscip, cons) );
            }
         }

         if ( probdata->flowConsMixingNode )
         {
            if ( SCIPstrAtStart( SCIPconsGetName( cons ), "TYPE-6-MixtureFlowCons", 22 ) )
            {
               SCIP_CALL( SCIPdelCons( subscip, cons ) );
            }
         }

         /* if hybrid model */
         /* TODO check if also linear cons!!!!*/
         if ( ! probdata->nodeMu && ! probdata->flowConsMixingNode )
         {
            if ( SCIPstrAtStart(SCIPconsGetName(cons), "TYPE-6-MixtureFlowCons", 22) )
            {
               SCIP_CALL( SCIPdelCons(subscip, cons) );
            }
            else if ( SCIPstrAtStart(SCIPconsGetName(cons), "TYPE-6-NodeMixingRatioCalculation", 33) )
            {
               SCIP_CALL( SCIPdelCons(subscip, cons) );
            }
         }
      }

      if ( strcmp("linear", SCIPconshdlrGetName(SCIPconsGetHdlr(cons))) == 0 )
      {
         if ( SCIPstrAtStart(SCIPconsGetName(cons), "TYPE-3-EXIT-fixNodeMixRatio", 27) )
         {
            SCIP_CALL( SCIPdelCons( subscip, cons) );
         }
         else if ( SCIPstrAtStart(SCIPconsGetName(cons), "TYPE-1A-2Arcs-fixMuMassCon", 26) )
         {
            SCIP_CALL(SCIPdelCons(subscip, cons));
         }
         else if ( SCIPstrAtStart(SCIPconsGetName(cons), "TYPE-1N-2Arcs-fixNodeMu", 23) )
         {
            SCIP_CALL(SCIPdelCons(subscip, cons));
         }
      }
   }

   SCIPdebugMessage("Finished deleting constraints in subscip");

   /* print problem file for debugging */
#ifndef NDEBUG
   SCIP_CALL( SCIPwriteOrigProblem(subscip, "deletedMixingCons.cip", "cip", FALSE) );
#endif
   SCIP_CALL( SCIPwriteOrigProblem(subscip, "deletedMixingCons.cip", "cip", FALSE) );
   return SCIP_OKAY;
}

/** This function gets the solution of a gas flow problem and a topological ordering of the underlying graph.
 *  It calculates the mixing ratio based on the flow from the solution.
 */
static
void updateMixingRatioTopo(
   SCIP_PROBDATA*        probdata,           /**< probdata of the original scip: used to access the nodes of the network */
   GAS_Node**            nodeTopo,           /**< Topological ordering of the nodes */
   SCIP_Real*            mixingRatios,       /**< List of mixing ratios which is updated */
   SCIP_Real*            mixingRatiosCopy,   /**< Copy of the list of mixing ratios which is updated */
   SCIP_Real*            diff,               /**< max diff between old and new mixing ratio: calculated in this function */
   SCIP_HASHMAP*         varmapfw,           /**< Hashmap to get the variables in the subcip which correspond to the original scip variables */
   SCIP*                 subscip,            /**< Subcip instance */
   SCIP_SOL*             sol                 /**< Solution of subscip */
   )
{
   /* copy mixing ratios before update */
   for (int k = 0; k < probdata->network->numnodes; k++)
      mixingRatiosCopy[k] = mixingRatios[k];

   /* loop through nodes and then select the one next in the topological ordering */
   for (int k = 0; k < probdata->network->numnodes; ++k)
   {
      SCIP_VAR* fbv;
      SCIP_VAR* flowvar;
      GAS_Node* node;
      GAS_Arc* arc;
      SCIP_Real mu = -1.0;
      SCIP_Real Q = 0.0;
      SCIP_Real Qmu = 0.0;
      SCIP_Real mr = 1.0;

      /* process nodes in topological ordering */
      node = nodeTopo[k];
      assert( node != NULL );

      /* check if node is an entry and handle accodingly */
      if ( node->type == ENTRY || node->flow > 0.0 )
      {
         SCIP_Real muNode;

         /* calculate the mass% mu */
         muNode = (node->molarMass - probdata->molarMass2) / (probdata->molarMass1 - probdata->molarMass2);

         Q = Q + node->flow;
         Qmu = Qmu + muNode * Q;
      }

      /* Determine the arcs with inflow and their corresponding nodes */
      arc = node->inarcs;
      while ( arc != NULL )
      {
         fbv = SCIPhashmapGetImage(varmapfw, arc->positiveFlowBinvar);
         flowvar = SCIPhashmapGetImage(varmapfw, arc->flowvar);
         assert( fbv != NULL );
         assert( flowvar != NULL );

         if ( SCIPisEQ(subscip, SCIPgetSolVal(subscip, sol, fbv), 1.0) && SCIPgetSolVal(subscip, sol, flowvar) > 0.0)
         {
            Q = Q + SCIPgetSolVal(subscip, sol, flowvar);
            mr = mixingRatios[arc->sourcenode->nodeposition];
            Qmu = Qmu + (SCIPgetSolVal(subscip, sol, flowvar) * mr);
         }
         arc = arc->next_inarc;
      }

      arc = node->outarcs;
      while ( arc != NULL )
      {
         fbv = SCIPhashmapGetImage(varmapfw, arc->negativeFlowBinvar);
         flowvar = SCIPhashmapGetImage(varmapfw, arc->flowvar);
         assert( fbv != NULL );
         assert( flowvar != NULL );

         if ( SCIPisEQ(subscip, SCIPgetSolVal( subscip, sol, fbv ), 1.0) && SCIPgetSolVal(subscip, sol, flowvar) < 0 )
         {
            Q = Q - SCIPgetSolVal(subscip, sol, flowvar);
            mr = mixingRatios[arc->targetnode->nodeposition];
            Qmu = Qmu - (SCIPgetSolVal(subscip, sol, flowvar) * mr);
         }
         arc = arc->next_outarc;
      }

      /* calculate and set node mixing ratio at current node -- only if flow over node */
      if ( Q > 0.0 )
      {
         mu = Qmu / Q;
         assert( 0 <= mu );
         assert( mu <= 1 );

         /* set mixing ratio */
         mixingRatios[node->nodeposition] = mu;

         /* update diff */
         if ( REALABS(mixingRatiosCopy[node->nodeposition] - mixingRatios[node->nodeposition]) > *diff )
         {
#ifdef SCIP_MORE_DEBUG
             SCIPinfoMessage(subscip, NULL, "Node with current max change [%f] in mixing ratio: %s Position in toplogical ordering: %i\n", *diff, node->id, k);
            /*NOTE: it can happen that the same max change in mr happens in the first two consecutive iterations:
             If after the first iteration one node has no flow over it then its mr=1 if then after the second iteration it
             gets flow over it then its mr changes. Example CYCLE.net
             */
#endif
            *diff = ABS( mixingRatiosCopy[node->nodeposition] - mixingRatios[node->nodeposition] );
         }

#ifdef SCIP_MORE_DEBUG
         SCIPinfoMessage(subscip, NULL, "Node: <%s> Mixing ratio <%f>\n", node->id, mu);
#endif
      }
      else
      {
#ifdef SCIP_MORE_DEBUG
         SCIPinfoMessage(subscip, NULL, "Node: <%s> has no flow over it! Set mixing ratio to <%f>\n", node->id, 1.0);
#endif
         /* set mu to something such that no assert flies */
         mu = 1.0;
         /* set mixing ratio */
         mixingRatios[node->nodeposition] = mu;
      }
   }
}


/** This function gets a solution of a subscip and builds a solution of the original scip.
 *  If the subscip solution is feasible it is added otherwise newsol is an infeasible sol in the original scip space   */
static
SCIP_RETCODE transferSolToOrigSCIP(
   SCIP*                 scip,               /**< original scip */
   SCIP*                 subscip,            /**< sub scip */
   SCIP_SOL*             sol,                /**< solution of sub scip */
   SCIP_SOL**            newsol,             /**< solution of original scip which is generated here */
   SCIP_HEUR*            heur,               /**< heur because needed */
   int                   nvars,              /**< number of variables in subscip */
   SCIP_VAR**            subvars,            /**< variables in subscip */
   SCIP_VAR**            vars,               /**< original scip variables */
   SCIP_RESULT*          result,             /**< scip result */
   SCIP_Bool*            transfered          /**< pointer to bool: true if solution is transfered */
   )
{
   SCIP_Bool success = FALSE;

   SCIP_CALL( SCIPcreateOrigSol(scip, newsol, heur) );

   for (int j = 0; j < nvars; ++j)
   {
      SCIP_Real val = SCIPgetSolVal(subscip, sol, subvars[j]);
      //SCIP_Real LB = SCIPvarGetLbGlobal(vars[j]);
      //SCIP_Real UB = SCIPvarGetUbGlobal(vars[j]);
      //assert( ! SCIPisEQ(scip, LB, UB) || SCIPisFeasEQ(scip, val, LB) ); / why necessary?

      SCIP_CALL( SCIPsetSolVal(scip, *newsol, vars[j], val) );
   }

   SCIP_CALL( SCIPcheckSolOrig(scip, *newsol, &success, FALSE, FALSE) );

   if ( success )
   {
      *transfered = TRUE;
      SCIP_CALL( SCIPaddSol(scip, *newsol, &success) );
      SCIPdebugMsg(scip, "Sub-SCIP found a feasible solution of value %g (stored = %u).\n", SCIPgetSolOrigObj(scip, *newsol), success);
      *result = SCIP_FOUNDSOL;
      SCIPinfoMessage(scip, NULL, "Transfered solution of ADM to original Problem \n");
   }
   else
   {
      SCIPdebugMsg(scip, "Sub-SCIP solution infeasible.\n");
   }

   return SCIP_OKAY;
}


static
SCIP_RETCODE changeObjective(
   SCIP* subscip                 /**< sub scip */
)
{

   SCIP_VAR** vars;
   int nvars;

   /* get variables */
   nvars = SCIPgetNOrigVars(subscip);
   vars = SCIPgetOrigVars(subscip);


   for (int i = 0; i < nvars; ++i)
   {
      SCIP_VAR* var = vars[i];

      if ( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      {

         if ( SCIPstrAtStart( SCIPvarGetName( var ), "CS_binvar", 9 ) )
         {
            SCIP_CALL( SCIPchgVarObj(subscip, var, 1.0) );
         }
         
      }
      else if ( SCIPstrAtStart( SCIPvarGetName( var ), "p#", 2 ) )
      {
         SCIP_CALL( SCIPchgVarObj(subscip, var, 0.0) );
      }
   }

   return SCIP_OKAY;
 
}




/** Gets the solution of ADM, creates a new subscip of the original.
 *  Afterwards creates and solves a NLP by fixing the integer variables to the ADM integer values.
*/
static
SCIP_RETCODE solveNLP(
   SCIP*                 scip,               /**< original scip instance */
   SCIP_SOL*             solADM,             /**< subscip solution (of ADM) */
   SCIP_SOL**            NLPsol,             /**< pointer so the NLP solution which is obtained here and then transfered to the orig scip.  */
   SCIP_HEUR*            heur,               /**< heur in order to get name */
   SCIP_RESULT*          result,             /**< pointer to result  */
   SCIP_Bool*            NLPsuccess,         /**< pointer to bool which stores if the NLP solving and transferring of the solution was successfull */
   SCIP_Real*            remaining_time,     /**< time limit for solving the NLP */
   SCIP_Real*            elapsed_time        /**< elapsed time in NLP */
   )
{
   SCIP_HASHMAP* varmapfw;
   char probname[SCIP_MAXSTRLEN];
   SCIP_RETCODE retcode;
   SCIP_Real timelimit;
   SCIP_Bool success;
   SCIP_VAR** subvars;
   SCIP_VAR** vars;
   SCIP* subscip;
   int nvars;
   SCIP_STATUS stat;
   const char* statusstr = NULL;

   assert( scip != NULL);

   /* check whether there is enough time and memory left */
   SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

   if ( ! success )
      return SCIP_OKAY;

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* get variables */
   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(scip), nvars) );

   /* create copy of original problem - we use the original in order to avoid aggregations */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "heur_NLP_ADM");
   SCIP_CALL( SCIPcopyOrig(scip, subscip, varmapfw, NULL, probname, FALSE, FALSE, TRUE, &success) );
   if ( ! success )
   {
      SCIPerrorMessage("Was not able to copy SCIP in solveNLP in ADM.\n");
      return SCIP_ERROR;
   }

   /* copy subproblem variables into the same order as the source SCIP variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   for (int i = 0; i < nvars; ++i)
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* first debugging stuff*/
#ifdef SCIP_MORE_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 10000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* copy time and memory limits */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );

#ifdef NLP_TL
   SCIP_CALL( SCIPgetRealParam(subscip, "limits/time", &timelimit) );
   if ( timelimit > *remaining_time )
   {
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", *remaining_time) );
   }
#endif

   /* stop after reaching a gap of ...*/
#ifdef NLP_GAP
   SCIP_CALL( SCIPsetRealParam( subscip, "limits/gap", 0.05 ) );
#endif

   /* stop after not improving the primal solution for x nodes */
#ifdef NLP_STALLNODE
   SCIP_CALL( SCIPsetLongintParam( subscip, "limits/stallnodes", 5000) );
#endif

   /* stop if primal solution is found */
#ifdef NLP_SOL
   SCIP_CALL( SCIPsetIntParam( subscip, "limits/solutions", 1 ) );
#endif   

   SCIPdebugMessage("Start fixing integer variables in subscip to create NLP\n");

   /* get integer variables from ADM solution and fix in copy of orig scip */
   for (int i = 0; i < nvars; ++i)
   {
      SCIP_VAR* var = vars[i];
      SCIP_VAR* subvar = subvars[i];
      SCIP_Real val;
      SCIP_Bool infeasible;
      SCIP_Bool fixed;

      /* Only fix integer variables */
      if ( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      {
         infeasible = TRUE;
         fixed = FALSE;

         /* Get value from ADM solution */
         val = SCIPgetSolVal(scip, solADM, var);

         /* If binary, round to 0 or 1 */
         if ( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
         {
            if ( SCIPisFeasEQ(scip, val, 1.0) )
               val = 1.0;
            else
               val = 0.0;
         }

         /* Fix variable in subscip */
         SCIP_CALL( SCIPfixVar(subscip, subvar, val, &infeasible, &fixed) );
         assert( ! infeasible );
         assert( fixed );
      }
   }

   /* Now the integer variables in the subscip should be all fixed */
   SCIPdebugMessage("Start solving NLP subscip\n");

#ifndef NDEBUG
   SCIP_CALL( SCIPwriteOrigProblem(subscip, "NLPprobADM.cip", NULL, FALSE) );
#endif

   /* solve subscip */
   retcode = SCIPsolve(subscip);
   *elapsed_time += SCIPgetSolvingTime(subscip); // solving time printed to console with status down below

#ifdef SCIP_MORE_DEBUG
   SCIPinfoMessage(scip, NULL, "NLP solved with code: %d \n", retcode);
#endif

   if ( retcode != SCIP_OKAY )
   {
      SCIPerrorMessage("Error while solving NLP subproblem in ADM; sub-SCIP terminated with code <%d>.\n", retcode);
      SCIP_CALL( SCIPfree(&subscip) );
      return SCIP_ERROR;
   }

   /* get best solution and check that the problem is not infeasible */
   stat = SCIPgetStatus(subscip);

   switch (stat)
   {
   case SCIP_STATUS_OPTIMAL: statusstr = "OPTIMAL"; break;
   case SCIP_STATUS_INFEASIBLE: statusstr = "INFEASIBLE"; break;
   case SCIP_STATUS_STALLNODELIMIT: statusstr = "STALLNODELIMIT"; break;
   case SCIP_STATUS_GAPLIMIT: statusstr = "GAPLIMIT"; break;
   case SCIP_STATUS_TIMELIMIT: statusstr = "TIMELIMIT"; break;
   case SCIP_STATUS_SOLLIMIT: statusstr = "SOLLIMIT"; break;
      // add other cases as needed
   default: statusstr = "UNKNOWN"; break;
   }
   assert( statusstr != NULL );
   SCIPinfoMessage(scip, NULL, "Solved NLP with status %s and gap <%f> in <%f> seconds\n", statusstr, SCIPgetGap(subscip), SCIPgetSolvingTime(subscip));

   if ( stat == SCIP_STATUS_OPTIMAL || stat == SCIP_STATUS_STALLNODELIMIT || stat == SCIP_STATUS_GAPLIMIT || stat == SCIP_STATUS_SOLLIMIT )
   {
      SCIP_Bool transfered = FALSE;
      SCIP_SOL* sol = SCIPgetBestSol(subscip);

      nvars = SCIPgetNOrigVars(scip);

      /* the solution obtained above is a subscip solution. Transfer it to the original scip:*/
      SCIP_CALL( transferSolToOrigSCIP(scip, subscip, sol, NLPsol, heur, nvars, subvars, vars, result, &transfered) );

      /* NLPsol now is the solution of the subcip but in the original scip */
      if ( transfered )
      {
         *NLPsuccess = TRUE;
      }
      else
      {
         SCIPwarningMessage(scip, "BIG PROBLEM");
      }

      /* free everything */
      SCIPhashmapFree(&varmapfw);
      /* free buffer */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
   }
   else if ( stat == SCIP_STATUS_INFEASIBLE )
   {
      /* warning messsage */
      SCIPwarningMessage( scip, "NLP-Problem in ADM is infeasible\n" );  // this is not displayed anymore!!! use infoMessage instead
      SCIPinfoMessage( scip, NULL, "NLP-Problem in ADM is infeasible\n" );

#ifndef NDEBUG
      /* write problem in order to check what happend*/
      SCIP_CALL( SCIPwriteOrigProblem( subscip, "infNLPprobADM.cip", NULL, FALSE ) );
#endif

      /* free everything */
      SCIPhashmapFree(&varmapfw);
      /* free buffer */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
   }
   else
   {
      /* warning messsage */
      SCIPwarningMessage(scip, "Unexpected SCIP status of NLP solution: Abort");
      SCIPinfoMessage( scip, NULL, "Unexpected SCIP status of NLP solution: Abort\n" );
      /* free everything */
      SCIPhashmapFree(&varmapfw);
      /* free buffer */
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
   }

   return SCIP_OKAY;
}


/** Main execution method of the ADM algorithm */
static
SCIP_DECL_HEUREXEC(heurExecADM)
{
   SCIP_PROBDATA* probdata = SCIPgetProbData( scip );
   GAS_Arc* arc;
   GAS_Node* node;
   int k, i;
   char probname[SCIP_MAXSTRLEN];
   SCIP_HEURDATA* heurdata;
   SCIP_HASHMAP* varmapfw;
   SCIP_RETCODE retcode;
   SCIP_Real* Lambdas;
   SCIP_Real* mixingRatios;
   SCIP_Real* mixingRatiosCopy;
   SCIP_Real* prevvals = NULL;
   SCIP_Bool success;
   SCIP_Bool transfered = FALSE;
   SCIP_VAR** subvars;
   SCIP_VAR** vars;
   SCIP_Real* flowBW;
   SCIP_Real* flowBWcopy;
   SCIP* subscip;
   int nvars;
   SCIP_SOL* sol = NULL;
   SCIP_STATUS stat;
   SCIP_Bool termination = FALSE;
   int counter = 0;
   int totalChangedFlowDir = 0;
   int totalChangedFlowDirCopy = 0;
   int changedFlowDir = 0;
   int NLPcounter = 0;
   SCIP_Real diff = 1.0;
   GAS_Node** nodeTopo;
   SCIP_SOL* newsol;
   SCIP_Real elapsed_time = 0.0;
   SCIP_Real global_timelimit = 0.0;
   SCIP_Real remaining_time;
   const char* statusstr = NULL;

   assert( scip != NULL);
   assert( result != NULL );

   /* The algorithm is implemented to work with the mass% model */
   if ( ! (probdata->linearMixing && (probdata->nodeMu || probdata->flowConsMixingNode || probdata->flowConsMixing)) )
     return SCIP_OKAY;

   assert( probdata->linearMixing );


   *result = SCIP_DIDNOTRUN;

   SCIP_CALL( SCIPallocBufferArray(scip, &Lambdas, probdata->network->numarcs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mixingRatios, probdata->network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mixingRatiosCopy, probdata->network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nodeTopo, probdata->network->numnodes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &flowBW, 2 * probdata->network->numarcs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &flowBWcopy, 2 * probdata->network->numarcs) );

   /* initialize mixingRatios */
   for (k = 0; k < probdata->network->numnodes; k++)
   {
      mixingRatios[k] = -1.0;
      mixingRatiosCopy[k] = -1.0;
   }

   /* initialize Lambdas */
   for (k = 0; k < probdata->network->numarcs; k++)
      Lambdas[k] = -1.0;

   /* initialize flowBW */
   for (k = 0; k < 2 * probdata->network->numarcs; k++)
   {
      flowBW[k] = -1.0;
      flowBWcopy[k] = -1.0;
   }

   /* do not run recursively */
   if ( SCIPgetSubscipDepth( scip ) > 0 )
     return SCIP_OKAY;

  /* check whether there is enough time and memory left */
  SCIP_CALL( SCIPcheckCopyLimits(scip, &success) );

  if ( ! success )
     return SCIP_OKAY;

  SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "running heuristic <%s> ...\n", SCIPheurGetName(heur));
  SCIPdebugMsg(scip, "Running <%s> ...\n", SCIPheurGetName(heur));
  *result = SCIP_DIDNOTFIND;

  /* get data */
  heurdata = SCIPheurGetData(heur);
  assert( heurdata != NULL );

  /* Check if heur was already run - if yes do not rerun */
  if ( heurdata->alreadyRun )
  {
     /* free everything */
     SCIPfreeBufferArray(scip, &mixingRatios);
     SCIPfreeBufferArray(scip, &mixingRatiosCopy);
     SCIPfreeBufferArray(scip, &nodeTopo);
     SCIPfreeBufferArray(scip, &Lambdas);
     SCIPfreeBufferArray(scip, &flowBW);
     SCIPfreeBufferArray(scip, &flowBWcopy);

     SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Abort ADM already run <%s> ...\n", SCIPheurGetName(heur));

     return SCIP_OKAY;
  }
  else
     heurdata->alreadyRun = TRUE;

  /* initialize the subproblem */
  SCIP_CALL( SCIPcreate(&subscip) );

  /* get variables */
  nvars = SCIPgetNOrigVars(scip);
  vars = SCIPgetOrigVars(scip);

  SCIP_CALL( SCIPallocBufferArray(scip, &prevvals, nvars) );

  /* create the variable mapping hash map */
  SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(scip), nvars) );

  /* create copy of original problem */
  (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "heur_ADM");
  SCIP_CALL( SCIPcopyOrig(scip, subscip, varmapfw, NULL, probname, FALSE, FALSE, TRUE, &success) );
  if ( ! success )
  {
     SCIPerrorMessage("Was not able to copy SCIP during <%s>.\n", SCIPheurGetName(heur));
     return SCIP_ERROR;
  }

  /* copy subproblem variables into the same order as the source SCIP variables */
  SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
  for (i = 0; i < nvars; ++i)
     subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

  /* Initilization of ADM | only done once at the beginning
     1. Delete Mixing ratio constraints
     2. set/ fix all mixing ratio to 1 (natural gas) except for already fixed values
   */

#ifdef CHG_OBJ  
  SCIP_CALL( changeObjective(subscip) );
#endif
  
  /* first debugging stuff*/
#ifdef SCIP_MORE_DEBUG
   /* for debugging, enable full output */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
   SCIP_CALL( SCIPsetIntParam(subscip, "display/freq", 10000) );
#else
   /* disable statistic timing inside sub SCIP and output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
   SCIP_CALL( SCIPsetBoolParam(subscip, "timing/statistictiming", FALSE) );
#endif

   /* copy time and memory limits */
   SCIP_CALL( SCIPcopyLimits(scip, subscip) );

   // get main SCIP global timelimit
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &global_timelimit) );
   /* Set ADM timelimit to 70% of global TL */
   global_timelimit = 0.7 * global_timelimit;
   //SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", global_timelimit) ); // In the while loop below the time limit for each iteration is additionally bounded by 10min


   /* set strong node limit */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->nodelimit) );

   /* delete constraints for calculation of the mixing ratio: TYPE-6-NodeMixingRatioCalculation,
     TYPE-7Mass-linkNodeMixRatioToArcMixRatio and constraints that link node to arc mixing ratio
   */
   SCIP_CALL( deleteMixingCons(subscip, probdata) );

   /* Initilaize mixing ratio variables in arrary */
   /* loop through all nodes and fix the variables */
   for (k = 0; k < probdata->network->numnodes; ++k)
   {
      int numInarcs = 0;
      int numOutarcs = 0;
      int sum;

      node = &probdata->network->nodes_ptr[k];

      assert( node != NULL );
      assert( node->type != UNKNOWNNODETYPE );

      /* Determine the number of in and outarcs */
      arc = node->inarcs;
      while ( arc != NULL )
      {
         numInarcs++;
         arc = arc->next_inarc;
      }
      arc = node->outarcs;
      while ( arc != NULL )
      {
         numOutarcs++;
         arc = arc->next_outarc;
      }

      sum = numInarcs + numOutarcs;

      /* In the code below all mixing ratio variables except for sources with 1 arc are fixed
         Alternatively one could ignore that we already know the correct mixing ratio for those
         and set those to 1.
      */

      /* mu is already fixed to the correct value for entries with only one arc*/
      if ( node->type == ENTRY && sum == 1 )
      {
         assert( SCIPisEQ(scip, SCIPvarGetLbLocal(node->nodeMixingRatio), SCIPvarGetUbLocal(node->nodeMixingRatio)) );
         mixingRatios[node->nodeposition] = SCIPvarGetLbLocal(node->nodeMixingRatio);
      }
      else
      {
         mixingRatios[node->nodeposition] = 1.0;
      }
   }

   /* check if  all mixing ratios have been set */
   for (k = 0; k < probdata->network->numnodes; k++)
   {
      assert( mixingRatios[k] >= 0.0 );
   }

   SCIPinfoMessage(scip, NULL, "Initialization finished\n");

   /* INITILIZATION FINISHED */

   /* Next Steps:
      1. solve subscip
      2. get topological ordering, calculate and fix mixing ratio
      3. repeat
   */

   /* Termination condition of ADM*/
   while ( ! termination )
   {
      SCIP_CALL( SCIPfreeTransform(subscip) );

      /* Loop through nodes and fix mixing ratio to previous calculated values in the mixingRatios array */
      for (k = 0; k < probdata->network->numnodes; ++k)
      {
         SCIP_Bool inf = FALSE;
         SCIP_Bool fixed = FALSE;
         SCIP_VAR* var;

         node = &probdata->network->nodes_ptr[k];
         assert( node != NULL );
         assert( node->type != UNKNOWNNODETYPE );
         var = SCIPhashmapGetImage(varmapfw, node->nodeMixingRatio );
         assert( var != NULL );
         SCIP_CALL( SCIPfixVar(subscip, var, mixingRatios[node->nodeposition], &inf, &fixed) );
         assert( ! inf );
         assert( fixed );
      }

#ifndef NDEBUG
      SCIP_CALL( SCIPwriteOrigProblem(subscip, "OrigProbADM.cip", NULL, FALSE) );
#endif

      SCIP_CALL( SCIPpresolve(subscip) );

      /* scip stops solving as soon as a primal solution is found*/
#ifdef STOP_IF_SUB_SCIP_FOUND_SOL
      SCIP_CALL( SCIPsetIntParam( subscip, "limits/solutions", 1 ) );
#endif

      /* stop after reaching a gap of ...*/
#ifdef STOP_IF_SUB_SCIP_SMALL_GAP
      SCIP_CALL( SCIPsetRealParam( subscip, "limits/gap", 0.1 ) );
#endif

      /* stop after not improving the primal solution for x nodes */
#ifdef STALLNODE_LIMIT
      SCIP_CALL( SCIPsetLongintParam( subscip, "limits/stallnodes", 1000) );
#endif

      SCIP_CALL( SCIPsetRealParam( subscip, "numerics/dualfeastol", 0.001 ) ); // seems to help

      //SCIP_CALL( SCIPsetRealParam(subscip, "numerics/feastol", 0.001) );
      //SCIP_CALL( SCIPsetParam( subscip, "lp/initalgorithm", "c" ) );  // seems to help
      //SCIP_CALL( SCIPsetParam( subscip, "lp/resolvealgorithm", "c" ) );

      // compute remaining time
      remaining_time = global_timelimit - elapsed_time;
      if ( remaining_time <= 0.0 )
      {
         SCIPinfoMessage( scip, NULL, "No time remaining for ADM: Abort\n" );
         /* free everything */
         SCIPfreeBufferArray(scip, &mixingRatios);
         SCIPfreeBufferArray(scip, &mixingRatiosCopy);
         SCIPfreeBufferArray(scip, &nodeTopo);
         SCIPfreeBufferArray(scip, &Lambdas);
         SCIPfreeBufferArray(scip, &flowBW);
         SCIPfreeBufferArray(scip, &flowBWcopy);
         return SCIP_OKAY;
      }

      SCIPinfoMessage(scip, NULL, "Remaining time for ADM: %f \n", remaining_time);
      /* each iteration of at most 10 min */
      if ( remaining_time > 600.0 )
      {
         SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", 600.0) );
      }
      else
         SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", remaining_time) );


      /* solve subscip */
      SCIPinfoMessage(scip, NULL, "Iteration %i: Start solving subscip\n", counter);
      retcode = SCIPsolve(subscip);

      // elapsed time of sub SCIP
      elapsed_time += SCIPgetSolvingTime(subscip); // Is printed to console below the solving status

      if ( retcode != SCIP_OKAY )
      {
         SCIPerrorMessage("Error while solving subproblem in <%s>; sub-SCIP terminated with code <%d>.\n", SCIPheurGetName(heur), retcode);
         SCIP_CALL( SCIPfree(&subscip) );
         return SCIP_ERROR;
      }

      /* get best solution and check that the problem is not infeasible */
      stat = SCIPgetStatus(subscip);

      switch ( stat )
      {
      case SCIP_STATUS_OPTIMAL: statusstr = "OPTIMAL"; break;
      case SCIP_STATUS_INFEASIBLE: statusstr = "INFEASIBLE"; break;
      case SCIP_STATUS_STALLNODELIMIT: statusstr = "STALLNODELIMIT"; break;
      case SCIP_STATUS_GAPLIMIT: statusstr = "GAPLIMIT"; break;
      case SCIP_STATUS_TIMELIMIT: statusstr = "TIMELIMIT"; break;
      case SCIP_STATUS_SOLLIMIT: statusstr = "SOLLIMIT"; break;   
         // add other cases as needed
      default: statusstr = "UNKNOWN"; break;
      }
      assert( statusstr != NULL );
      SCIPinfoMessage(scip, NULL, "Solved subscip with status %s and gap <%f> in <%f> seconds\n", statusstr, SCIPgetGap(subscip), SCIPgetSolvingTime(subscip));
      SCIPinfoMessage(scip, NULL, "Total elapsed time in ADM: %f \n", elapsed_time);

      if ( stat != SCIP_STATUS_INFEASIBLE && SCIPgetNSols(subscip) != 0)
      {
         sol = SCIPgetBestSol(subscip);
      }
      else if ( stat == SCIP_STATUS_INFEASIBLE || SCIPgetNSols(subscip) == 0 )
      {
#ifndef NDEBUG
         /* print infeasible problem for debugging */
         SCIP_CALL( SCIPwriteOrigProblem(subscip, "infProbADM.cip", NULL, FALSE) );
         SCIP_CALL( SCIPwriteTransProblem(subscip, "infTransProbADM.cip", NULL, FALSE) );
#endif

         /* free everything */
         SCIPfreeBufferArray(scip, &mixingRatios);
         SCIPfreeBufferArray(scip, &mixingRatiosCopy);
         SCIPfreeBufferArray(scip, &nodeTopo);
         SCIPfreeBufferArray(scip, &Lambdas);
         SCIPfreeBufferArray(scip, &flowBW);
         SCIPfreeBufferArray(scip, &flowBWcopy);
         SCIPfreeBufferArray(scip, &prevvals);
         /* free hash map */
         SCIPhashmapFree(&varmapfw);
         /* free buffer */
         SCIPfreeBufferArray(scip, &subvars);
         SCIP_CALL(SCIPfree(&subscip));

         /* warning messsage */
         SCIPwarningMessage(scip, "Problem is either infeasible or timelimit reached: Abort\n");

         *result = SCIP_DIDNOTFIND;

         return SCIP_OKAY;
      }
      else
      {
         /* warning messsage */
         SCIPwarningMessage(scip, "Unexpected SCIP status for ADM solution: Abort");

         *result = SCIP_DIDNOTFIND;

         return SCIP_OKAY;
      }

      /* Calculate topological ordering based on the flow directions of the current solution */
      {
         int n = probdata->network->numnodes;
         int depth[n];
         int index = n - 1;

         /* initilization: all nodes univisited */
         for (i = 0; i < n; ++i)
            depth[i] = -1;

         /* call topSort for each node */
         for (i = 0; i < n; ++i)
         {
            node = &probdata->network->nodes_ptr[i];

            if (depth[i] < 0) // (depth[i] < 0 && node->type == ENTRY)
               topSort(probdata->network, node, depth, nodeTopo, &index, sol, subscip, varmapfw);
         }
      }

      /* check if flow directions changed */
      {
         /* copy values for comparison */
         for (k = 0; k < probdata->network->numarcs; k++)
            flowBWcopy[k] = flowBW[k];

         /* update values to current solution */
         for (k = 0; k < probdata->network->numarcs; k++)
         {
            arc = &probdata->network->arcs_ptr[k];
            flowBW[k] = SCIPgetSolVal(scip, sol, arc->positiveFlowBinvar);
         }

         totalChangedFlowDirCopy = totalChangedFlowDir;

         if ( counter > 0 )
         {
            for (k = 0; k < probdata->network->numarcs; k++)
            {
               if ( ! SCIPisFeasEQ(scip, flowBW[k], flowBWcopy[k]) )
               {
                  ++totalChangedFlowDir;
#ifdef PRINT_ARC_IDS_WITH_CHANGED_FLOW_DIRECTION
                  arc = &probdata->network->arcs_ptr[k];
                  SCIPinfoMessage( scip, NULL, "Flow direction changed on arc %s\n", arc->id );
#endif
               }
            }

            changedFlowDir = totalChangedFlowDir - totalChangedFlowDirCopy;
#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage( scip, NULL, "#Flow direction changed in iteration %i: %i\n", counter, changedFlowDir );
#endif
         }
      }

      /* calculate max difference between old and new mixing ratios */
      diff = 0.0;

      /* Calculate and update mixing ratios in topological order */
      updateMixingRatioTopo(probdata, nodeTopo, mixingRatios, mixingRatiosCopy, &diff, varmapfw, subscip, sol);

      SCIPinfoMessage(scip, NULL, "Iteration %i completed\n", counter);
      SCIPinfoMessage(scip, NULL, "Maximum difference between old and new mixing ratio: %f\n\n", diff);
      counter++;

      /* Write each solution of every iteration to an lsf file - for visualization */
#ifdef VISUALIZE
      if ( counter > 0 )
      {
         const char* name = "ADM";
         const char* nameOrig = probdata->netname;
         const char* scnOrig = probdata->scnname;
         char iter[SCIP_MAXSTRLEN];

         SCIP_CALL( transferSolToOrigSCIP( scip, subscip, sol, &newsol, heur, nvars, subvars, vars, result, &transfered ) );

         /* Write solution of ADM to LSF file with name ADM-iteration */
         snprintf(iter, SCIP_MAXSTRLEN, "%d", counter);
         strcpy(probdata->netname, name);
         strcpy(probdata->scnname, iter); //Iteration of ADM
         SCIP_CALL( writeLSF(scip, newsol, "bin/gasa01") );
         SCIP_CALL( SCIPfreeSol( scip, &newsol) );
      }
#endif

#ifdef NORMAL_TERMINATION
      if ( diff < 1e-8 || counter > 1000 )
      {
         termination = TRUE;
      }

      /* try to transfer solution to original scip and write solution file */
      if ( termination )
      {
         FILE* file = fopen("ADMsol.sol", "w");
         assert( file != NULL );
         SCIP_CALL( SCIPprintBestSol(subscip, file, FALSE) );
         fclose(file);

         SCIP_CALL( transferSolToOrigSCIP(scip, subscip, sol, &newsol, heur, nvars, subvars, vars, result, &transfered) );

         /* Write solution of ADM to LSF file with name ADM-iteration */
         snprintf(iter, SCIP_MAXSTRLEN, "%d", counter);
         strcpy(probdata->netname, name);
         strcpy(probdata->scnname, iter); //Iteration of ADM
         SCIP_CALL( writeLSF(scip, newsol, "bin/gasa01") );

         /* reset of netname and scnname */
         strcpy(probdata->netname, nameOrig);
         strcpy(probdata->scnname, scnOrig);

         SCIP_CALL( SCIPfreeSol(scip, &newsol) );
      }
#endif

      /* different termination criteria: Stop as soon as a primal feasible solution is found */
#ifdef STOP_IF_PRIMAL_SOL_FOUND
      /* try to check if sol is primal feas after second iteration - always infeasible in the first one */
      if ( counter > 1 )
      {
         SCIP_CALL( transferSolToOrigSCIP(scip, subscip, sol, &newsol, heur, nvars, subvars, vars, result, &transfered) );

         if ( ! transfered )
            SCIP_CALL( SCIPfreeSol(scip, &newsol) );
      }

      /* Stop if the solution was successfully transfered */
      if ( transfered )
      {
         FILE* file = fopen("ADMsol.sol", "w");
         assert(file != NULL);
         SCIP_CALL(SCIPprintBestSol(subscip, file, FALSE));
         fclose(file);

         /* Write solution of ADM to LSF file with name ADM-iteration */
         snprintf(iter, SCIP_MAXSTRLEN, "%d", counter);
         strcpy(probdata->netname, name);
         strcpy(probdata->scnname, iter); //Iteration of ADM
         SCIP_CALL( writeLSF(scip, newsol, "bin/gasa01") );
         SCIP_CALL( SCIPfreeSol(scip, &newsol) );

         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "ADM found primal solution. Stoping ADM now.\n");
         /* set termination */
         termination = TRUE;
      }
#endif

      /* If already close to convergence: Take binary variable values from ADM sol
         and try to find a feasible solution for the remaining NLP
         The last criterion: counter > 1 condition is for GL582 since this is successful for many scenarios.
         Should be disabled / removed to simply test the algorithm without looking at performance
      */
      if ( (diff < 0.1 && ! transfered) || (changedFlowDir == 0 && counter > 1 && ! transfered) || (counter > 1 && ! transfered) )
      {
         SCIP_SOL* NLPsol = NULL;
         SCIP_Bool NLPsuccess = FALSE;
         SCIP_Real elapsed_time_NLP = 0.0;
         FILE* file;

         /* increae NLP counter */
         NLPcounter++;

         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Convergence of ADM whithin tolerance [diff = %f] [#chgFlowDir = %i] : Solving NLP now.\n",
            diff, changedFlowDir);

         file = fopen("ADMsol.sol", "w");
         assert(file != NULL);
         SCIP_CALL( SCIPprintBestSol(subscip, file, FALSE) );
         fclose(file);

         /* transfer ADM solution (sol) to original scip solution (newsol) in order to use it in solveNLP */
         SCIP_CALL( transferSolToOrigSCIP(scip, subscip, sol, &newsol, heur, nvars, subvars, vars, result, &transfered) );
         SCIP_CALL( solveNLP(scip, newsol, &NLPsol, heur, result, &NLPsuccess, &remaining_time, &elapsed_time_NLP ) );
         // add elapsed time of NLP to total elapsed time of ADM
         elapsed_time += elapsed_time_NLP;

         SCIP_CALL( SCIPfreeSol(scip, &newsol) );

         if ( NLPsuccess )
         {
            termination = TRUE;
            SCIP_CALL( SCIPfreeSol(scip, &NLPsol) );

            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "NLP found primal feasible solution .\n");
         }
         else if ( NLPcounter > 10 && diff < 1e-5 )
         {
            /* this condition might give a false positive but it ts usefull to ensure not
               spending unneccessary time if the algorithm does not find a solution.
            */
            termination = TRUE;
            SCIPverbMessage( scip, SCIP_VERBLEVEL_NORMAL, NULL, "Problem is likely to be infeasible.\n" );

            /* write the solution to a file in order to compare with the solution in case the problem was not infeasible */
            file = fopen("ADMsolDiff0ButNotFeasible.sol", "w");
            assert( file != NULL );
            SCIP_CALL( SCIPprintBestSol(subscip, file, FALSE) );
            fclose(file);
         }
      }

      /* stop if iteration or time limit is hit */
      if ( counter > 50 || ( elapsed_time / global_timelimit ) > 0.7 )
      {
         termination = TRUE;
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Hit iteration or time limit without finding a solution .\n");
      }

      if ( termination )
      {
         SCIPinfoMessage( scip, NULL, "Total number of changed flow directions: %i\n", totalChangedFlowDir );
      }
   }

   /* free everything */
   SCIPfreeBufferArray(scip, &Lambdas);
   SCIPfreeBufferArray(scip, &mixingRatios);
   SCIPfreeBufferArray(scip, &mixingRatiosCopy);
   SCIPfreeBufferArray(scip, &nodeTopo);
   SCIPfreeBufferArray(scip, &flowBW);
   SCIPfreeBufferArray(scip, &flowBWcopy);
   SCIPfreeBufferArray(scip, &prevvals);
   SCIPhashmapFree(&varmapfw);
   SCIPfreeBufferArray(scip, &subvars);
   SCIP_CALL( SCIPfree(&subscip) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Finished solving <%s>.\n", SCIPheurGetName(heur));

   return SCIP_OKAY;
}



/*
 * primal heuristic specific interface methods
 */

/** creates the ADM primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurADM(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create Gas primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecADM, heurdata) );

   assert( heur != NULL );

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeADM) );

    /* add gas primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodelimit",
      "maximal number of nodes to use in sub-SCIP",
      &heurdata->nodelimit, TRUE, DEFAULT_NODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam( scip, "heuristics/" HEUR_NAME "/transfersols",
      "Transfer solutions of main SCIP to sub-SCIP?",
      &heurdata->transfersols, TRUE, DEFAULT_TRANSFERSOLS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam( scip, "heuristics/" HEUR_NAME "/alreadyrun",
      "Did the heur already run?",
      &heurdata->alreadyRun, FALSE, DEFAULT_ALREADYRUN, NULL, NULL) );

   return SCIP_OKAY;
}
