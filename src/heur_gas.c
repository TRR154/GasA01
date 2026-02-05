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

/**@file   heur_gas.c
 * @ingroup PRIMALHEURISTICS
 * @brief  LNS heuristic that is adapted to gas transportation problems
 * @author Marc Pfetsch
 *
 * We copy the original problem and remove all constraints that contain the positive/negative direction variables
 * (presolving will then eliminate these variables). Then we run the sub-SCIP with a node limit. The hope is that the
 * heuristics of the sub-SCIP will be more successful. In particular, the sub-NLP heuristic is likely to have a larger
 * success rate.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "heur_gas.h"
#include <scip/scip.h>
#include <scip/cons_linear.h>
#include <scip/pub_misc.h>
#include <scip/pub_heur.h>

#include <elements_gas.h>

/* default values for standard parameters that every primal heuristic has in SCIP */
#define HEUR_NAME             "gas"
#define HEUR_DESC             "LNS heuristic adapated to gas transport"
#define HEUR_DISPCHAR         '~'
#define HEUR_PRIORITY         -1100000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERLPNODE
#define HEUR_USESSUBSCIP      TRUE      /**< does the heuristic use a secondary SCIP instance? */

#define DEFAULT_NODELIMIT     3LL       /**< maximal number of nodes to use in sub-SCIP */
#define DEFAULT_TRANSFERSOLS  TRUE      /**< Transfer solutions of main SCIP to sub-SCIP? */

/*
 * Data structures
 */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          nodelimit;          /**< maximal number of nodes to use in sub-SCIP */
   SCIP_Bool             transfersols;       /**< transfer solutions of main SCIP to sub-SCIP */
};


/*
 * Callback methods of primal heuristic
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyGas)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurGas(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeGas)
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

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecGas)
{  /*lint --e{715}*/
   char probname[SCIP_MAXSTRLEN];
   SCIP_HEURDATA* heurdata;
   SCIP_HASHMAP* varmapfw;
   SCIP_RETCODE retcode;
   SCIP_Real timelimit;
   SCIP_Bool success;
   SCIP_CONS** conss;
   SCIP_VAR** subvars;
   SCIP_VAR** vars;
   SCIP* subscip;
   int nconss;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* do not run recursively */
   if ( SCIPgetSubscipDepth(scip) > 0 )
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

   /* initialize the subproblem */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* get variables */
   nvars = SCIPgetNOrigVars(scip);
   vars = SCIPgetOrigVars(scip);

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );

   /* create copy of original problem - we use the original in order to avoid aggregations which would confuse our deletion below */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "heur_gas");
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

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* transfer bounds */
   for (i = 0; i < nvars; ++i)
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(vars[i]);
      ub = SCIPvarGetUbLocal(vars[i]);

      if ( lb > SCIPvarGetLbLocal(subvars[i]) )
      {
         SCIP_CALL( SCIPchgVarLb(subscip, subvars[i], lb) );
      }
      if ( ub < SCIPvarGetUbLocal(subvars[i]) )
      {
         SCIP_CALL( SCIPchgVarUb(subscip, subvars[i], ub) );
      }
   }

   /* transfer solutions */
   if ( heurdata->transfersols )
   {
      SCIP_SOL** sols;
      SCIP_SOL* sol;
      int nsols;
      int j;

      nsols = SCIPgetNSols(scip);
      sols = SCIPgetSols(scip);
      for (i = 0; i < nsols; ++i)
      {
         SCIP_SOL* subsol;
         SCIP_CALL( SCIPcreateSol(subscip, &subsol, NULL) );
         sol = sols[i];

         for (j = 0; j < nvars; ++j)
         {
            SCIP_Real val;

            val = SCIPgetSolVal(scip, sol, vars[j]);
            SCIP_CALL( SCIPsetSolVal(subscip, subsol, subvars[j], val) );
         }
         SCIP_CALL( SCIPaddSolFree(subscip, &subsol, &success) );
      }
   }

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

   SCIP_CALL( SCIPgetRealParam(subscip, "limits/time", &timelimit) );
   if ( timelimit > 600.0 )
   {
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", 600.0) );
   }

   /* set strong node limit */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", heurdata->nodelimit) );

   /* delete constraints involving binary flow variables in order to improve chances of NLP */
   conss = SCIPgetConss(subscip);
   nconss = SCIPgetNConss(subscip);

   for (i = nconss-1; i >= 0; --i)
   {
      SCIP_CONS* cons;

      cons = conss[i];
      assert( cons != NULL );

      if ( strcmp("linear", SCIPconshdlrGetName(SCIPconsGetHdlr(cons)) ) == 0 )
      {
         SCIP_VAR** linvars;
         int nlinvars;
         int l;

         linvars = SCIPgetVarsLinear(scip, cons);
         nlinvars = SCIPgetNVarsLinear(scip, cons);
         for (l = 0; l < nlinvars; ++l)
         {
            if ( SCIPstrAtStart(SCIPvarGetName(linvars[l]), "posFlowBinvar#", 14) )
               break;
            if ( SCIPstrAtStart(SCIPvarGetName(linvars[l]), "negFlowBinvar#", 14) )
               break;
         }

         if ( l < nlinvars )
         {
            SCIP_CALL( SCIPdelCons(subscip, cons) );
         }
      }
      else if ( strcmp("bounddisjunction", SCIPconshdlrGetName(SCIPconsGetHdlr(cons)) ) == 0 )
      {
         /* delete all bound disjunction constraints, since it might couple binary variables */
         SCIP_CALL( SCIPdelCons(subscip, cons) );
      }
   }

#ifdef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPwriteOrigProblem(subscip, "heurgas.cip", "cip", FALSE) );

   SCIP_CALL( SCIPpresolve(subscip) );
   SCIP_CALL( SCIPwriteTransProblem(subscip, "heurgastrans.cip", "cip", FALSE) );
#endif

   /* solve subscip */
   retcode = SCIPsolve(subscip);

   if ( retcode != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Error while solving subproblem in <%s>; sub-SCIP terminated with code <%d>.\n", SCIPheurGetName(heur), retcode);
      SCIPABORT();
   }
   else
   {
      SCIP_SOL** subsols;
      int nsubsols;

#ifdef SCIP_MORE_DEBUG
      /* print solving statistics of subproblem */
      SCIP_CALL( SCIPprintStatistics(subscip, NULL) );
#endif

      /* check, whether a solution was found */
      nsubsols = SCIPgetNSols(subscip);
      if ( nsubsols > 0 )
      {
         GAS_Network* network;
         SCIP_PROBDATA* probdata;
         SCIP_SOL* newsol;
         int j;

         subsols = SCIPgetSols(subscip);
         success = FALSE;

         probdata = SCIPgetProbData(scip);
         assert( probdata != NULL );
         network = probdata->network;

         for (i = 0; i < nsubsols && ! success; ++i)
         {
            SCIP_CALL( SCIPcreateOrigSol(scip, &newsol, heur) );
            for (j = 0; j < nvars; ++j)
            {
               /* avoid binary flow variables */
               if ( ! SCIPstrAtStart(SCIPvarGetName(vars[j]), "posFlowBinvar#", 14) )
               {
                  if ( ! SCIPstrAtStart(SCIPvarGetName(vars[j]), "negFlowBinvar#", 14) )
                  {
                     SCIP_Real val;

                     val = SCIPgetSolVal(subscip, subsols[i], subvars[j]);

                     assert( ! SCIPisEQ(scip, SCIPvarGetLbGlobal(vars[j]), SCIPvarGetUbGlobal(vars[j])) || SCIPisFeasEQ(scip, val, SCIPvarGetLbGlobal(vars[j])) );
                     SCIP_CALL( SCIPsetSolVal(scip, newsol, vars[j], val) );
                  }
               }
            }

            /* try to set flow direction values */
            for (j = 0; j < network->numarcs; ++j)
            {
               GAS_Arc* arc;
               SCIP_Real flowval;

               arc = &(network->arcs_ptr[j]);
               if ( arc->flowvar != NULL )
               {
                  flowval = SCIPgetSolVal(scip, newsol, arc->flowvar);

                  /* try to decide on binary variables in compressor stations */
                  if ( arc->type == CS )
                  {
                     GAS_CS* cs = (GAS_CS*) arc->detailed_info;

                     assert( SCIPisGE(scip, flowval, 0.0) );

                     if ( SCIPisFeasPositive(scip, flowval) )
                     {
                        SCIP_CALL( SCIPsetSolVal(scip, newsol, cs->compressor_binvar, 1.0) );
                        if ( cs->bypass != NULL )
                        {
                           GAS_Valve* valve = (GAS_Valve*) cs->bypass->detailed_info;
                           assert( cs->bypass->type == VALVE );

                           SCIP_CALL( SCIPsetSolVal(scip, newsol, valve->valve_binvar, 0.0) );
                        }
                     }
                     else
                     {
                        SCIP_CALL( SCIPsetSolVal(scip, newsol, cs->compressor_binvar, 0.0) );
                        if ( cs->bypass != NULL )
                        {
                           GAS_Valve* valve = (GAS_Valve*) cs->bypass->detailed_info;
                           assert( cs->bypass->type == VALVE );

                           SCIP_CALL( SCIPsetSolVal(scip, newsol, valve->valve_binvar, 1.0) );
                        }
                     }
                  }
                  else if ( arc->type == CONTROLVALVE )
                  {
                     GAS_Controlvalve* cv;

                     cv = (GAS_Controlvalve*) arc->detailed_info;

                     if ( SCIPisFeasZero(scip, flowval) )
                     {
                        SCIP_CALL( SCIPsetSolVal(scip, newsol, cv->binvar, 0.0) );
                     }
                     else
                     {
                        SCIP_CALL( SCIPsetSolVal(scip, newsol, cv->binvar, 1.0) );
                     }
                  }
                  else if ( arc->type == VALVE )
                  {
                     GAS_Valve* valve = (GAS_Valve*) arc->detailed_info;

                     if ( SCIPisFeasZero(scip, flowval) )
                     {
                        SCIP_CALL( SCIPsetSolVal(scip, newsol, valve->valve_binvar, 0.0) );
                     }
                     else
                     {
                        SCIP_CALL( SCIPsetSolVal(scip, newsol, valve->valve_binvar, 1.0) );
                     }
                  }
                  else
                  {
                     assert( arc->type == PIPE || arc->type == RESISTOR || arc->type == SHORTPIPE );
                     if ( arc->positiveFlowBinvar != NULL )
                     {
                        assert( arc->negativeFlowBinvar != NULL );

                        /* determine direction variables */
                        if ( SCIPisFeasPositive(scip, flowval) )
                        {
                           SCIP_CALL( SCIPsetSolVal(scip, newsol, arc->positiveFlowBinvar, 1.0) );
                           SCIP_CALL( SCIPsetSolVal(scip, newsol, arc->negativeFlowBinvar, 0.0) );
                        }
                        else if ( SCIPisFeasNegative(scip, flowval) )
                        {
                           SCIP_CALL( SCIPsetSolVal(scip, newsol, arc->negativeFlowBinvar, 1.0) );
                           SCIP_CALL( SCIPsetSolVal(scip, newsol, arc->positiveFlowBinvar, 0.0) );
                        }
                     }
                  }
               }
            }

#ifdef SCIP_DEBUG
            SCIP_CALL( SCIPcheckSolOrig(scip, newsol, &success, TRUE, FALSE) );
#else
            SCIP_CALL( SCIPcheckSolOrig(scip, newsol, &success, FALSE, FALSE) );
#endif
            if ( success )
            {
               SCIP_CALL( SCIPaddSolFree(scip, &newsol, &success) );
               SCIPdebugMsg(scip, "Sub-SCIP found a feasible solution of value %g (stored = %u).\n", SCIPgetSolOrigObj(scip, newsol), success);
               *result = SCIP_FOUNDSOL;
            }
            else
            {
               SCIP_CALL( SCIPfreeSol(scip, &newsol) );
               SCIPdebugMsg(scip, "Sub-SCIP solution infeasible.\n");
            }
         }
      }
      else
         SCIPdebugMsg(scip, "Sub-SCIP did not find feasible solution.\n");
   }

   /* free buffer */
   SCIPfreeBufferArray(scip, &subvars);

   /* finally free subscipdata */
   SCIP_CALL( SCIPfree(&subscip) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Finished solving heuristic <%s>.\n", SCIPheurGetName(heur));

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the gas primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurGas(
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
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecGas, heurdata) );

   assert( heur != NULL );

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyGas) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeGas) );

   /* add gas primal heuristic parameters */
   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/nodelimit",
         "maximal number of nodes to use in sub-SCIP",
         &heurdata->nodelimit, TRUE, DEFAULT_NODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "heuristics/" HEUR_NAME "/transfersols",
         "Transfer solutions of main SCIP to sub-SCIP?",
         &heurdata->transfersols, TRUE, DEFAULT_TRANSFERSOLS, NULL, NULL) );

   return SCIP_OKAY;
}
