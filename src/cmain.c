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

/**@file   cmain.c
 * @brief  main file
 * @author Marc Pfetsch
 * @author Oliver Habeck
 */

#define FLOWBOUNDS_OUTPUT

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include <scip/scip.h>
#include <scip/scipdefplugins.h>
#include "elements_gas.h"
#include "probdata_gas.h"
#include "generateModel_gas.h"
#include "cons_pipe_ode.h"
#include "LSFiles_gas.h"
#include "disp_gradient.h"
#include "heur_gas.h"
#include "heur_ADM.h"
#include "presol_flowobbt.h"
#include "heur_lpround.h"

/** define macro to print error message and exit */
#define SCIP_CALL_ERROR(x)   do                                                                               \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPprintError(_restat_);                                                        \
                             return -1;                                                                       \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

/** read comand line arguments */
static
SCIP_Bool readArguments(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_PROBDATA*        probdata,           /**< problem data structure */
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char**          netfilename,        /**< file name for network from arguments */
   const char**          scnfilename,        /**< file name for scenario from arguments */
   const char**          csfilename,         /**< file name for cs from arguments */
   const char**          lsfilename,         /**< ls-file name */
   SCIP_Bool*            writelsf,           /**< whether an LSF-file should be written */
   const char**          settingsname,       /**< name of settings file */
   SCIP_Real*            timelimit,          /**< time limit read from arguments */
   SCIP_Real*            memlimit,           /**< memory limit read from arguments */
   SCIP_Longint*         nodelimit,          /**< node limit read from arguments */
   int*                  dispfreq            /**< display frequency */
   )
{
   int i;
   char usage[5000];

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( argv != NULL );
   assert( netfilename != NULL );
   assert( scnfilename != NULL);
   assert( csfilename != NULL );
   assert( lsfilename != NULL );
   assert( writelsf != NULL );
   assert( settingsname != NULL );
   assert( timelimit != NULL );
   assert( memlimit != NULL );
   assert( nodelimit != NULL );
   assert( dispfreq != NULL );

   /* init usage text */
   (void) snprintf(usage, (unsigned long)5000, "usage: %s <network file> <scenario file> [-a] [-b <cs file>] [-cD <cs file>] [-s <setting file>]", argv[0]);
   strcat(usage, "[-lsf <ls file>] [-writelsf] [-t <time limit>] [-mem <mem limit>] [-m] [-ml] [-nodeMu] [-fcm] [-n <node limit>] [-d <display frequency>]");
   strcat(usage, " [-h] [-minCompSum] [-minCompInc] [-minPressSum] [-noObjective] [-powerLoss] [-allCycles]");
   strcat(usage, " [-relaxUpperBounds <bar>] [-relaxLowerBounds <bar>] [-relaxCVbounds <bar>] [-minSlack] [-maxFlow]\n\n");
   strcat(usage, "  -a                 use the algebraic form\n");
   strcat(usage, "  -m                 use the molar speed of sound mixing model\n");
   strcat(usage, "  -ml                use the mass percatage convex comb. speed of sound mixing model\n");
   strcat(usage, "  -nodeMu            mixing ratio only at nodes\n");
   strcat(usage, "  -fcm               use the mixing model with flow conservation for the mixture\n");
   strcat(usage, "  -fcmn              use the node mixing model with flow conservation for the mixture\n");
   strcat(usage, "  -b <cs file>       use the boxConstraintModel for compressor stations\n");
   strcat(usage, "  -cD <cs file>      use detailed model for compressor stations\n");
   strcat(usage, "  -lsf <ls file>     check if solution specified in the file is feasible\n");
   strcat(usage, "  -allCycles         forbid cyclic flow for all cycles in the graph\n");
   strcat(usage, "  -relax*Bounds      relax lower or upper pressure bounds up to <bar>\n");
   strcat(usage, "  -relaxCVbounds     relax pressure bounds of controlvalves up to <bar>\n");
   strcat(usage, "\nOptions for the objective function:\n");
   strcat(usage, "  -minCompSum        minimize number of used compressors\n");
   strcat(usage, "  -minCompInc        minimize number the pressure increase by compressors\n");
   strcat(usage, "  -minPressSum       minimize the sum of pressures at exits\n");
   strcat(usage, "  -noObjective       do not use an objective\n");
   strcat(usage, "  -powerLoss         minimize the power loss\n");
   strcat(usage, "  -minSlack          minimize difference to original pressure bounds if they are relaxed\n");
   strcat(usage, "  -minSlackPerBound  minimize the sum of slack variables for relaxed pressure bounds\n");
   strcat(usage, "  -maxFlow           maximize the flow\n");
   strcat(usage, "Otherwise maximize the sum of all pressures.\n");

   assert( strlen(usage) < 5000 );

   /* init arguments */
   *timelimit    = 1e20;
   *memlimit     = 1e20;
   *nodelimit    = SCIP_LONGINT_MAX;
   *netfilename  = NULL;
   *scnfilename  = NULL;
   *csfilename   = NULL;
   *lsfilename   = NULL;
   *writelsf     = FALSE;
   *settingsname = NULL;
   *dispfreq     = -1;

   if ( argc <= 2 )
   {
      SCIPinfoMessage(scip, NULL, "%s\n", usage);
      return FALSE;
   }

   *netfilename = argv[1];

   if ( *netfilename == NULL || strstr(*netfilename, ".net") == NULL)
   {
      SCIPinfoMessage(scip, NULL, "No net-file supplied.\n\n");
      SCIPinfoMessage(scip, NULL, "%s\n", usage);
      return FALSE;
   }

   *scnfilename = argv[2];

   if ( *scnfilename == NULL || strstr(*scnfilename, ".scn") == NULL )
   {
      SCIPinfoMessage(scip, NULL, "No scn-file supplied.\n\n");
      SCIPinfoMessage(scip, NULL, "%s\n", usage);
      return FALSE;
   }

   /* check all remaining arguments */
   /*lint -e{850}*/
   for (i = 3; i < argc; ++i)
   {
      if ( strcmp(argv[i], "-a") == 0 )
      {
         probdata->algebraic = TRUE;
      }
      else if ( strcmp(argv[i], "-m") == 0 )
      {
         probdata->mixing = TRUE;
         probdata->algebraic = TRUE;
      }
      else if ( strcmp(argv[i], "-ml") == 0 )
      {
         probdata->linearMixing = TRUE;
         probdata->algebraic = TRUE;
      }
      else if ( strcmp( argv[i], "-nodeMu" ) == 0 )
      {
         probdata->nodeMu = TRUE;
      }
      else if ( strcmp(argv[i], "-fcm") == 0 )
      {
         probdata->flowConsMixing = TRUE;
         probdata->linearMixing = TRUE;
         probdata->algebraic = TRUE;
      }
       else if ( strcmp(argv[i], "-fcmn") == 0 )
      {
         probdata->flowConsMixingNode = TRUE;
         probdata->linearMixing = TRUE;
         probdata->algebraic = TRUE;
      }
      else if ( strcmp( argv[i], "-mol" ) == 0 )
      {
         probdata->mol = TRUE;
      }
      else if ( strcmp( argv[i], "-minLambda" ) == 0 )
      {
         probdata->minLambda = TRUE;
      }
      else if ( strcmp( argv[i], "-approx" ) == 0 )
      {
         probdata->approx = TRUE;
      }
      else if ( strcmp( argv[i], "-b" ) == 0 )
      {
         if ( probdata->charDiagram )
         {
            SCIPerrorMessage("Can not use both boxConstraintModel and detailed compressor station model.\n");
            return FALSE;
         }

         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "No cs-filename supplied.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         probdata->boxConstraintModel = TRUE;
         *csfilename = argv[i];
      }
      else if ( strcmp(argv[i], "-cD") == 0 )
      {
         if ( probdata->boxConstraintModel )
         {
            SCIPerrorMessage("Can not use both boxConstraintModel and detailed compressor station model.\n");
            return FALSE;
         }

         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "No cs-filename supplied.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         probdata->charDiagram = TRUE;
         *csfilename = argv[i];
      }
      else if ( strcmp(argv[i], "-writelsf") == 0 )
      {
         *writelsf = TRUE;
      }
      else if ( strcmp(argv[i], "-lsf") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == NULL ||  strstr(argv[i], ".lsf") == NULL )
         {
            SCIPinfoMessage(scip, NULL, "No solution file supplied.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         *lsfilename = argv[i];
      }
      else if ( strcmp(argv[i], "-minCompSum") == 0 )
      {
         probdata->minCompSum = TRUE;
         if ( probdata->minPressSum || probdata->noObjective || probdata->powerLoss || probdata->minSlack  || probdata->minSlackPerBound || probdata->maxFlow )
         {
            SCIPerrorMessage("More than one objective function given.\n");
            return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-binVarEQone") == 0 )
      {
         probdata->binVarEQone = TRUE;
      }
      else if ( strcmp(argv[i], "-minPressSum") == 0 )
      {
         probdata->minPressSum = TRUE;
         if ( probdata->minCompSum || probdata->noObjective || probdata->powerLoss || probdata->minSlack  || probdata->minSlackPerBound || probdata->maxFlow )
         {
            SCIPerrorMessage("More than one objective function given.\n");
            return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-minCompInc") == 0 )
      {
         probdata->minCompInc = TRUE;
         if ( probdata->minCompSum || probdata->noObjective || probdata->powerLoss || probdata->minSlack  || probdata->minSlackPerBound || probdata->maxFlow )
         {
            SCIPerrorMessage("More than one objective function given.\n");
            return FALSE;
         }
      }
      else if (strcmp(argv[i], "-noObjective") == 0)
      {
         probdata->noObjective = TRUE;
         if ( probdata->minCompSum || probdata->minPressSum || probdata->powerLoss || probdata->minSlack || probdata->minSlackPerBound || probdata->maxFlow )
         {
             SCIPerrorMessage("More than one objective function given.\n");
             return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-powerLoss") == 0 )
      {
         probdata->powerLoss = TRUE;
         if ( probdata->minCompSum || probdata->minPressSum || probdata->noObjective || probdata->minSlack || probdata->minSlackPerBound || probdata->maxFlow )
         {
             SCIPerrorMessage("More than one objective function given.\n");
             return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-minSlack") == 0 )
      {
         probdata->minSlack = TRUE;
         if ( probdata->minCompSum || probdata->minPressSum || probdata->noObjective || probdata->powerLoss ||  probdata->minSlackPerBound )
         {
             SCIPerrorMessage("More than one objective function given.\n");
             return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-minSlackPerBound") == 0 )
      {
         probdata->minSlackPerBound = TRUE;
         if ( probdata->minCompSum || probdata->minPressSum || probdata->noObjective || probdata->powerLoss || probdata->minSlack )
         {
             SCIPerrorMessage("More than one objective function given.\n");
             return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-maxFlow") == 0 )
      {
	 probdata->maxFlow = TRUE;
	 if ( probdata->minCompSum ||  probdata->minPressSum || probdata->noObjective || probdata->powerLoss )
	 {
	    SCIPerrorMessage("More than one objective function given.\n");
	    return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-relaxLowerBounds") == 0 )
      {
         probdata->relaxLowerBounds = TRUE;

         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "Error: No value for relaxing lower pressure bounds given.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         probdata->relaxLBvalue = atof(argv[i]);
         if ( SCIPisNegative (scip, probdata->relaxLBvalue) )
         {
            SCIPinfoMessage(scip, NULL, "Error: Value for relaxing lower pressure bounds has to be positive.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-relaxUpperBounds") == 0 )
      {
         probdata->relaxUpperBounds = TRUE;

         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "Error: No value for relaxing upper pressure bounds given.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         probdata->relaxUBvalue = atof(argv[i]);
         if ( SCIPisNegative(scip, probdata->relaxUBvalue) )
         {
            SCIPinfoMessage(scip, NULL, "Error: Value for relaxing upper pressure bounds has to be positive.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-relaxCVbounds") == 0 )
      {
         probdata->relaxCVbounds = TRUE;

         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "Error: No value for relaxing pressure bounds of controlvalves given.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         probdata->relaxCVvalue = atof(argv[i]);
         if ( SCIPisNegative (scip, probdata->relaxCVvalue) )
         {
            SCIPinfoMessage(scip, NULL, "Error: Value for relaxing controlvalve pressure bounds has to be positive.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
      }
      else if ( strcmp(argv[i], "-s") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "Error: No setting file name supplied.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         *settingsname = argv[i];
      }
      /* check for time limit */
      else if ( strcmp(argv[i], "-t") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "No time limit supplied.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         *timelimit = atof(argv[i]);
      }
      /* check for memory limit */
      else if ( strcmp(argv[i], "-mem") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "No memory limit supplied.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         *memlimit = atof(argv[i]);
      }
      /* check for node limit */
      else if ( strcmp(argv[i], "-n") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "No node limit supplied.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         *nodelimit = atol(argv[i]);
      }
      /* check for display frequency */
      else if ( strcmp(argv[i], "-d") == 0 )
      {
         ++i;
         if ( i >= argc || argv[i] == NULL )
         {
            SCIPinfoMessage(scip, NULL, "No display frequency supplied.\n\n");
            SCIPinfoMessage(scip, NULL, "%s\n", usage);
            return FALSE;
         }
         *dispfreq = atoi(argv[i]);
      }
      else if ( strcmp(argv[i], "-allCycles") == 0 )
      {
         probdata->allCycles = TRUE;
      }
      else
      {
	SCIPinfoMessage(scip, NULL, "Unknown argument %s.\n\n", argv[i]);
	SCIPinfoMessage(scip, NULL, "%s\n", usage);
        return FALSE;
      }
   }
   /*lint +e850*/

   if ( !(probdata->relaxLowerBounds) && !(probdata->relaxUpperBounds) && !(probdata->relaxCVbounds) && (probdata->minSlack || probdata->minSlackPerBound) )
   {
      SCIPinfoMessage(scip, NULL, "Error: Some pressure bounds have to be relaxed, such that one can minimize the slack.\n\n");
      SCIPinfoMessage(scip, NULL, "%s\n", usage);
      return FALSE;
   }

   return TRUE;
}


/** output statistics on constraint types */
static
SCIP_RETCODE outputConssTypes(
   SCIP*                 scip                /**< SCIP pointer */
   )
{
   SCIP_CONS** origconss;
   SCIP_CONSHDLR** conshdlrs;
   int norigconss;
   int nconshdlrs;
   int* nconss;
   int c;
   int h;

   norigconss = SCIPgetNOrigConss(scip);
   origconss = SCIPgetOrigConss(scip);

   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nconss, nconshdlrs) );

   /* loop over all constraints and constraint-handlers to count for each type the amount of original constraints */
   for (c = norigconss - 1; c >= 0; --c)
   {
      for (h = nconshdlrs - 1; h >= 0; --h)
      {
         if ( SCIPconsGetHdlr(origconss[c]) == conshdlrs[h] )
         {
            ++(nconss[h]);
            break;
         }
      }
      /* constraint handler should be found */
      assert( h >= 0 );
   }

   /* loop over all constraints handlers for printing the number of original constraints */
   for (h = 0; h < nconshdlrs; ++h)
   {
      if ( nconss[h] > 0 )
      {
         SCIPinfoMessage(scip, NULL, "%7d constraints of type <%s>\n", nconss[h], SCIPconshdlrGetName(conshdlrs[h]));
      }
   }
   SCIPinfoMessage(scip, NULL, "\n");

   SCIPfreeBufferArray(scip, &nconss);

   return SCIP_OKAY;
}

/** main function, which starts the solution of the gas problem */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP* scip = NULL;
   SCIP_PROBDATA* probdata = NULL;
   const char* netfilename;
   const char* scnfilename;
   const char* csfilename;
   const char* lsfilename;
   const char* settingsname;
   SCIP_Bool writelsf;
   SCIP_Real timelimit;
   SCIP_Real memlimit;
   SCIP_Longint nodelimit;
   int dispfreq;
   int errorcode = 0;

   /* initialize SCIP */
   SCIP_CALL_ERROR( SCIPcreate(&scip) );

   /* initialize probdata */
   SCIP_CALL_ERROR( GASinitProb(scip, &probdata) );

   /* parse command line arguments */
   if ( ! readArguments(scip, probdata, argc, argv, &netfilename, &scnfilename, &csfilename, &lsfilename, &writelsf,
         &settingsname, &timelimit, &memlimit, &nodelimit, &dispfreq) )
   {
      SCIP_CALL_ERROR( SCIPfree(&scip) );
      BMScheckEmptyMemory();
      return -1;
   }

   assert( netfilename != NULL );
   assert( scnfilename != NULL);

   /* include default SCIP plugins */
   SCIP_CALL_ERROR( SCIPincludeDefaultPlugins(scip) );

   /* include PipeODE constraint handler */
   if ( ! probdata->algebraic )
   {
      SCIP_CALL_ERROR( SCIPincludeConshdlrPipeODE(scip));
   }

   /* include heuristics */
   SCIP_CALL_ERROR( SCIPincludeHeurLPRound(scip) );
   SCIP_CALL_ERROR( SCIPincludeHeurGas(scip) );
   SCIP_CALL_ERROR( SCIPincludeHeurADM(scip) );

   /* include flow OBBT presolver */
   SCIP_CALL( SCIPincludePresolFlowOBBT(scip) );

   /* output version information */
   SCIPprintVersion(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPprintExternalCodes(scip, NULL);
   SCIPinfoMessage(scip, NULL, "\n");

   /* check for parameters */
   if ( settingsname != NULL )
   {
      if ( SCIPfileExists(settingsname) )
      {
         SCIPinfoMessage(scip, NULL, "Reading parameter file <%s> ...\n\n", settingsname);
         SCIP_CALL_ERROR( SCIPreadParams(scip, settingsname) );
      }
      else
      {
         SCIPwarningMessage(scip, "parameter file <%s> not found - using default parameters.\n", settingsname);
      }
   }

   if ( probdata->noFlowBinvars )
   {
      probdata->binaryFlowCons = FALSE;
      /* make sure no cycles are used */
      probdata->allCycles = FALSE;
      probdata->noCycles = TRUE;
   }

   if ( ! probdata->algebraic )
   {
      /* include display */
      SCIP_CALL_ERROR( SCIPincludeDispGradient(scip) );
   }

   /* output detailed information on model options and objective */
   SCIPinfoMessage(scip, NULL, "Solving the stationary gas transport problem.\n\n");

   if ( probdata->algebraic || probdata->boxConstraintModel || probdata->charDiagram || probdata->relaxLowerBounds
      || probdata->relaxUpperBounds || probdata->relaxCVbounds || probdata->noCycles || probdata->allCycles )
   {
      SCIPinfoMessage(scip, NULL, "Model options:\n");
      if ( probdata->algebraic )
         SCIPinfoMessage(scip, NULL, "\tUsing algebraic model for gas flow.\n");

      if ( probdata->mixing )
      {
         if ( probdata->nodeMu )
            SCIPinfoMessage( scip, NULL, "\tMixing algebraic model for gas flow with node mixing ratio.\n");
         else
            SCIPinfoMessage( scip, NULL, "\tMixing algebraic model for gas flow with arc mixing ratio.\n");
      }

      if ( probdata->linearMixing )
      {
         if ( probdata->nodeMu )
            SCIPinfoMessage( scip, NULL, "\tLinear mixing algebraic model for gas flow with node mixing ratio.\n");
         else
            SCIPinfoMessage( scip, NULL, "\tLinear mixing algebraic model for gas flow with arc mixing ratio.\n");
      }

      if ( probdata->approx )
         SCIPinfoMessage( scip, NULL, "\tUsing the algebraic model with only molar mass 1 in order to approximate the mixing ratio.\n");

      if ( probdata->binVarEQone )
         SCIPinfoMessage(scip, NULL, "\tFlow direction constraint z+ + z- == 1.\n");

      if (probdata->boxConstraintModel)
      {
         assert( csfilename != NULL);
         SCIPinfoMessage(scip, NULL, "\tUsing box contraint model for compressor stations.\n");
      }
      else if ( probdata->charDiagram )
      {
         assert( csfilename != NULL);
         SCIPinfoMessage(scip, NULL, "\tUsing characteristic diagrams for compressor machines.\n");
      }

      if ( probdata->relaxLowerBounds )
         SCIPinfoMessage(scip, NULL, "\tRelax lower pressure bounds by %f bar.\n", probdata->relaxLBvalue);

      if ( probdata->relaxUpperBounds )
         SCIPinfoMessage(scip, NULL, "\tRelax upper pressure bounds by %f bar.\n", probdata->relaxUBvalue);

      if ( probdata->relaxCVbounds )
         SCIPinfoMessage(scip, NULL, "\tRelax controlvalve pressure bounds by %f bar.\n", probdata->relaxCVvalue);

      if ( probdata->noFlowBinvars )
         SCIPinfoMessage(scip, NULL, "\tDo not use binary variables to indicate flow directions.\n");

      if ( probdata->allCycles && probdata->noCycles )
      {
         SCIPerrorMessage("Can not use both all cycles and no cycles to create noCircularFlowConstraints!");
         SCIP_CALL_ERROR( SCIPfree(&scip) );
         return -1;
      }
      else if ( probdata->allCycles )
         SCIPinfoMessage(scip, NULL, "\tFinding all cycles in the graph and add noCircularFlowConstraints.\n");
      else if ( probdata->noCycles )
         SCIPinfoMessage(scip, NULL, "\tDo not use cycles in the graph to add noCircularFlowConstraints.\n");

      SCIPinfoMessage(scip, NULL, "\n");
   }
   else
      SCIPinfoMessage(scip, NULL, "Using the default model.\n\n");

   if ( lsfilename != NULL )
      SCIPinfoMessage(scip, NULL, "Checking solution specified in file %s.\n\n", lsfilename);

   SCIPinfoMessage(scip, NULL, "Objective:\n\t");
   if ( probdata->minCompSum )
      SCIPinfoMessage(scip, NULL, "Minimizing sum of compressors.\n");
   else if ( probdata->minCompInc )
     SCIPinfoMessage(scip, NULL, "Minimizing the pressure increase by compressors.\n");
   else if (probdata->minPressSum)
      SCIPinfoMessage(scip, NULL, "Minimizing sum of pressure values.\n");
   else if ( probdata->noObjective )
      SCIPinfoMessage(scip, NULL, "Solving the feasibility problem.\n");
   else if ( probdata->powerLoss )
      SCIPinfoMessage(scip, NULL, "Minimizing the power loss.\n");
   else if ( probdata->minSlack )
      SCIPinfoMessage(scip, NULL, "Minimizing the difference to original pressure bounds.\n");
   else if ( probdata->minSlackPerBound )
      SCIPinfoMessage(scip, NULL, "Minimizing the sum of slack variables for relaxed pressure bounds.\n");
   else if ( probdata->maxFlow )
      SCIPinfoMessage(scip, NULL, "Maximizing the flow.\n");
   else if ( probdata->minLambda )
      SCIPinfoMessage(scip, NULL, "Minimizing the resistance on each pipe equal to maximizing the pressure.\n");
   else
      SCIPinfoMessage(scip, NULL, "Maximizing the sum of pressure values.\n");

   /* set time, node, memory, and display limit */
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      SCIP_CALL_ERROR( SCIPsetRealParam(scip, "limits/time", timelimit) );
   }

   if ( ! SCIPisInfinity(scip, memlimit) )
   {
      SCIP_CALL_ERROR( SCIPsetRealParam(scip, "limits/memory", memlimit) );
   }

   if ( nodelimit < SCIP_LONGINT_MAX )
   {
      SCIP_CALL_ERROR( SCIPsetLongintParam(scip, "limits/nodes", nodelimit) );
   }

   if ( dispfreq >= 0 )
   {
      SCIP_CALL_ERROR( SCIPsetIntParam(scip, "display/freq", dispfreq) );
   }

   /* print changed paramters */
   SCIPinfoMessage(scip, NULL, "\nChanged settings:\n");
   SCIP_CALL_ERROR( SCIPwriteParams(scip, NULL, FALSE, TRUE) );
   SCIPinfoMessage(scip, NULL, "\n");

   /* output tolerances for the pipe constraint handler */
   if ( ! probdata->algebraic )
   {
      SCIPinfoMessage(scip, NULL, "Feasibility tolerances:\n");
      SCIP_CALL_ERROR( SCIPwriteParam(scip, SCIPgetParam(scip, "pressure_feastol"), NULL, FALSE, FALSE) );
      SCIP_CALL_ERROR( SCIPwriteParam(scip, SCIPgetParam(scip, "flow_feastol"), NULL, FALSE, FALSE) );
      SCIPinfoMessage(scip, NULL, "\n");
   }

   /* read problem data */
   SCIP_CALL_ERROR( GAScreateProb(scip, probdata, netfilename, scnfilename, csfilename) );

   /* generate gas problem */
   SCIP_CALL_ERROR( GASgenerateModel(scip) );

   if ( lsfilename != NULL )
   {
      SCIPinfoMessage(scip, NULL, "\nSolution file name:\t%s\n", lsfilename);
      SCIP_CALL_ERROR( readLSF(scip, lsfilename) );
   }

   /* we don't need the hashtables any more */
   SCIPhashtableFree(&(probdata->network->Nodeshashtable));
   SCIPhashtableFree(&(probdata->network->Arcshashtable));

   probdata->network->Nodeshashtable = NULL;
   probdata->network->Arcshashtable = NULL;

   SCIPenableDebugSol(scip);

   /* print problem file for debugging */
#ifndef NDEBUG
   SCIP_CALL_ERROR( SCIPwriteOrigProblem(scip, "gas.cip", "cip", FALSE) );
   SCIP_CALL_ERROR( GASwriteNetwork(scip, probdata, "network.net", FALSE) );
#endif
   
   SCIP_CALL_ERROR( SCIPwriteOrigProblem( scip, "gas.cip", "cip", FALSE ) );

   SCIPinfoMessage(scip, NULL, "\nOriginal problem has %d variables (%d bin, %d int, %d impl, %d cont) and %d constraints.\n",
      SCIPgetNOrigVars(scip), SCIPgetNOrigBinVars(scip), SCIPgetNOrigIntVars(scip), SCIPgetNOrigImplVars(scip), SCIPgetNOrigContVars(scip), SCIPgetNOrigConss(scip));

   /* output constraint types */
   SCIP_CALL_ERROR( outputConssTypes(scip) );

   /* print presolve file for debugging*/
   SCIP_CALL_ERROR( SCIPpresolve(scip) );
#ifndef NDEBUG
   SCIP_CALL_ERROR( SCIPwriteTransProblem(scip, "pregas.cip", "cip", FALSE) );
#endif
 SCIP_CALL_ERROR( SCIPwriteTransProblem(scip, "pregas.cip", "cip", FALSE) );
#ifdef FLOWBOUNDS_OUTPUT
   {
      int nfixedflow = 0;
      int nfixeddir = 0;
      int counter = 0;
      int counter2 = 0;
      SCIP_Real qlbsum = 0.0;
      SCIP_Real qubsum = 0.0;
      SCIP_Real qlbsum2 = 0.0;
      SCIP_Real qubsum2 = 0.0;
      int i;

      /* SCIPinfoMessage(scip, NULL, "\n"); */
      for (i = 0; i < probdata->network->numarcs; ++i)
      {
         GAS_Arc* arc;
         SCIP_VAR* flowvar;
         SCIP_VAR* transflowvar;
         SCIP_Real qlb;
         SCIP_Real qub;

         arc = &probdata->network->arcs_ptr[i];

         if ( arc->type != PIPE )
            continue;

         flowvar = arc->flowvar;
         transflowvar = SCIPvarGetTransVar(flowvar);

         qlb = SCIPcomputeVarLbGlobal(scip, transflowvar);
         qub = SCIPcomputeVarUbGlobal(scip, transflowvar);

         qlbsum += qlb;
         qubsum += qub;
         ++counter;

         /* SCIPinfoMessage(scip, NULL, "flow bounds on arc %30s\t LB %8.4f \t UB %8.4f\n", arc->id, qlb, qub); */
         if ( SCIPisFeasEQ(scip, qlb, qub) )
         {
            ++nfixedflow;
         }
         else
         {
            qlbsum2 += qlb;
            qubsum2 += qub;
            ++counter2;

            if ( SCIPisFeasLE(scip, qub, 0.0) || SCIPisFeasGE(scip, qlb, 0.0) )
               ++nfixeddir;
         }
      }

      SCIPinfoMessage(scip, NULL, "\nFlow bound statistics:\n");
      SCIPinfoMessage(scip, NULL, "\tNumber of pipes:            %3d\n", counter);
      SCIPinfoMessage(scip, NULL, "\tPipes with fixed flow:      %3d\n", nfixedflow);
      SCIPinfoMessage(scip, NULL, "\tPipes with fixed direction: %3d\n", nfixeddir);
      SCIPinfoMessage(scip, NULL, "\tArithmetic mean lower bound:        %7.2f\n", qlbsum / ((SCIP_Real) counter));
      SCIPinfoMessage(scip, NULL, "\tArithmetic mean upper bound:        %7.2f\n", qubsum / ((SCIP_Real) counter));
      SCIPinfoMessage(scip, NULL, "\tArithmetic mean lower bound without fixed flows: %7.2f\n", qlbsum2 / ((SCIP_Real) counter2));
      SCIPinfoMessage(scip, NULL, "\tArithmetic mean upper bound without fixed flows: %7.2f\n", qubsum2 / ((SCIP_Real) counter2));
      SCIPinfoMessage(scip, NULL, "\n");
   }
#endif

   /* solve the model */
   SCIP_CALL_ERROR( SCIPsolve(scip) );

   /* print statistics */
   SCIP_CALL_ERROR( SCIPprintStatistics(scip, NULL) );

#if 0
   /* print branching statistics */
   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL_ERROR( SCIPprintBranchingStatistics(scip, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL_ERROR( SCIPprintBestSol(scip, NULL, FALSE) );

   /* possibly save solution */
   if ( writelsf )
   {
      SCIP_SOL* sol = SCIPgetBestSol(scip);
      SCIP_CALL_ERROR( writeLSF(scip, sol) );
   }

   SCIP_CALL_ERROR( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return errorcode;
}
