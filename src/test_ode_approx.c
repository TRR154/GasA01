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

/**@file   test_ode_approx.c
 * @brief  test different ODE approximation codes
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "ode_functions_gas.h"
#include "data_structs_cons_pipe_ode.h"


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


/** main function, which starts the solution of the gas problem */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP* scip = NULL;
   struct SCIP_ConsData consdata;
   struct SCIP_ConshdlrData conshdlrdata;
   SCIP_Real eps = 1e-3;
   SCIP_Real val;
   SCIP_Real valeps;
   SCIP_Real dp;
   SCIP_Real dq;

   /* initialize SCIP */
   SCIP_CALL_ERROR( SCIPcreate(&scip) );

   /* initialize fake cons(hdlr)data */
   conshdlrdata.flow_feastol = 1e-6;
   conshdlrdata.machbound = INT_MAX;

   consdata.L = 1000.0;
   consdata.N = 10000;
   consdata.A = 0.5;
   consdata.D = 1.0;
   consdata.lambda = 0.0105409;
   consdata.c = 330.0;

   consdata.Asq = consdata.A * consdata.A;
   consdata.csq = 330.0 * 330.0;
   consdata.lcsqd = consdata.lambda * consdata.csq / consdata.D;

   printf("Results for input pressure %g and flow %g with N = %d steps:\n", 10.0, 1000.0, consdata.N);
   val = computeModifiedEulerScheme(scip, &conshdlrdata, &consdata, 10.0, 1000.0, -1);
   printf("Modified Euler output value:\t %.12g\n", val);

   val = computeTrapezoidalScheme(scip, &conshdlrdata, &consdata, 10.0, 1000.0, -1);
   printf("Trapazoidal output value:\t %.12g\n", val);

   val = computeRK4(scip, &conshdlrdata, &consdata, 10.0, 1000.0, -1);
   printf("RK4 output value:\t\t %.12g\n", val);

   /* -----------------------------------------------------------------*/
   printf("\n\n");
   printf("Results for input flow %g:\n", 1000.0);

   val = computeModifiedEulerSchemeWithGradients(scip, &conshdlrdata, &consdata, 10.0, 1000.0, &dp, &dq);
   printf("Modified Euler output value at p = % g:\t %.12g, gradient: (%g, %g)\n", 10.0, val, dp, dq);

   valeps = computeModifiedEulerScheme(scip, &conshdlrdata, &consdata, 10.0 - eps, 1000.0, -1);
   printf("Modified Euler output value at p = %g:\t %.12g\n", 10.0 - eps, valeps);
   val = computeModifiedEulerScheme(scip, &conshdlrdata, &consdata, 10.0, 1000.0, -1);
   printf("Modified Euler output value at p = %g:\t\t %.12g\n", 10.0, val);
   printf("Approximated gradient w.r.t. p: %g\n", (val - valeps) / eps);

   printf("\n");
   printf("Results for input pressure %g:\n", 10.0);

   val = computeModifiedEulerSchemeWithGradients(scip, &conshdlrdata, &consdata, 10.0, 1000.0, &dp, &dq);
   printf("Modified Euler output value at p = % g:\t %.12g, gradient: (%g, %g)\n", 10.0, val, dp, dq);

   valeps = computeModifiedEulerScheme(scip, &conshdlrdata, &consdata, 10.0, 1000.0 - eps, -1);
   printf("Modified Euler output value at q = %g:\t %.12g\n", 1000.0 - eps, valeps);
   val = computeModifiedEulerScheme(scip, &conshdlrdata, &consdata, 10.0, 1000.0, -1);
   printf("Modified Euler output value at q = %g:\t %.12g\n", 1000.0, val);
   printf("Approximated gradient w.r.t. q: %g\n", (val - valeps) / eps);

   /* -----------------------------------------------------------------*/
   printf("\n\n");
   printf("Results for input flow %g:\n", 1000.0);

   val = computeRK4WithGradients(scip, &conshdlrdata, &consdata, 10.0, 1000.0, &dp, &dq);
   printf("RK4 output value at p = %g:\t\t %.12g, gradient: (%g, %g)\n", 10.0, val, dp, dq);

   valeps = computeRK4(scip, &conshdlrdata, &consdata, 10.0 - eps, 1000.0, -1);
   printf("RK4 output value at p = %g:\t\t %.12g\n", 10.0 - eps, valeps);
   val = computeRK4(scip, &conshdlrdata, &consdata, 10.0, 1000.0, -1);
   printf("RK4 output value at p = %g:\t\t %.12g\n", 10.0, val);
   printf("Approximated gradient w.r.t. p: %g\n", (val - valeps) / eps);


   printf("\n");
   printf("Results for input pressure %g:\n", 10.0);

   val = computeRK4WithGradients(scip, &conshdlrdata, &consdata, 10.0, 1000.0, &dp, &dq);
   printf("RK4 output value at q = %g:\t\t %.12g, gradient: (%g, %g)\n", 1000.0, val, dp, dq);

   valeps = computeRK4(scip, &conshdlrdata, &consdata, 10.0, 1000.0 - eps, -1);
   printf("RK4 output value at q = %g:\t %.12g\n", 1000.0 - eps, valeps);
   val = computeRK4(scip, &conshdlrdata, &consdata, 10.0, 1000.0, -1);
   printf("RK4 output value at q = %g:\t\t %.12g\n", 1000.0, val);
   printf("Approximated gradient w.r.t. q: %g\n", (val - valeps) / eps);

   SCIP_CALL_ERROR( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return 0;
}
