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

/**@file   ode_functions_gas.c
 * @brief  functions used for the relaxation of the stationary, isothermal Euler equations
 * @author Alexandra Stern
 * @author Oliver Habeck
 */

#include "ode_functions_gas.h"
#include <scip/misc.h>
#include "comparisons_gas.h"
#include "data_structs_cons_pipe_ode.h"

/* #define SCIP_DEBUG */

/** compute function for right hand side of ODE */
static
SCIP_Real compute_h(
   SCIP_Real             p,                  /**< pressure in Pascal */
   SCIP_Real             q,                  /**< massflow in kg per second */
   SCIP_Real             Asq,                /**< squared area of pipe in square meter */
   SCIP_Real             lcsqd,              /**< lambda*c^2/D */
   SCIP_Real             csq                 /**< squared speed of sound  */
   )
{
   SCIP_Real h;
   SCIP_Real qsq;
   SCIP_Real psp;

   psp =  p * p;
   qsq =  q * q;

   h = (- lcsqd * qsq * p)/(2.0 * (Asq * psp - csq * qsq));

   if ( h >= 0 )
   {
      SCIPerrorMessage("Die Funktion h hat das falsche Vorzeichen.\n");
   }

   return h;
}

/** compute derivative of h w.r.t. p */
static
SCIP_Real compute_dph(
   SCIP_Real             p,                  /**< pressure in Pascal */
   SCIP_Real             q,                  /**< massflow in kg per second */
   SCIP_Real             Asq,                /**< squared area of pipe in square meter */
   SCIP_Real             lcsqd,              /**< lambda*c^2/D */
   SCIP_Real             csq                 /**< squared speed of sound  */
   )
{
   SCIP_Real dph;
   SCIP_Real qsq;
   SCIP_Real psp;

   psp =  p * p;
   qsq =  q * q;

   dph = (lcsqd * qsq * ( Asq * psp + csq * qsq))/(2.0 * (Asq * psp - csq * qsq) * (Asq * psp - csq * qsq));

   return dph;
}

/** compute derivative of h w.r.t. q */
static
SCIP_Real compute_dqh(
   SCIP_Real             p,                  /**< pressure in Pascal */
   SCIP_Real             q,                  /**< massflow in kg per second */
   SCIP_Real             Asq,                /**< squared area of pipe in square meter */
   SCIP_Real             lcsqd,              /**< lambda*c^2/D */
   SCIP_Real             csq                 /**< squared speed of sound  */
   )
{
   SCIP_Real dqh;
   SCIP_Real qsq;
   SCIP_Real psp;

   psp =  p * p;
   qsq =  q * q;

   dqh = - ((lcsqd * Asq * psp * p *  q) / (Asq * psp - csq * qsq)) / (Asq * psp - csq * qsq);

   return dqh;
}

/** compute step length
 *
 *  Note: all delta_x have the same length
 */
static
SCIP_Real compute_delta_x(
   SCIP_Real             L,                  /**< length of pipe */
   int                   N                   /**< discretization number */
   )
{
   SCIP_Real delta_x;

   delta_x = L/N;

   return delta_x;
}

/** computes the Modified Euler Scheme (= Explicit Midpoint Rule) either in direction of flow (dir==1) or against
 *  (dir==-1).
 *
 *  Returns -1 if either the pressure falls below the
 *  critical bound  machbound * c * q/A (dir == 1) or if something goes terribly wrong and
 *  the pressure decreases instead of increases when dir == -1.
 */
SCIP_Real computeModifiedEulerScheme(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of the constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   int                   dir                 /**< direction of computation */
   )
{
   SCIP_Real p_i;
   SCIP_Real p_me;
   SCIP_Real delta_x;
   SCIP_Real p_e;
   SCIP_Real p_bound;
   int       i;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( consdata != NULL );
   assert( dir == -1 || dir == 1 );

   if ( flowIsZero(q, conshdlrdata->flow_feastol) )
   {
      return p;
   }

   assert( flowIsPositive(q, conshdlrdata->flow_feastol) );

   if ( dir == 1 )
   {
      p_bound = 5.0 * consdata->c * q / ( ((SCIP_Real) conshdlrdata->machbound) * consdata->A);
   }
   else
   {
      p_bound = p * 1e5;
   }

   /* initialize pressure in pascal */
   p_me = p * 1e5;

   delta_x = dir * compute_delta_x(consdata->L, consdata->N);

   for (i = 0; i < consdata->N; ++i)
   {
      p_i = p_me;

      /* euler step within h */
      p_e = p_i + (delta_x / 2) * compute_h(p_i, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      /* SCIPisLE(scip, p_e, p_bound)? */
      if ( SCIPisFeasLT(scip, p_e, p_bound) )
      {
         return -1.0;
      }

      /* modified euler step */
      p_me = p_i + delta_x * compute_h(p_e, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      if ( SCIPisFeasLT(scip, p_me, p_bound) )
      {
         return -1.0;
      }
   }

   /* turn unit of pressure back to bar */
   p_me /= 1e5;

   return p_me;
}

/** computes the Modified Euler Scheme (= Explicit Midpoint Rule) against the direction of the flow and the gradient of
 *  the "incoming" pressure with respect to the flow and the "outgoing" pressure
 */
SCIP_Real computeModifiedEulerSchemeWithGradients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   SCIP_Real*            dp,                 /**< pointer to the partial derivative dp with no unit */
   SCIP_Real*            dq                  /**< pointer to the partial derivative dq in  bar * (kg / s)^{-1} */
   )
{
   SCIP_Real p_i;
   SCIP_Real p_me;
   SCIP_Real delta_x;
   SCIP_Real p_e;
   SCIP_Real dp_step;
   SCIP_Real dq_step;
   int i;
   SCIP_Real dph_pe;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( consdata != NULL );
   assert( dp != NULL );
   assert( dq != NULL );

   if ( flowIsZero(q, conshdlrdata->flow_feastol) )
   {
      *dp = 1.0;
      *dq = 0.0;
      return p;
   }

   assert( flowIsPositive(q, conshdlrdata->flow_feastol) );

   /* initialize with pressure in pascal */
   p_me = p * 1e5;
   /* initialize gradients */
   *dp = 1.0;
   *dq = 0.0;

   delta_x = - compute_delta_x(consdata->L, consdata->N);

   for (i = 0; i < consdata->N; ++i)
   {
      p_i = p_me;

      /* Euler step within h */
      p_e = p_i + (delta_x / 2) * compute_h(p_i, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      /* modified Euler step */
      p_me = p_i + delta_x * compute_h(p_e, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      dph_pe = compute_dph(p_e, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      /* derivatives of one step */
      dp_step = 1.0 + delta_x * dph_pe * (1.0 + (delta_x / 2.0) * compute_dph(p_i, q, consdata->Asq, consdata->lcsqd, consdata->csq));
      dq_step = delta_x * compute_dqh(p_e, q, consdata->Asq, consdata->lcsqd, consdata->csq) + (delta_x *delta_x / 2.0) * dph_pe * compute_dqh(p_i, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      /* update dp */
      (*dp) = dp_step * (*dp);

      /* update dq */
      (*dq) = dq_step + dp_step * (*dq);
   }

   /* turn back to pressure in bar */
   p_me /= 1e5;
   /* *dq has unit Pa * (kg / s)^{-1} */
   (*dq) /= 1e5;

   return p_me;
}

/** newton method for forward computation of the trapezoidal scheme
 *  this method is special tailored to produce underestimators of the pressure.
 *
 *  Therefore write the scheme as function
 *  \f[
 *  F(p_{i-1}, p_i, q) = p_{i-1} - p_i + \frac{{\Delta} x}{2} ( h(p_{i-1}, q) + h(p_i, q) )
 *  \f]
 *  only for positive q.
 */
static
SCIP_Real forwardNewton(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in Pascal*/
   SCIP_Real             q,                  /**< massflow in kg per second */
   SCIP_Real             tol,                /**< feasibility tolerance */
   int*                  tcount              /**< counter of total newton iterations per call of computeTrapezoidalScheme */
   )
{
   SCIP_Real p_out;
   SCIP_Real F;
   SCIP_Real F_p;
   SCIP_Real h_p;
   SCIP_Real dF;
   SCIP_Real delta_p;
   SCIP_Real delta_x;
   int counter = 0;
   SCIP_Real numerator;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( consdata != NULL );
   assert( flowIsPositive(q, conshdlrdata->flow_feastol) );

   numerator = (SCIP_Real) conshdlrdata->machbound;

   p_out = 5.0 * consdata->c * q / (numerator * consdata->A );
   assert( SCIPisFeasGE(scip, p, p_out) );

   delta_x = compute_delta_x(consdata->L, consdata->N);

   /* the p_{i-1} part of F(p_{i-1}, p_i, q) does not change; define it as F_p =>; don't recompute h(p_{i-1}, q) */
   h_p = compute_h(p, q, consdata->Asq, consdata->lcsqd, consdata->csq);
   F_p = p + delta_x * h_p / 2.0;

   /* first we check if there is a solution of F = 0 between 5cq/(numerator * A) and p */
   F = (F_p - p_out) + delta_x * compute_h(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq) / 2.0;
   if ( F <= 0.0 )
   {
      return -1.0;
   }

   /* we start the newton method with the previous pressure value */
   p_out = p;
   F = delta_x * h_p;
   dF = - 1.0 + (delta_x / 2) * compute_dph(p, q, consdata->Asq, consdata->lcsqd, consdata->csq);
   delta_p = - F / dF;

   /* since p is bigger than the solution should be, F has to be negative */
   assert( F <= 0.0 );

   /* now we first use the complete newton steps on the right side of the solution to get close to F(p_out) = 0 */
   while ( F < - 1.0 )
   {
      p_out += delta_p;
      F = (F_p - p_out) + delta_x * compute_h(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq) / 2.0;
      dF = - 1.0 + (delta_x / 2.0) * compute_dph(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      delta_p = - F / dF;
      ++(*tcount);
   }

   /* next we push out current p_out to the left side of the solution */
   while ( F < 0.0 )
   {
      /* if delta_p is very small we may get into an infinite loop, therefore we accept making a bigger error by subtracting .5 bar to p_out */
      if ( counter == 5 )
      {
         assert( SCIPisGT(scip, delta_p, - 0.5 ) );
         delta_p = - 0.5;
      }

      p_out += delta_p;
      F = (F_p - p_out) + delta_x *  compute_h(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq) / 2.0;
      ++counter;
      ++(*tcount);
   }

   if ( SCIPisFeasLT(scip, numerator * consdata->A * p_out, 5.0 * consdata->c * q) )
   {
      p_out = 5.0 * consdata->c * q / (numerator * consdata->A );
      F = (F_p - p_out) + delta_x * compute_h(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq) / 2.0;
   }

   /* next we use the newton method with half of the newton steps to assure, that we keep the sign of F */
   dF = - 1.0 + (delta_x / 2.0) * compute_dph(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq);
   delta_p = - F / dF;
   while ( F > tol )
   {
      p_out += (delta_p / 2.0);
      F = (F_p - p_out) + delta_x * compute_h(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq) / 2.0;
      assert( F >= 0.0 );
      dF = - 1.0 + (delta_x / 2.0) * compute_dph(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      delta_p = - F / dF;
      ++(*tcount);
   }

   assert( F >= 0.0 );
   assert( p_out <= p );

   return p_out;
}


/** newton method for backward computation of the trapezoidal scheme
 *  this method is special tailored to produce overestimators of the pressure
 *
 *  Therefore write the scheme as function
 *  \f[
 *  F(p_{i-1}, p_i, q) = p_{i-1} - p_i - \frac{{\Delta} x}{2} ( h(p_{i-1}, q) + h(p_i, q) )
 *  \f]
 *  only for positive q and positive stepsizes.
 */
static
SCIP_Real backwardNewton(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   SCIP_Real             tol,                /**< feasibility tolerance */
   int*                  tcount              /**< counter of total newton iterations per call of computeTrapezoidalScheme */
   )
{
   SCIP_Real p_out = p;
   SCIP_Real F;
   SCIP_Real h_p;
   SCIP_Real F_p;
   SCIP_Real dF;
   SCIP_Real delta_p;
   SCIP_Real delta_x;
   int counter = 0;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( consdata != NULL );
   assert( flowIsPositive(q, conshdlrdata->flow_feastol) );
   assert( SCIPisFeasGE(scip, ((SCIP_Real) conshdlrdata->machbound) * consdata->A * p, 5.0 * consdata->c * q) );

   delta_x = compute_delta_x(consdata->L, consdata->N);

   /* the p_{i-1} part of F(p_{i-1}, p_i, q) does not change; define it as F_p =>; don't recompute h(p_{i-1}, q) */
   h_p = compute_h(p, q, consdata->Asq, consdata->lcsqd, consdata->csq);
   F_p = p - delta_x * h_p / 2.0;

   /* we start the newton method with the previous pressure value */
   F = - delta_x * h_p;
   dF = - 1.0 - (delta_x / 2) * compute_dph(p, q, consdata->Asq, consdata->lcsqd, consdata->csq);
   delta_p = - F / dF;

   /* since p is smaller than the solution should be, F has to be positive */
   assert( F >= 0.0 );

   /* now we first use the complete newton steps "left" of the solution to get close to F(p_out) = 0 */
   while ( F > 1.0 )
   {
      p_out += delta_p;
      F = (F_p - p_out) - delta_x * compute_h(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq) / 2.0;
      dF = - 1.0 - (delta_x / 2.0) * compute_dph(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      delta_p = - F / dF;
      ++(*tcount);
   }

   /* next we push out current p_out to the right side of the solution */
   while ( F > 0.0 )
   {
      /* if delta_p is very small we may get into an infinite loop, therefore we accept making a bigger error by adding .5 bar to p_out */
      if ( counter == 5 )
      {
         assert( SCIPisLT(scip, delta_p, 0.5 ) );
         delta_p = 0.5;
      }

      p_out += delta_p;
      F = (F_p - p_out) - delta_x * compute_h(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq) / 2.0;
      ++counter;
      ++(*tcount);
   }

   /* next we use the newton method with two-third of the newton steps to assure, that we keep the sign of F */
   dF = - 1.0 - (delta_x / 2.0) * compute_dph(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq);
   delta_p = - F / dF;
   while ( - F > tol )
   {
      p_out += (delta_p / 3.0) * 2.0;
      F = (F_p - p_out) - delta_x * compute_h(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq) / 2.0;
      assert( F <= 0.0 );
      dF = - 1.0 - (delta_x / 2.0) * compute_dph(p_out, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      delta_p = - F / dF;
      ++(*tcount);
   }

   assert( F <= 0.0 );
   assert( p_out >= p );

   return p_out;
}

/** compute the pressure bound via the trapezoidal scheme */
SCIP_Real computeTrapezoidalScheme(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   int                   dir                 /**< direction of computation */
   )
{
   SCIP_Real p_out;
   SCIP_Real tol = 1e-5;
   int i;
   int tcount = 0;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( consdata != NULL );
   assert( dir == -1 || dir == 1 );

   if ( flowIsZero(q, conshdlrdata->flow_feastol) )
   {
      return p;
   }

   assert( flowIsPositive(q, conshdlrdata->flow_feastol) );

   if ( SCIPisFeasLT(scip, (((SCIP_Real) conshdlrdata->machbound) * consdata->A * 1e5 * p),(5.0 * consdata->c * q)) )
   {
      return -1.0;
   }

   /* initialize with pressure in pascal */
   p_out = p * 1e5;

   for (i = 0; i < consdata->N; ++i)
   {
      if ( dir == 1 )
      {
         p_out = forwardNewton(scip, conshdlrdata, consdata, p_out, q, tol, &tcount);
      }
      else
      {
         p_out = backwardNewton(scip, conshdlrdata, consdata, p_out, q, tol, &tcount);
      }
      /* if the newton step fails to produce a solution,
       * that is, either there is no solution or the discretization
       * is too coarse. */
      if ( p_out < -0.5 )
      {
         return -1.0;
      }
   }

   assert( i > 0 );
   if ( dir == -1 )
   {
      SCIPdebugMessage("average #iterations in backwardNewton: %f \n", (double)tcount/ (double)i);
   }
   else
   {
      SCIPdebugMessage("average #iterations in forwardNewton: %f \n", (double)tcount/ (double)i);
   }

   /* turn back to bar */
   p_out /= 1e5;

   return p_out;
}

/** computes RK4 solutions either in direction of flow (dir==1) or against (dir==-1). */
SCIP_Real computeRK4(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of the constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   int                   dir                 /**< direction of computation */
   )
{
   SCIP_Real p_i;
   SCIP_Real delta_x;
   int       i;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( consdata != NULL );
   assert( dir == -1 || dir == 1 );

   if ( flowIsZero(q, conshdlrdata->flow_feastol) )
   {
      return p;
   }

   assert( flowIsPositive(q, conshdlrdata->flow_feastol) );

   /* initialize pressure in pascal */
   p_i = p * 1e5;

   delta_x = (SCIP_Real) dir * compute_delta_x(consdata->L, consdata->N);

   for (i = 0; i < consdata->N; ++i)
   {
      SCIP_Real k1;
      SCIP_Real k2;
      SCIP_Real k3;
      SCIP_Real k4;

      /* compute RK4 iteration */
      k1 = compute_h(p_i, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      k2 = compute_h(p_i + delta_x / 2.0 * k1, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      k3 = compute_h(p_i + delta_x / 2.0 * k2, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      k4 = compute_h(p_i + delta_x * k3, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      p_i = p_i + (delta_x / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
   }

   /* turn unit of pressure back to bar */
   p_i /= 1e5;

   return p_i;
}

/** computes RK4 solutions against the direction of the flow and the gradient of the "incoming" pressure with respect to
 *  the flow and the "outgoing" pressure.
 */
SCIP_Real computeRK4WithGradients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< data of constraint handler */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Real             p,                  /**< pressure in bar */
   SCIP_Real             q,                  /**< massflow in kg per second */
   SCIP_Real*            dp,                 /**< pointer to the partial derivative dp with no unit */
   SCIP_Real*            dq                  /**< pointer to the partial derivative dq in  bar * (kg / s)^{-1} */
   )
{
   SCIP_Real p_i;
   SCIP_Real delta_x;
   int i;

   assert( scip != NULL );
   assert( conshdlrdata != NULL );
   assert( consdata != NULL );
   assert( dp != NULL );
   assert( dq != NULL );

   /* initialize gradients */
   *dp = 1.0;
   *dq = 0.0;

   if ( flowIsZero(q, conshdlrdata->flow_feastol) )
      return p;

   assert( flowIsPositive(q, conshdlrdata->flow_feastol) );

   /* initialize with pressure in pascal */
   p_i = p * 1e5;

   delta_x = - compute_delta_x(consdata->L, consdata->N);

   for (i = 0; i < consdata->N; ++i)
   {
      SCIP_Real k1;
      SCIP_Real k2;
      SCIP_Real k3;
      SCIP_Real k4;
      SCIP_Real dpk1;
      SCIP_Real dpk2;
      SCIP_Real dpk3;
      SCIP_Real dpk4;
      SCIP_Real dqk1;
      SCIP_Real dqk2;
      SCIP_Real dqk3;
      SCIP_Real dqk4;
      SCIP_Real dph1;
      SCIP_Real dph2;
      SCIP_Real dph3;
      SCIP_Real dph4;
      SCIP_Real dqh1;
      SCIP_Real dqh2;
      SCIP_Real dqh3;
      SCIP_Real dqh4;
      SCIP_Real ddp;
      SCIP_Real ddq;

      /* compute RK4 iteration */
      k1 = compute_h(p_i, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      k2 = compute_h(p_i + delta_x / 2.0 * k1, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      k3 = compute_h(p_i + delta_x / 2.0 * k2, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      k4 = compute_h(p_i + delta_x * k3, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      /* compute derivative of rhs w.r.t. pressure */
      dph1 = compute_dph(p_i, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      dph2 = compute_dph(p_i + delta_x / 2.0 * k1, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      dph3 = compute_dph(p_i + delta_x / 2.0 * k2, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      dph4 = compute_dph(p_i + delta_x * k3, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      /* compute derivative of rhs w.r.t. flow */
      dqh1 = compute_dqh(p_i, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      dqh2 = compute_dqh(p_i + delta_x / 2.0 * k1, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      dqh3 = compute_dqh(p_i + delta_x / 2.0 * k2, q, consdata->Asq, consdata->lcsqd, consdata->csq);
      dqh4 = compute_dqh(p_i + delta_x * k3, q, consdata->Asq, consdata->lcsqd, consdata->csq);

      /* compute derivative of pressure iteration w.r.t. initial pressure via chain rule */
      dpk1 = dph1 * (*dp);
      dpk2 = dph2 * ((*dp) + delta_x / 2.0 * dpk1);
      dpk3 = dph3 * ((*dp) + delta_x / 2.0 * dpk2);
      dpk4 = dph4 * ((*dp) + delta_x * dpk3);

      ddp = (*dp) + delta_x / 6.0 * (dpk1 + 2.0 * dpk2 + 2.0 * dpk3 + dpk4);

      /* compute derivative of pressure iteration w.r.t. initial flow via chain rule */
      dqk1 = dph1 * (*dq) + dqh1;
      dqk2 = dph2 * ((*dq) + delta_x / 2.0 * dqk1) + dqh2;
      dqk3 = dph3 * ((*dq) + delta_x / 2.0 * dqk2) + dqh3;
      dqk4 = dph4 * ((*dq) + delta_x * dqk3) + dqh4;

      ddq = (*dq) + delta_x / 6.0 * (dqk1 + 2.0 * dqk2 + 2.0 * dqk3 + dqk4);

      /* update */
      p_i = p_i + (delta_x / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
      (*dp) = ddp;
      (*dq) = ddq;
   }

   /* turn back to pressure in bar */
   p_i /= 1e5;

   /* *dq has unit Pa * (kg / s)^{-1} */
   (*dq) /= 1e5;

   return p_i;
}
