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

/**@file   myxmldefs.h
 * @brief  glue code to match xml functions between SCIP 10 and lower versions
 * @author Marc Pfetsch
 */

#ifndef __MYXMLDEFS__
#define __MYXMLDEFS__

/* some code to make code compatible with SCIP 10 */
#if ( SCIP_VERSION < 1000 )
#define SCIPxmlGetAttrval xmlGetAttrval
#define SCIPxmlFindNodeMaxdepth xmlFindNodeMaxdepth
#define SCIPxmlFirstChild xmlFirstChild
#define SCIPxmlNextSibl xmlNextSibl
#define SCIPxmlGetName xmlGetName
#define SCIPxmlFreeNode xmlFreeNode
#define SCIPxmlProcess xmlProcess
#endif

#endif
