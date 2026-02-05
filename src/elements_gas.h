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

/**@file   elements_gas.h
 * @brief  handling of data needed for solving stationary gas transport problems
 * @author Marc Pfetsch
 * @author Alexandra Stern
 * @author Imke Joormann
 * @author Oliver Habeck
 * @author PAscal Boerner
 */

#include <scip/def.h>

struct GAS_arc;


/** types of nodes */
enum node_type
{
   UNKNOWNNODETYPE = 0,                      /**< unknown type (for debugging) */
   ENTRY           = 1,                      /**< entry */
   EXIT            = 2,                      /**< exit */
   INNODE          = 3                       /**< innode */
};
typedef enum node_type Node_type;

/** node of gas network */
struct GAS_node
{
   char*                 id;                           /**< id of node */
   int                   nodeposition;                 /**< position in the nodes_ptr array in the network struct */
   SCIP_Real             geoWGS84Long;                 /**< geographic longitude coordinate of node */
   SCIP_Real             geoWGS84Lat;                  /**< geographic latitude coordinate of node */
   SCIP_Real             x;                            /**< x coordinate of node */
   SCIP_Real             y;                            /**< y coordinate of node */
   SCIP_Real             height;                       /**< height of node */
   SCIP_Real             pressureMin;                  /**< minimal pressure */
   SCIP_Real             pressureMax;                  /**< maximal pressure */
   SCIP_Real             netPressureMin;               /**< minimal pressure given in net-file */
   SCIP_Real             netPressureMax;               /**< maximal pressure given in net-file */
   SCIP_Real             arcPressureMax;               /**< maximal allowed pressure of incident pipes */
   struct GAS_arc*       inarcs;                       /**< ingoing arcs: depends on the orientation of a pipe */
   struct GAS_arc*       outarcs;                      /**< outgoing arcs: depends on the orientation of a pipe */
   int                   numneighbors;                 /**< number of adjacent nodes */
   int                   numincidentarcs;              /**< number of incident arcs */
   Node_type             type;                         /**< "entry", "exit" or "innode": information from scn file */
   SCIP_Real             scnPressureMin;               /**< information from scn file: */
   SCIP_Real             scnPressureMax;               /**< information from scn file */
   SCIP_Real             flow;                         /**< information from scn file */
   SCIP_Real             scnFlowMin;                   /**< information from scn file */
   SCIP_Real             scnFlowMax;                   /**< information from scn file */
   SCIP_VAR*             pressurevar;                  /**< pressure variable */
   SCIP_Bool             needPIvar;                    /**< whether or not PIvar is needed */
   SCIP_VAR*             PIvar;                        /**< squared pressure variable */
   SCIP_VAR*             pSlackVar;                    /**< slack variable for minimizing relaxation of pressure bounds */
   SCIP_VAR*             nodeMixingRatio;              /**< node mixing ratio variable in Mass % needed to calc and link the arc mixing ratios in Mass % */
   SCIP_VAR*             nodeMixingRatioMol;           /**< node mixing ratio variable in Mol % */
   SCIP_Real             specificHeatCap;              /**< specific heat capacity ratio of the entering gas if node is an entry */
   SCIP_Real             molarMass;                    /**< molar mass of entering gas if node is an entry */
};
typedef struct GAS_node GAS_Node;

/** types of arcs */
enum arc_type
{
   UNKNOWNARCTYPE,                                     /**< unknown type (for debugging) */
   PIPE,                                               /**< pipe */
   VALVE,                                              /**< valve */
   CONTROLVALVE,                                       /**< control valve */
   CS,                                                 /**< compressor station */
   RESISTOR,                                           /**< resistor */
   SHORTPIPE                                           /**< short pipe */
};
typedef enum arc_type Arc_type;

/** arc of gas network */
struct GAS_arc
{
   GAS_Node*             sourcenode;                   /**< source of pipe */
   GAS_Node*             targetnode;                   /**< target of pipe */
   struct GAS_arc*       next_inarc;                   /**< next pipe in list with the same target or NULL */
   struct GAS_arc*       next_outarc;                  /**< next pipe in list with the same source or NULL */
   char*                 id;                           /**< id of pipe */
   SCIP_Real             flowMin;                      /**< minimal flow */
   SCIP_Real             flowMax;                      /**< maximal flow */
   Arc_type              type;                         /**< type of arc; e.g., pipe, valve, ... */
   int                   arcposition;                  /**< position in arcs_ptr */
   void*                 detailed_info;                /**< pointer to typed data */
   SCIP_VAR*             flowvar;                      /**< flow variable on arc */
   SCIP_VAR*             mixingRatio;                  /**< mixing ratio on the current arc in Mol % */
   SCIP_VAR*             mixingRatioMass;              /**< mixing ratio on the current arc in Mass % */
   SCIP_VAR*             Lambda;                       /**< Resistance variable for Weymouth equation depending on speed of sound */
   SCIP_VAR*             positiveFlowBinvar;           /**< binary variable that indicates positive flow */
   SCIP_VAR*             negativeFlowBinvar;           /**< binary variable that indicates negative flow */
};
typedef struct GAS_arc GAS_Arc;

/** data structure for a pipe */
struct GAS_pipe
{
   SCIP_Real             length;                       /**< length of pipe */
   SCIP_Real             diameter;                     /**< diameter of pipe */
   SCIP_Real             roughness;                    /**< roughness of pipe */
   SCIP_Real             pressureMax;                  /**< maximum pressure in pipe */
   SCIP_Real             slope;                        /**< slope of pipe */
   SCIP_VAR*             PIdiffVar;                    /**< variable Delta_P for pressure difference p_u^2 p_v^2 in algebraic model */
};
typedef struct GAS_pipe GAS_Pipe;

/** data structure for a valve */
struct GAS_valve
{
   SCIP_Real             pressureDifferentialMax;      /**< maximal pressure difference of a valve */
   SCIP_VAR*             valve_binvar;                 /**< binary variable for the mode of a valve (open or closed) */
   GAS_Arc*              bypassed_cv;                  /**< pointer to an controlvalve arc, if the valve is the bypass of an cv; else NULL */
   GAS_Arc*              bypassed_cs;                  /**< pointer to an controlvalve arc, if the valve is the bypass of an cv; else NULL */
};
typedef struct GAS_valve GAS_Valve;

/** data structure for control valve */
struct GAS_controlvalve
{
   SCIP_Real             pressureLossIn;               /**< pressure loss at controlvalve entry */
   SCIP_Real             pressureLossOut;              /**< pressure loss at controlvalve exit */
   SCIP_Real             dragFactorIn;                 /**< drag factor of a nonlinear resistor at the entry */
   SCIP_Real             dragFactorOut;                /**< drag factor of a nonlinear resistor at the exit */
   SCIP_Real             pressureInMin;                /**< minimal pressure at controlvalve entry */
   SCIP_Real             pressureOutMax;               /**< maximal pressure at controlvalve exit */
   SCIP_Real             pressureDifferentialMin;      /**< minimal pressure difference of controlvalve */
   SCIP_Real             pressureDifferentialMax;      /**< maximal pressure difference of controlvalve */
   SCIP_VAR*             binvar;                       /**< binary variable for active/inactive cv */
   int                   internalBypassRequired;       /**< 0/1 if we need a bypass */
   GAS_Arc*              bypass;                       /**< pointer to the gas_arc of the bypass */
};
typedef struct GAS_controlvalve GAS_Controlvalve;

/** types of resistor */
enum resistor_type {
   LINEAR,               /**< linear resistor behavior */
   NONLINEAR             /**< nonlinear resistor behavior */
};
typedef enum resistor_type Resistor_type;

/** data structure for resistor */
struct GAS_resistor
{
   SCIP_Real             dragFactor;                   /**< dragfactor of the resistor */
   SCIP_Real             diameter;                     /**< diameter of the resistor */
   SCIP_Real             pressureLoss;                 /**< pressure loss parameter */
   Resistor_type         type;                         /**< type of resistor */
   int                   LR_position;                  /**< position of the variables in the LR-Variable arrays; -1 if NLR */
   int                   NLR_position;                 /**< position of the variables in the NLR-Variable arrays; -1 if LR */
};
typedef struct GAS_resistor GAS_Resistor;

/** data structure for a configuration of a compressor station with data from the cs file*/
struct GAS_csboxcons
{
   char*                 id;                           /**< name of the config */
   SCIP_Real             massFlowMin;                  /**< minimal mass flow from cs file */
   SCIP_Real             massFlowMax;                  /**< maximal mass flow from cs file */
   SCIP_Real             pressureInMin;                /**< minimal pressure at entry */
   SCIP_Real             pressureInMax;                /**< maximal pressure at entry */
   SCIP_Real             pressureOutMin;               /**< minimal pressure at exit */
   SCIP_Real             pressureOutMax;               /**< maximal pressure at exit */
   SCIP_Real             pressureIncAbsMin;            /**< minimal bound on absolute pressure increase */
   SCIP_Real             pressureIncAbsMax;            /**< maximal bound on absolute pressure increase */
   SCIP_Real             pressureIncRelMin;            /**< minimal bound on relative pressure increase */
   SCIP_Real             pressureIncRelMax;            /**< maximal bound on relative pressure increase */
   SCIP_Real**           facetcoeff;                   /**< matrix with coefficents for facets*/
   SCIP_Real*            rhs;                          /**< right hand side of facet*/
   char*                 a;                            /**< information about variable */
   char*                 a_unit;                       /**< information about variable */
   char*                 b;                            /**< information about variable */
   char*                 b_unit;                       /**< information about variable */
   char*                 c;                            /**< information about variable */
   char*                 c_unit;                       /**< information about variable */
   int                   numfacets;                    /**< number of facets */
   SCIP_VAR*             config_binvar;                /**< binary variable for config */
   SCIP_VAR*             BCM_flow;                     /**< flow on the compressor arc */
   SCIP_VAR*             BCM_pin;                      /**< pressure before the compressor */
   SCIP_VAR*             BCM_pout;                     /**< pressure after the compressor */
};
typedef struct GAS_csboxcons GAS_CSBoxCons;

/** data structure for piston compressor */
struct GAS_pistoncs
{
   char*                 id;                           /**< compressor name */
   SCIP_Real             speedMin;                     /**< minimal speed of cs */
   SCIP_Real             speedMax;                     /**< maximal speed of cs */
   /* etc. */
};
typedef struct GAS_pistoncs GAS_PistonCS;

/** data structure for turbo compressor */
struct GAS_turbocs
{
   char*                 id;                           /**< compressor name */
   SCIP_Real             speedMin;                     /**< minimal speed of cs */
   SCIP_Real             speedMax;                     /**< maximal speed of cs */
   SCIP_Real             speedMatrix[9];               /**< entries of matrix for bi-quadratic ansatz */
   SCIP_Real             surgeline[3];                 /**< coefficients of surgeline */
   SCIP_Real             chokeline[3];                 /**< coefficients of chokeline */
   /* etc. */
};
typedef struct GAS_turbocs GAS_TurboCS;

/* data structure for a compressor station configuration */
struct GAS_csconfig
{
   char*                 id;                           /**< name of the configuration */
   int                   numSerialStages;              /**< number of cs units connected in serial */
   int                   numCsStageOne;                /**< number of compressors used in parallel in stage one */
   char*                 stageOneCsOne;                /**< first compressor in stage one */
   char*                 stageOneCsTwo;                /**< second compressor in stage one */
};
typedef struct GAS_csconfig GAS_CSConfig;

/** data structure for compressor station */
struct GAS_cs
{
   SCIP_Real             pressureLossIn;               /**< pressure loss at compressor entry */
   SCIP_Real             pressureLossOut;              /**< pressure loss at compressor exit */
   SCIP_Real             dragFactorIn;                 /**< drag factor of a nonlinear resistor at the entry */
   SCIP_Real             dragFactorOut;                /**< drag factor of a nonlinear resistor at the exit */
   SCIP_Real             diameterIn;                   /**< diameter of a nonlinear resistor at the entry */
   SCIP_Real             diameterOut;                  /**< diameter of a nonlinear resistor at the exit */
   SCIP_Real             pressureInMin;                /**< minimal pressure at entry */
   SCIP_Real             pressureOutMax;               /**< maximal pressure at exit */
   int                   numPistonCS;                  /**< number of piston compressors */
   GAS_PistonCS*         pistoncs;                     /**< array of piston compressor structs */
   int                   numTurboCS;                   /**< number of turbo compressors */
   GAS_TurboCS*          turbocs;                      /**< array of turbo compressor structs */
   int                   numconfigurations;            /**< number of configurations */
   GAS_CSConfig*         configurations;               /**< pointer to configuration struct(s) */
   GAS_CSBoxCons*        boxcons;                      /**< pointer to the configuration struct */
   SCIP_VAR*             compressor_binvar;            /**< binary variable for compressor station only used if there are no configurations */
   SCIP_VAR*             NLRin_pressure;               /**< pressure variable for the nlr at entry */
   SCIP_VAR*             NLRin_pressureDiff;           /**< pressure difference for nlr at entry */
   SCIP_VAR*             NLRout_pressure;              /**< pressure variable for the nlr at entry */
   SCIP_VAR*             NLRout_pressureDiff;          /**< pressure difference for nlr at entry */
   int                   internalBypassRequired;       /**< information about existing bypass: 0 or 1 */
   int                   bypassPosition;               /**< position in valve variables or -1 */
   GAS_Arc*              bypass;                       /**< pointer to the gas_arc of the bypass */
   SCIP_VAR*             mixingRatio;                  /**< mixing ratio on the current arc in Mol % */
   SCIP_VAR*             mixingRatioMass;              /**< mixing ratio on the current arc in Mass % */
};
typedef struct GAS_cs GAS_CS;

/** data for the spanning tree of the network
 *  @note: We assume that the network is still connected without the compressors
 */
struct GAS_spanningtree
{
   GAS_Node*             root_node;                    /**< root node of the tree */
   SCIP_Bool             foundCycles;                  /**< false if network is a tree, or not connected */
   int**                 wayToNode;                    /**< numarcs x numnodes 0/1 matrix, entries correspond to path from root to node without orientation */
   int*                  arcsInTree;                   /**< array with entries 0,1; 1 if arc is part of the tree, else 0 */
};
typedef struct GAS_spanningtree GAS_SpanningTree;

/** struct for cycles found in the graph and generating noCircularFlowCuts */
struct GAS_cycle
{
   int*                  combinedBasisCycles;          /**< 0/1 array with length numBasisCycles */
   int*                  arcsInCycle;                  /**< 0/+-1 array with the arcs in the cycle */
   int                   nr_cycle;                     /**< number of the cycle */
   struct GAS_cycle*     next_cycle;                   /**< pointer to the next cycle found */
   struct GAS_cycle*     previous_cycle;               /**< pointer to the previous cycle found */
};
typedef struct GAS_cycle GAS_Cycle;

/** gas network */
struct GAS_network
{
   int                   numnodes;                     /**< number of nodes */
   int                   numarcs;                      /**< number of arcs */
   int                   maxnumincidentarcs;           /**< maximal number of arcs incident to one node */
   int                   numflowdirvars;               /**< number of arcs with a flow direction binvars */
   int                   numpipes;                     /**< number of pipes */
   int                   numshortpipes;                /**< number of shortpipes */
   int                   numvalves;                    /**< number of valves */
   int                   numcontrolvalves;             /**< number of controlvalves */
   int                   numcompressor;                /**< number of compressor stations */
   int                   numcsconfigurations;          /**< number of configurations for the compressor stations */
   int                   numcsvarsfornlrs;             /**< number of variables for nonlinear resistors in compressor stations */
   int                   numlinresistor;               /**< number of linear resistors */
   int                   numnonlinresistor;            /**< number of nonlinear resistors */
   int                   numpivars;                    /**< number of variables Pivars */
   int                   numBasisCycles;               /**< number of fundamental cycles found */
   GAS_Node*             nodes_ptr;                    /**< pointer to nodes of graph */
   GAS_Arc*              arcs_ptr;                     /**< pointer to arcs of graph */
   GAS_SpanningTree*     spanningTree;                 /**< pointer to the spanning tree */
   GAS_Cycle*            firstCycle;                   /**< pointer to the first cycle in the network */
   GAS_Cycle*            lastCycle;                    /**< pointer to the last cycle found */
   GAS_Cycle*            firstCombinedCycle;           /**< first cycle that is not a fundamental cycle */
   SCIP_HASHTABLE*       Nodeshashtable;               /**< pointer to a hashtable for nodes */
   SCIP_HASHTABLE*       Arcshashtable;                /**< pointer to a hashtable for arcs */
   SCIP_Real             gasTemperature;               /**< in Kelvin */
   SCIP_Real             pseudocriticalPressure;       /**< in bar */
   SCIP_Real             pseudocriticalTemperature;    /**< in Kelvin */
   SCIP_Real             normDensity;                  /**< in kg/(m^3) */
   SCIP_Real             molarMass1;                   /**< molar mass of gas 1 */
   SCIP_Real             molarMass2;                   /**< molar mass of gas 2 */
   SCIP_Real             specificHeatCap1;             /**< specific heat capacity ratio of gas 1 */
   SCIP_Real             specificHeatCap2;             /**< specific heat capacity ratio of gas 2 */
};
typedef struct GAS_network GAS_Network;

/** problem data */
struct SCIP_ProbData
{
   GAS_Network*          network;                      /**< gas network */
   char*                 netname;                      /**< name of the network */
   char*                 scnname;                      /**< name of the scn file */
   SCIP_VAR*             objective;                    /**< objective value */
   SCIP_VAR**            flowvars;                     /**< flow variables for all arcs */
   SCIP_VAR**            Lambda;                       /**< variable used for speed of sound mixing on all arcs */
   SCIP_VAR**            mixingRatio;                  /**< variable which describes the ratio in Mol % of a certain gas (natural Gas) on the arc in [0,1] */
   SCIP_VAR**            mixingRatioMass;              /**< variable which describes the ratio in Mass % of a certain gas (natural Gas) on the arc in [0,1] */
   SCIP_VAR**            nodeMixingRatio;              /**< mixing ratio at the node in Mass % */
   SCIP_VAR**            nodeMixingRatioMol;           /**< mixing ratio at the node in Mol % */
   SCIP_VAR**            netflowvars;                  /**< net flow in entries */
   SCIP_VAR**            positiveFlowBinvars;          /**< binary variables that indicate positive flow */
   SCIP_VAR**            negativeFlowBinvars;          /**< binary variables that indicate negative flow */
   SCIP_VAR**            pressurevars;                 /**< pressure variables for all nodes */
   SCIP_VAR**            PIvars;                       /**< pressure variables squared */
   SCIP_VAR**            pSlackVars;                   /**< slack variables for minimizing bound relaxation */
   SCIP_VAR**            PIdiffvars;                   /**< variable used for pressure constraint */
   SCIP_VAR**            VALVE_binvars;                /**< binary variables for valve */
   SCIP_VAR**            CV_binvars;                   /**< binary variables for controlvalve */
   SCIP_VAR**            NLR_DeltaVars;                /**< continuous variables used in resistor constraint */
   SCIP_VAR**            NLR_AbsDeltaVars;             /**< continuous variables used in resistor constraint */
   SCIP_VAR**            NLR_AbsFlowVars;              /**< continuous variables used in resistor constraint */
   SCIP_VAR**            LR_posFlowDir;                /**< binary variables used in linear resistor constraint */
   SCIP_VAR**            LR_smoothingFlow;             /**< continuous variables used in linear resistor constraint */
   SCIP_VAR**            LR_negFlowDir;                /**< binary variables used in linear resistor constraint */
   SCIP_VAR**            CS_binvars;                   /**< binary variables used in compressor constraint */
   SCIP_VAR**            CS_resistorvars;              /**< variables for nonlinear resistors in compressor stations */
   SCIP_VAR**            BCM_flowCS;                   /**< continuous flow variables for compressor station with box constraint model */
   SCIP_VAR**            BCM_pressureIn;               /**< continuous pressure variables for incoming pressure at compressor stations with box constraint model */
   SCIP_VAR**            BCM_pressureOut;              /**< continuous pressure variables for outgoing pressure at compressor stations with box constraint model */
   SCIP_Bool             noFlowBinvars;                /**< if no binary variables for flowdirections should be used */
   SCIP_Bool             binaryFlowCons;               /**< if binary flow conservation constraints should be used */
   SCIP_Bool             noCycles;                     /**< false if noCircularFlow constraints should be used  */
   SCIP_Bool             allCycles;                    /**< if all cycles in the graph should be found */
   SCIP_Bool             algebraic;                    /**< use algebraic model */
   SCIP_Bool             mixing;                       /**< use mixing model */
   SCIP_Bool             flowConsMixing;               /**< use flow cons mixing model */
   SCIP_Bool             flowConsMixingNode;               /**< use node flow cons mixing model */
   SCIP_Bool             linearMixing;                 /**< use linear mixing model mass% convex comb of speed of sounds */
   SCIP_Bool             nodeMu;                       /**< use the node mu mixing model*/
   SCIP_Bool             mol;                          /**< directly calculate Mol % without use of Mass % */
   SCIP_Bool             binVarEQone;                  /**< set z+ +z- == 1. This will force flow over every arc. */
   SCIP_Bool             papay;                        /**< if formula of papay is used */
   SCIP_Bool             boxConstraintModel;           /**< use box constraint model */
   SCIP_Bool             additionalFacets;             /**< use additional facets box constraint model */
   SCIP_Bool             charDiagram;                  /**< use characteristic diagrams for compressors */
   SCIP_Bool             relaxLowerBounds;             /**< if true the lower bounds of the pressure variables should be relaxed */
   SCIP_Bool             relaxUpperBounds;             /**< if true the upper bounds of the pressure variables should be relaxed by relaxUBvalue */
   SCIP_Bool             relaxCVbounds;                /**< if pressure bounds of controlvalves should be relaxed by relaxCVvalue */
   SCIP_Real             relaxLBvalue;                 /**< maximum value the lower bound can be violated */
   SCIP_Real             relaxUBvalue;                 /**< maximum value the upper bound can be violated  */
   SCIP_Real             relaxCVvalue;                 /**< relaxation of controlvalve pressure bounds */
   SCIP_Bool             resistorCVInternal;           /**< whether resistors of a controlvalve are treated to be internal */
   SCIP_Bool             minCompSum;                   /**< if true the objective function minimizes the sum of the opened compressors */
   SCIP_Bool             minCompInc;                   /**< if true the objective function minimizes the pressure increase of the opened compressors */
   SCIP_Bool             minPressSum;                  /**< if true the objective function minimizes the sum of the pressure values */
   SCIP_Bool             noObjective;                  /**< if true the problem is a pure feasibility problem */
   SCIP_Bool             powerLoss;                    /**< if true the objective function minimizes the power loss */
   SCIP_Bool             minSlack;                     /**< if true, minimize maximal difference to original variable bounds */
   SCIP_Bool             minSlackPerBound;             /**< minimize sum of slack variables for all relaxed pressure bounds */
   SCIP_Real             idealCSminIncrease;           /**< minimal pressure increase by an ideal compressor */
   SCIP_Real             idealCSmaxIncrease;           /**< maximal pressure increase by an ideal compressor */
   SCIP_Bool             csusemaxincrease;             /**< whether a maximal increase for compressors should be used */
   SCIP_Bool             maxFlow;                      /**< if true the objective function maximizes the flow */
   SCIP_Bool             addComponentCuts;             /**< whether component cuts should be added */
   SCIP_Bool             minLambda;                    /**< if the overall resistance in the pipes should be minimized (only for mixing) */
   SCIP_Bool             approx;                       /**< approximate mixing */
   SCIP_Bool             addParallelCSCuts;            /**< whether cuts for equal parallel compressors should be added */
   SCIP_Bool             reduceFlowConservation;       /**< skip on flow conservation for one node, to make matrix full rank */
   SCIP_Real             molarMass1;                   /**< molar mass of gas 1 */
   SCIP_Real             molarMass2;                   /**< molar mass of gas 2 */
   SCIP_Real             specificHeatCap1;             /**< specific heat capacity ratio of gas 1 */
   SCIP_Real             specificHeatCap2;             /**< specific heat capacity ratio of gas 2 */
   SCIP_Bool             aggregateParallelFlowdirVars; /**< whether flow direction variables of parallel arcs should be aggregated */
   SCIP_Bool             preprocessPassiveComponents;  /**< preprocess passive components to fix flows */
   SCIP_Real             scalescn;                     /**< factor to scale scenario flow values with */
};
