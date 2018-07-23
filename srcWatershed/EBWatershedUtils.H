#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBWATERSHEDUTILS_H_
#define _EBWATERSHEDUTILS_H_

#include "EBIndexSpace.H"
#include "BaseIF.H"

#include "NamespaceHeader.H"

/**
   Agglomeration of the many parameters of an AMR calc
*/
class WatershedParameters
{
public:
    ///almost everything gets a default
    WatershedParameters();

    ~WatershedParameters();

    //Initial values as string for muParser
    std::string m_initial_pressurehead;

    //switched for on/off overland flow solver
    bool m_includeSurfaceSolver;

    //switched for on/off overland flow solver
    bool m_includeSUBSurfaceSolver;
    
    //source term for surface solver
    std::string m_scalarSurfaceSourceFunc;

    //use cc or centroid data for irregular cells
    bool m_irregular_value_cc;

    //Domain parameters
    RealVect              m_domainLength;// length of the domain in the 0 direction
    int                   m_maxLevel; // number of  highest AMR Level (0 - numLevels -1)
    int                   m_maxBoxSize;   //  1-D maximum size of boxes in grid hierarchy

    //time integrator parameters
    Real                  m_dtGrowFactor;//maximum growth factor for time step
    Real                  m_maxDt;//maximum time step
    bool                  m_useSmallerDt;// use smallerDt if not converged
    // bool                  m_useDivergedSolutionRestart;  //try use diverged solution as approximation when restart with smaller step (experimental)
 
    int                   m_checkpointInterval; // number of steps between checkpoints (-1 means no checkpoints)
    int                   m_plotInterval;   // number of steps between plotfiles (-1 means no plotfiles)
    int                   m_orderTimeIntegration;  //time intergrator type

    bool                m_printOutputRate; //print surface flow output rate
    //grid parameters
    Real                  m_fillRatio; // minimum percentage of tagged cells per box ala B-R
    int                   m_blockFactor; // blocking factor ala B-R
    int                   m_regridInterval; //time steps between regrids
    Real                  m_refineThreshold; // factor on which to refine
    int                   m_nestingRadius; //???
    int                   m_tagBuffer; //???
    bool                  m_refineAllIrreg; //refine all irregular cells
    bool                  m_refineUnderresolved;//refine underresolved (if m_refineAllIrreg=false)

    int                   m_tagShrinkDomain;  //number of cells to shrink domain by so that there is no refinement at side walls
    Vector<int>           m_refRatio;   // refinement ratio between levels.  0 is the refrat between 0 and 1.

    int                   m_verbosity;  // verbosity of amr class
    
    BaseIF*               m_implicitBaseIF; //implicit function from GeometryShop

    //function to tag cells for refinment
    std::string m_tagCellsFunc;

private:
    bool isInitialized;
};

//get stuff from input file
extern void getWatershedParameters(WatershedParameters&  a_params,
				   ProblemDomain&   a_coarsestInputDomain);

///define ebindexspace from parameters (and parmparse)
extern void defineWatershedGeometry(WatershedParameters&  a_params,
				    const ProblemDomain&   a_coarsestInputDomain, 
				    int a_whichGeomForced = -1,
				    EBIndexSpace* a_ebisPtr = Chombo_EBIS::instance());

#include "NamespaceFooter.H"
#endif
