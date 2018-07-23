#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif


#include <iomanip>
#include <cmath>
#include <cstdio>
#include "ParmParse.H"
#include "REAL.H"
#include "memusage.H"
#include "memtrack.H"
#include "WatershedIBC.H"
#include "EBAMRWatershed.H"
#include "EBAMRWatershedSolver.H"
#include "EBSUBSurfaceOps.H"
#include "MeshRefine.H"
#include "BRMeshRefine.H"
#include "BCFunc.H"
//#include "EBPhysIBCFactory.H"
#include "EBArith.H"
#include "EBPWLFineInterp.H"
#include "EBFluxFactory.H"
#include "EBCellFactory.H"
#include "EBAMRIO.H"
#include "EBAMRDataOps.H"
#include "EBLevelDataOps.H"
#include "EBEllipticLoadBalance.H"

#include "UsingNamespace.H"

static const Real Tol = 1e-15;

/**********************/
EBAMRWatershed::EBAMRWatershed(const WatershedParameters&      a_params,
			       const WatershedIBCFactory&       a_ibcfact,
			       const ProblemDomain&      a_coarsestDomain,
			       const EBIndexSpace* const a_ebisPtr)
    :m_ebisPtr(a_ebisPtr),
    m_richardsCompGrid(NULL)
{
    CH_TIME("EBAMRWatershed::EBAMRWatershed");
    if (a_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::EBAMRWatershed" << endl;
    }

    //set parameters of the run
    m_params    = a_params;

    //create initial and boundary condition object
    m_ibc    =   a_ibcfact.create(m_params.m_includeSurfaceSolver);

    //setting ghostcells
    m_numGhostCells = 2;

    //resize vectors and set them where we can
    RealVect coarsestDx;
    for(int idir=0; idir<3;idir++)
    {
	coarsestDx[idir] = m_params.m_domainLength[idir]/Real(a_coarsestDomain.size(idir));
    }
    int nlevels = m_params.m_maxLevel + 1;
    m_domain.resize(nlevels); 

    m_vectDx.resize(nlevels);
    m_grids.resize(nlevels);
    m_ebisl.resize(nlevels);
    m_eblg.resize(nlevels);
    m_psi_current.resize(nlevels, NULL);
    
    //allocate CompGrid richards equation solver
    m_richardsCompGrid  = new  EBPetscCompGridRichards(nlevels-1, m_params.m_includeSurfaceSolver);
  
    if (m_params.m_includeSurfaceSolver)
    {
	    m_EBSource.resize(nlevels);
    }
    //allocate Petsc Watershed solver
    m_watershedSolver = new EBAMRWatershedSolver(m_richardsCompGrid, m_ibc);

    allocateDataHolders();

    m_domain[0] = a_coarsestDomain;
    if (m_params.m_includeSurfaceSolver)
    { 
	// m_domain_surface.resize(nlevels); 
	// Box surfaceBox = m_domain[0].domainBox();
	// surfaceBox.setBig(2, (surfaceBox.smallEnd(2) + 1 ));
	// m_domain_surface[0].define(surfaceBox);
    }
    m_vectDx[0]     =   coarsestDx;
    for (int ilev = 1; ilev < nlevels; ilev++)
    {
	CH_assert(m_params.m_refRatio[ilev-1] > 0);
	m_domain[ilev] = refine(m_domain[ilev-1], m_params.m_refRatio[ilev-1]);
	m_vectDx[ilev] = m_vectDx[ilev-1]/Real(m_params.m_refRatio[ilev-1]);
	if (m_params.m_includeSurfaceSolver)
	{
	    // IntVect refRatio = IntVect::Unit* m_params.m_refRatio[ilev-1];
	    // refRatio[2] = 1;
	    //m_domain_surface[ilev] = refine(m_domain_surface[ilev-1], refRatio);
	}
    }
  
    m_time = 0.0;
    m_curStep = 0;
    m_dt = -1.0;
    //setup still needs to get called
    m_isSetup  = false;

    m_pointsUpdated = 0;
}
/**********/
void
EBAMRWatershed::allocateDataHolders()
{
    CH_TIME("EBAMRWatershed::allocateDataHolders");
    for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
	m_psi_current[ilev]     = new LevelData<EBCellFAB>();
	if (m_params.m_includeSurfaceSolver)
	{
	    m_EBSource[ilev] = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >());
	}
    }
}
/**********/
EBAMRWatershed::~EBAMRWatershed()
{
    CH_TIME("EBAMRWatershed::~EBAMRWatershed");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::~EBAMRWatershed" << endl;
    }
    for (int ilev = 0; ilev <= m_params.m_maxLevel; ilev++)
    {
	delete m_psi_current[ilev];
	if (m_params.m_includeSurfaceSolver) 
	{

	}
    }
}
/**********/
void
EBAMRWatershed::setupForAMRRun()
{
    CH_TIME("EBAMRWatershed::setupForAMRMRun");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::setupForAMRRun" << endl;
    }
    m_isSetup= true;
    m_doRestart    = false;

    if (m_params.m_verbosity > 3 && m_params.m_includeSurfaceSolver)
    {
	pout()<<"EBAMRWatershed::setupForAMR: Overland surface flow solver is included"<<endl;
    }

    //we're generating the hiearachy dyamically
    //modelled on AMR::initialGrid()
    Vector<Vector<Box> > old_boxes(1);
    Vector<Vector<Box> > new_boxes;
  
    // //also keep old tags around
    // Vector<IntVectSet> oldTags;
  
    //define base mesh
    //chop up base level into grids to satisfy box size requirements
    domainSplit(m_domain[0], old_boxes[0], m_params.m_maxBoxSize,
		m_params.m_blockFactor);
  
  
    //now intialize data for existing hierarchy
    initialGrid(old_boxes);
      
    initialData();

    m_finestLevel = 0;
      
    //now generate more levels if necessary
    int top_level = 0;
  
    bool moreLevels = (m_params.m_maxLevel > 0);

  
    //create grid generation object
    BRMeshRefine meshrefine;
  
    if (moreLevels)
    {
	meshrefine.define(m_domain[0],              m_params.m_refRatio,
			  m_params.m_fillRatio,     m_params.m_blockFactor,
			  m_params.m_nestingRadius, m_params.m_maxBoxSize);
    }
  
    while (moreLevels)
    {
	//default is moreLevels = false
	//(only repeat loop in the case where a new level
	//is generated which is still coarser than maxLevel)
	moreLevels = false;
      
	int base_level = 0;
	int old_top_level = top_level;
      
	Vector<IntVectSet> tagsVect(top_level+1);
	tagCells(tagsVect);
      
	int new_finest = meshrefine.regrid(new_boxes, tagsVect,
					   base_level, top_level,
					   old_boxes);
	if (new_finest > top_level) top_level++;
      
	old_boxes = new_boxes;
      
	//now see if we need another pass through grid generation
	if ((top_level<m_params.m_maxLevel) && (top_level > old_top_level))
	    moreLevels = true;
      
      
	//if we added another level, reinintialize everything again
	if (top_level > old_top_level)
	{
	    initialGrid(new_boxes);
	  
	    initialData();
	}
    } //end loop over regridding passes
      
    defineIrregularData();
    //finally, call post-initialization
    postInitialize();
}
/**************************/
Real
EBAMRWatershed::run(Real a_maxTime, 
		    int a_maxStep)
{
    CH_TIME("EBAMRWatershed::run");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::run" << endl;
    }
    CH_assert(m_isSetup);

    //  ParmParse pp;

    if (m_params.m_verbosity > 0)
    {
	pout() << "EBAMRWatershed: starting run " << endl;
    }

    //skip computeDt if in the first step since it was already calculated above
    bool skipDtComputation = true;
   
    Real initDt = m_dt;
#ifdef CH_USE_HDF5
    if (m_params.m_plotInterval > 0)
    {
	if (m_params.m_verbosity > 0)
	{
	    pout() << "EBAMRWatershed: writing plot file" << endl;
	}
	writePlotFile();
    }
#endif

    pout() << "EBAMRWatershed: starting time advance  " << endl;
   
    //advance solution until done
    while ((a_maxTime > m_time) && (m_curStep < a_maxStep) )
    {
	//advance step number
	++m_curStep;
      
	pout() << "##############################################################################" << endl;
	pout() << "EBAMRWatershed: step " << m_curStep << endl;

	//do regridding if appropriate
	if ( (m_curStep%m_params.m_regridInterval == 0) &&
	     (m_params.m_regridInterval > 0) &&
	     (m_params.m_maxLevel > 0) )
        {
	    regrid();
        }

	//compute new dt
	if (skipDtComputation)//do not skip computation in future iterations
        {
	    skipDtComputation = false;
        }
	else
        {
	    m_dt = computeDt();
        }

	if (m_dt < 1e-5 * initDt)
	{
	    pout() << "EBAMRWatershed::run -- Time step too small" << endl;
	    break;
	}

	if ((m_time + m_dt) > a_maxTime)
	{
	    m_dt = a_maxTime - m_time;
	}

	pout() << "Beginning of time step " << m_curStep
	       << ", start time = " << m_time
	       << ", dt = " << m_dt << endl;

	//do timestep
	advance();

	postTimeStep();

	//advance time
	m_time += m_dt;

	//dump plotfile and checkpointfile before regridding
#ifdef CH_USE_HDF5
	if ((m_curStep%m_params.m_plotInterval == 0) && (m_params.m_plotInterval > 0))
        {
	    if (m_params.m_verbosity > 0)
            {
		pout() << "EBAMRWatershed: writing plot file" << endl;
            }
	    writePlotFile();
        }

	if ((m_curStep%m_params.m_checkpointInterval == 0) && (m_params.m_checkpointInterval > 0))
        {
	    if (m_params.m_verbosity > 0)
            {
		pout() << "EBAMRWatershed: writing checkpoint file" << endl;
            }
	    writeCheckpointFile();
        }
#endif

	pout() << "End of time step " << m_curStep << ", end time = " << m_time << endl;

    }//end while loop over timesteps
   
    pout() << "##############################################################################" << endl;
    pout() << "total number of points updated = " << m_pointsUpdated << endl;
    pout() << "" << endl;
   
    conclude();
   
    return m_time;
}
/******************/
void
EBAMRWatershed::conclude()
{
    CH_TIME("EBAMRWatershed::conclude");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::conclude" << endl;
    }
#ifdef CH_USE_HDF5
    if (m_params.m_plotInterval >= 0)
    {
	writePlotFile();
    }
  
    if (m_params.m_checkpointInterval >= 0)
    {
	writeCheckpointFile();
    }
#endif
}
/*****************/
void
EBAMRWatershed::tagCells(Vector<IntVectSet>& a_tags)
{
    CH_TIME("EBAMRWatershed::tagCells");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::tagCells" << endl;
    }
    int numlevels = a_tags.size();
    if (numlevels > m_finestLevel+1) numlevels = m_finestLevel+1;
  
    for (int lev=0; lev<numlevels; lev++)
    {
	IntVectSet& levelTags = a_tags[lev];
      
	tagCellsLevel(levelTags, lev);
#ifdef CH_MPI
	int gatherint = 0;
	if (!levelTags.isEmpty()) gatherint = 1;
	int itags;
	MPI_Allreduce(&gatherint, &itags, 1, MPI_INT,
		      MPI_MAX, Chombo_MPI::comm);
	bool thereAreTags = (itags==1);
	if (!thereAreTags)
        {
	    MayDay::Error("EBAMRWatershed::tagCells -- numlevels > 1 and no cells tagged");
        }
#endif
    }
}
/*****************/
void
EBAMRWatershed::regrid()
{
    CH_TIME("EBAMRWatershed::regrid");
    //smoothing to make regridding less heinous
    preRegrid();

    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::regrid" << endl;
    }
    //don't regrid base grid
    int lbase = 0;

    if (m_params.m_maxLevel > 0)
    {
	//first, construct tags
	int top_level = Min(m_finestLevel, m_params.m_maxLevel-1);
	Vector<IntVectSet> tagsVect(top_level+1);
	Vector<Vector<Box> > new_grids;
	Vector<Vector<Box> > vectBoxes(top_level+1);
	for (int ilev = 0; ilev <= top_level; ilev++)
        {
	    domainSplit(m_domain[ilev], vectBoxes[ilev], m_params.m_maxBoxSize,
			m_params.m_blockFactor);
        }
      
	// do tagging
	tagCells(tagsVect);

	int new_finest_level;
      
	BRMeshRefine meshrefine(m_domain[0], m_params.m_refRatio,
				m_params.m_fillRatio, m_params.m_blockFactor,
				m_params.m_nestingRadius, m_params.m_maxBoxSize);
      
	new_finest_level = meshrefine.regrid(new_grids,
					     tagsVect,
					     lbase,
					     top_level,
					     vectBoxes);
      

      
	//can only add one level at a time
	new_finest_level = Min(m_finestLevel+1, new_finest_level);

	if ((m_finestLevel != new_finest_level) && (m_params.m_verbosity >= 2))
        {
	    pout() << "finest level changes here from " << m_finestLevel << " to "<< new_finest_level << endl;
        }

	//allow for levels to change
	m_finestLevel = Min(new_finest_level, m_params.m_maxLevel);

      
	//now redefine grid hierarchy
	regrid(new_grids);
      
	defineIrregularData();
      
	//finish up
	postRegrid();
    } //end if max level > 0
}
/*********************/
void
EBAMRWatershed::preRegrid()
{
    CH_TIME("EBAMRWatershed::preRegrid");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::preRegrid" << endl;
    }
}
/*********************/
void
EBAMRWatershed::postRegrid()
{
    CH_TIME("EBAMRWatershed::postRegrid");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::postRegrid" << endl;
    }
    m_watershedSolver->deallocateSolver();
}

/*********************/
void
EBAMRWatershed::defineGrids(const Vector<Vector<Box> >& a_vectBoxes)
{
    CH_TIME("EBAMRWatershed::defineGrids");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::defineGrids" << endl;
    }
    m_finestLevel = 0;
    //now define all of the storage we need
    int start = 0;
    int end = a_vectBoxes.size()-1;

    for (int ilev = start; ilev<=end; ilev++)
    {
	if (a_vectBoxes[ilev].size() > 0)
	{
	    m_finestLevel = ilev;
	    //first do load balance
	    Vector<int> procAssign;
	    mortonOrdering((Vector<Box>&)(a_vectBoxes[ilev]));
	    LoadBalance(procAssign,  a_vectBoxes[ilev]);//, m_domain[ilev], false, m_ebisPtr );
	    //Vector<Box>& castBoxes = const_cast<Vector<Box> & >(a_vectBoxes[ilev]);
	    //EBEllipticLoadBalance(procAssign,  castBoxes, m_domain[ilev], false, m_ebisPtr );
	    // 	  DisjointBoxLayout newDBL(a_vectBoxes[ilev], procAssign, m_domain[ilev] );
	    m_grids[ilev] = DisjointBoxLayout();
	    m_grids[ilev].define(a_vectBoxes[ilev], procAssign);
	    if (m_params.m_includeSurfaceSolver)
	    {

	    }
	}
    }
}
/*********************/
void
EBAMRWatershed::defineEBISLs()
{
    CH_TIME("EBAMRWatershed::defineEBISLs");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::defineEBISLs" << endl;
    }
    int numEBGhost = 4;//number of ghost cells in EBISL 
    //define 3D ebisl 
    for (int ilev = 0; ilev<= m_finestLevel; ilev++)
    {
	m_eblg[ilev] = EBLevelGrid(m_grids[ilev], m_domain[ilev], numEBGhost, m_ebisPtr);
	m_ebisl[ilev] = m_eblg[ilev].getEBISL();
    }
}
/*********************/
void 
EBAMRWatershed::defineExtraEBISLs()
{
    //update 3D grid information in the subsurface
    m_richardsCompGrid->defineGrids(m_finestLevel, m_eblg,
				    m_params.m_refRatio, 
				    true /*include irregular cells into mapping*/,
				    m_params.m_includeSUBSurfaceSolver/*include regular cells*/);
    m_richardsCompGrid->setVectDXs(m_vectDx);
    
    //printing currentl ebisl info 
    if (m_params.m_verbosity > 0) printTotalGridPoints();

    return;
}

/*********************/

void
EBAMRWatershed::initialGrid(const Vector<Vector<Box> >& a_vectBoxes)
{
    CH_TIME("EBAMRWatershed::initialGrid");
    if (m_params.m_verbosity >= 3)
    {
	pout () << "EBAMRWatershed::initialGrid "  << endl;
    }
    defineGrids(a_vectBoxes);
    defineEBISLs();
    defineExtraEBISLs();
    defineVariables();
    defineExtraTerms();
}
/*********************/
void
EBAMRWatershed::defineVariables(const int a_startLevel)
{
    CH_TIME("EBAMRWatershed::defineVariables");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::defineVariable" << endl;
    }
    int startLevel = 0;
    if (a_startLevel > startLevel)
    {
	startLevel=a_startLevel;
    }
    for (int ilev = startLevel; ilev <= m_finestLevel; ilev++)
    {
	EBCellFactory ebcellfact(m_ebisl[ilev]);
	m_psi_current[ilev]->define(m_grids[ilev], 1,  m_numGhostCells*IntVect::Unit, ebcellfact);
    
	if (m_params.m_includeSurfaceSolver)
	{
	    LayoutData<IntVectSet> irregSets(m_grids[ilev]);
	    const EBISLayout& ebisl = m_eblg[ilev].getEBISL();
	    for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
	    {
		//has to correspond to number of ghost cells
		Box grownBox = grow(m_grids[ilev].get(dit()), m_numGhostCells);
		grownBox &= m_domain[ilev];
		irregSets[dit()] = ebisl[dit()].getIrregIVS(grownBox);
	    }
	    BaseIVFactory<Real>  baseivfact(ebisl, irregSets);
	    m_EBSource[ilev]->define(m_grids[ilev], 1,  m_numGhostCells*IntVect::Unit, baseivfact);
	}
    }
    EBAMRDataOps::setToZero(m_psi_current);
    EBAMRDataOps::setCoveredVal(m_psi_current,0.0);
      
    if (m_params.m_includeSurfaceSolver)
    {
	m_ibc->setBoundarySource(m_EBSource);
    }
}
/*********************/
void
EBAMRWatershed::
defineExtraTerms(const int a_startLevel)
{
    CH_TIME("EBAMRWatershed::defineExtraTerms");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::defineExtraTerms" << endl;
    }
    int startLevel = 0;
    if (a_startLevel > startLevel)
    {
	startLevel=a_startLevel;
    }

    //set BaseIF implicit function reference for depth calculation
    m_richardsCompGrid->setImplicitBaseIFPtr(m_params.m_implicitBaseIF);

    //define parameters for surface and subsurface solvers
    m_richardsCompGrid->defineParameters(startLevel, m_finestLevel, m_numGhostCells);
  

    if (m_params.m_includeSurfaceSolver)
    {
	//define exchange fluxes in BC
	// m_ibc->defineIrregularData(startLevel, m_finestLevel, m_numGhostCells);
    }

    //define solvers operators 
    //in richards solver we additionally define sonductivity operator
    m_richardsCompGrid->defineOperators(0, m_finestLevel, m_numGhostCells, m_ibc->getScalarBC(), m_ibc->getScalarEBBC(), m_params.m_includeSurfaceSolver);

  
    if (m_params.m_includeSurfaceSolver)
    {

    }
}
/*********************/
void
EBAMRWatershed::initialData()
{
    CH_TIME("EBAMRWatershed::initialData");
    if (m_params.m_verbosity >= 3)
    {
	pout () << "EBAMRWatershed::initialData "  << endl;
    }
  
    ParserFunc initPFunc(m_params.m_initial_pressurehead);
    setParserFuncValue(m_psi_current, initPFunc, 0, 0.0);
    if (m_params.m_includeSurfaceSolver)
    {
    }
    initialExtraData();
}
/*********************/
void
EBAMRWatershed::initialExtraData()
{
    CH_TIME("EBAMRWatershed::initialExtraData");
    m_richardsCompGrid->setInitialParameters();
  
    if (m_params.m_includeSurfaceSolver)
    {
    }
}

void
EBAMRWatershed::postInitialize()
{
    CH_TIME("EBAMRWatershed::postInitialize");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::postInitialize" << endl;
    }
}
/*********************/
void
EBAMRWatershed::setInitialDt(Real a_dt)
{
    CH_TIME("EBAMRWatershed::setInitialDt");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::setInitialDt" << endl;
    }
    m_dt = a_dt;
    if (m_dt>m_params.m_maxDt) m_dt=m_params.m_maxDt;
}
/*********************/
Real
EBAMRWatershed::computeDt()
{
    CH_TIME("EBAMRWatershed::computeDt");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::computeDt " << endl;
    }

    Real dt;
    dt = m_dt*m_params.m_dtGrowFactor;
    if (dt>m_params.m_maxDt) dt=m_params.m_maxDt;
  
    return dt;
}
/*********************/
void
EBAMRWatershed::regrid(const Vector<Vector<Box> >& a_newGrids)
{
    CH_TIME("EBAMRWatershed::regrid");

    if (m_params.m_verbosity >= 3)
    {
	pout() << "EBAMRWatershed::regrid " << endl;
    }

    Interval intervScalar(0, 0);
//    Interval intervScalar2(1, 1);

    Vector<LevelData<EBCellFAB>* > tempLDPtr;

    tempLDPtr.resize(m_finestLevel+1);
    for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
	EBCellFactory ebcellfact(m_ebisl[ilev]);
	tempLDPtr[ilev] = new LevelData<EBCellFAB>();
	tempLDPtr[ilev]->define(m_grids[ilev], 1,  m_numGhostCells*IntVect::Unit, ebcellfact);
	m_psi_current[ilev]->copyTo(intervScalar, *(tempLDPtr[ilev]), intervScalar);
    }

    cacheExtraTerms();
  
    //this changes m_grids and m_ebisl
    defineGrids(a_newGrids);
    defineEBISLs();
    defineExtraEBISLs();
  
    //this redefines new data with new set of grids
    //this can also change m_finestLevel
    defineVariables();
    defineExtraTerms();
    //now fill new data holders, copying from old over old grids and interpolating where there is new grid
    tempLDPtr[0]->copyTo(intervScalar, *m_psi_current[0], intervScalar);
  
  
    for (int ilev=1; ilev<= m_finestLevel; ilev++)
    {
	//interpolate everywhere
	EBPWLFineInterp ebInterpScalar(m_grids[ ilev  ],
				       m_grids[ ilev-1],
				       m_ebisl[ ilev  ],
				       m_ebisl[ ilev-1],
				       m_domain[ilev-1],
				       m_params.m_refRatio[ilev-1],
				       1,
				       m_ebisPtr);
	ebInterpScalar.interpolate(*m_psi_current[ilev  ],
				   *m_psi_current[ilev-1],
				   intervScalar);
	tempLDPtr[ilev]->copyTo(intervScalar, *m_psi_current[ilev], intervScalar);
      
    } //end loop over levels

    for (int ilev=0; ilev<= m_finestLevel; ilev++)
    {
	delete tempLDPtr[ilev];
    }
  
    interpolateExtraTerms();
}


void
EBAMRWatershed::cacheExtraTerms(const int a_endLevel)
{
    m_richardsCompGrid->cacheParameters();
}


void
EBAMRWatershed::interpolateExtraTerms(const int a_startLevel, const int a_endLevel)
{
    m_richardsCompGrid->interpolateParameters();
}
/**********/
void
EBAMRWatershed::
tagCellsLevel(IntVectSet& a_tags, const int a_level)
{
    CH_TIME("EBAMRWatershed::tagCellsLevel");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::tagCellsLevel" << endl;
    }
    CH_assert(m_isSetup);
    
    //make a very simple refinement criteria for now
    LevelData<EBCellFAB> grad;
    //  computeGradient(grad, a_level);
    a_tags.makeEmpty();
    IFData2<3> ifdata(*m_params.m_implicitBaseIF);

    RealVect domainRealSize = m_domain[0].size();
    
    for (int idir=0; idir<3; idir++)
    {
	domainRealSize[idir] *=m_vectDx[0][idir];
    }

    ParserFunc tagCellsFunc;
    if (m_params.m_tagCellsFunc !="")
    {
	tagCellsFunc.define(m_params.m_tagCellsFunc);
    }
    const EBISLayout& ebisl = m_eblg[a_level].getEBISL();
    for (DataIterator dit = m_grids[a_level].dataIterator(); dit.ok(); ++dit)
    {
	const EBISBox& ebisBox = ebisl[dit()];
	//      const EBCellFAB& gradFAB = grad[dit()];
	const Box& grid =        m_grids[a_level].get(dit());
	const EBGraph& ebgraph = ebisBox.getEBGraph();
	IntVectSet ivsTot(grid);
	//if (ivsTot.contains(IntVect(52,95,74))) pout()<<" contains_tag="<<ivsTot.contains(IntVect(52,95,74))<<endl;

	for (VoFIterator vofit(ivsTot, ebgraph); vofit.ok(); ++vofit)
	{
	    const VolIndex& vof = vofit();
	    const IntVect& iv = vof.gridIndex();

	    if (m_params.m_tagCellsFunc !="")
	    {
		RealVect point = 0.5*RealVect::Unit;
		point += iv;

		if (ebisBox.isIrregular(iv) && !m_params.m_irregular_value_cc)
		{
		    //point += ebisBox.centroid(vof);
		    point += ebisBox.bndryCentroid(vof);
		}
		point *= m_vectDx[a_level];
		
		Real depth;
		if (ebisBox.isIrregular(iv))
		{
		    // RealVect distVec = ebisBox.bndryCentroid(vof);
		    // distVec *=m_vectDx[ilev];
		    
		    // depth  = -distVec.dotProduct(getAnisotropicNormal(ebisBox.normal(vof), m_vectDx[ilev]));
		    if (ebisBox.bndryArea(vof)>1.e-12)
			depth = 0.0;
		    else
			depth = -m_params.m_implicitBaseIF->value(point/m_vectDx[0])*m_vectDx[0][2];
		}
		else
		{
		    depth = -m_params.m_implicitBaseIF->value(point/m_vectDx[0])*m_vectDx[0][2];
		}
		if (depth < 0.0) depth = 0.0;
		//	if (ebisBox.isIrregular(iv)) pout()<<vof<<"\t"<<point<<" \t depth="<<depth<<"\t centroid="<<ebisBox.centroid(vof)<<endl;

		ParserInput parserinput;
		parserinput.point = point;
		parserinput.time  = m_time;
		parserinput.depth = depth;
		parserinput.psi   = (*m_psi_current[a_level])[dit()](vof,0);
		parserinput.phi   = (*m_richardsCompGrid->m_Phi[a_level])[dit()](vof,0);
		if (ebisBox.isIrregular(iv) && m_params.m_includeSurfaceSolver)
		{
		    parserinput.slope_x = (*m_richardsCompGrid->m_domainSlopes[a_level])[dit()](vof,0);
		    parserinput.slope_y = (*m_richardsCompGrid->m_domainSlopes[a_level])[dit()](vof,1);
		    parserinput.manncoef = (*m_richardsCompGrid->m_manningCoeff[a_level])[dit()](vof,0); 
		}
		const bool isTagged = tagCellsFunc.Eval(parserinput);
		if (isTagged>0.5) a_tags |= iv;
	    }

	    if(ebisBox.isIrregular(iv) && !m_params.m_refineAllIrreg && m_params.m_refineUnderresolved)
	    {
		Vector<RealVect> points;
		Vector<IntVect>  sides;
		Vector<PairPoints> pairPoints;
		const int numPoints=m_richardsCompGrid->findIntersectionPoints(ifdata, points, sides,vof, a_level, m_params.m_refRatio);
		const int numPairs=m_richardsCompGrid->findPairPoints(pairPoints,points,sides,iv);
		if (numPairs!=numPoints)
		{
		    a_tags |= iv;
		} 
		else if (!numPairs)
		{
		    a_tags |= iv;
		    // pout()<< "Tag Level = "<<a_level<<"\t vof= "<< vof<<endl;
		    // pout()<< "isIrregular="<<ebisBox.isIrregular(iv)
		    // 	  <<"\t isCovered="<<ebisBox.isCovered(iv)
		    // 	  <<"\t isMultiValued="<<ebisBox.isMultiValued (iv)
		    // 	  <<"\t bndryAreaEbis="<<ebisBox.bndryArea(vof)
		    // 	  <<"\t normalEbis="<<ebisBox.normal(vof)
		    // 	  <<"\t volFrac="<<ebisBox.volFrac(vof)
		    // 	  <<endl;
		} 
	    }
	   
	    // if (iv==IntVect(52,95,74))
	    // 	pout()<< " check   isIrregular="<<ebisBox.isIrregular(iv)
	    // 	      <<"\t isCovered="<<ebisBox.isCovered(iv)
	    // 	      <<"\t isMultiValued="<<ebisBox.isMultiValued (iv)
	    // 	      <<"\t bndryAreaEbis="<<ebisBox.bndryArea(vof)
	    // 	      <<"\t normalEbis="<<ebisBox.normal(vof)
	    // 	      <<"\t volFrac="<<ebisBox.volFrac(vof)
	    // 	      <<endl;

	} //end loop over vofs
	
	if (m_params.m_refineAllIrreg)
	{      
	    //refine all irregular cells
	    IntVectSet irregIVS = ebgraph.getIrregCells(grid);
	    a_tags |= irregIVS;
	}

    }
    a_tags.grow(m_params.m_tagBuffer);

}
/*****************/

void
EBAMRWatershed::
defineIrregularData()
{
    CH_TIME("EBAMRWatershed::defineIrregularData");
    if (m_params.m_includeSurfaceSolver)
    {
	//setting m_domainSlopes  
	m_richardsCompGrid->setDomainSlopes();
	//define stencil for surface kin wave operator
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::defineUnstructuredSurfaceGrid" << endl;
    }
	m_richardsCompGrid->defineUnstructuredSurfaceGrid();
    }
}
/*****************/
void
EBAMRWatershed::
allocateTemporaries()
{
    CH_TIME("EBAMRWatershed::allocateTemporaries");
    //number of active levels
    int numLevels = m_finestLevel+1;
    if (m_params.m_includeSUBSurfaceSolver)
    {
	m_psi_new.resize(numLevels, NULL);
	m_rhs_richards.resize(numLevels, NULL);
    }
    else
    {
	m_psiIrr_new.resize(numLevels, NULL);
	m_rhsIrr_richards.resize(numLevels, NULL);
    }
    
    if (m_params.m_includeSurfaceSolver)
    {

    }
    
    //allocate storage
    if (m_params.m_includeSUBSurfaceSolver)
    {
	for (int ilev=0; ilev < numLevels; ilev++)
	{
	    // LevelData<EBFluxFAB> testFlux;
	    // EBFluxFactory ebfluxfact(m_ebisl[ilev]);
	    // testFlux.define(m_grids[ilev],  1,  IntVect::Unit, ebfluxfact);
	    
	    EBCellFactory ebcellfact(m_ebisl[ilev]);
	    m_psi_new[ilev]        = new LevelData<EBCellFAB>(m_grids[ilev], 1, m_numGhostCells*IntVect::Unit, ebcellfact);
	    m_rhs_richards[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], 1, m_numGhostCells*IntVect::Unit, ebcellfact);
	    //copy current to temporary variable
	    //      EBLevelDataOps::clone(*m_psi_new[ilev], *m_psi_current[ilev]);
	    m_psi_current[ilev]->copyTo(*m_psi_new[ilev]);
	    m_psi_new[ilev]->exchange();

	    EBLevelDataOps::setToZero(*m_rhs_richards[ilev]);
	    
	}
	//m_richardsCompGrid->applyExtrapolation(m_richardsCompGrid->m_bndryPsiold, m_psi_current);
    }
    else
    {
	for (int ilev = 0; ilev <= m_finestLevel; ilev++)
	{
	     m_psiIrr_new[ilev] = new LevelData<BaseIVFAB<Real> >();
	     m_rhsIrr_richards[ilev] = new LevelData<BaseIVFAB<Real> >();

	    LayoutData<IntVectSet> irregSets(m_grids[ilev]);
	    const EBISLayout& ebisl = m_eblg[ilev].getEBISL();
	    for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
	    {
		//has to correspond to number of ghost cells
		//Box grownBox = m_grids[ilev].get(dit());
		Box grownBox = grow(m_grids[ilev].get(dit()), m_numGhostCells);
		grownBox &= m_domain[ilev];
		irregSets[dit()] = ebisl[dit()].getIrregIVS(grownBox);
	    }
	    BaseIVFactory<Real>  baseivfact(ebisl, irregSets);
	    m_psiIrr_new[ilev]->define(m_grids[ilev], 1,  m_numGhostCells*IntVect::Unit, baseivfact);
	    m_rhsIrr_richards[ilev]->define(m_grids[ilev], 1,  m_numGhostCells*IntVect::Unit, baseivfact);

	    for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)     
	    {  
		(*m_rhsIrr_richards[ilev])[dit()].setVal(0.0);
	    }
	    
	}
	copyEBCellFAB2BIVF(m_psi_current, m_psiIrr_new, m_eblg);
	
    }
    if (!m_curStep)
    {
	if (m_params.m_includeSurfaceSolver && m_params.m_printOutputRate)
	{
	    m_richardsCompGrid->printOutputRate(m_psi_current, m_time);
	}
    }
 
	
    if (m_params.m_includeSurfaceSolver)
    {
	
    }
}
/*****************/
void
EBAMRWatershed::
deleteTemporaries()
{
    CH_TIME("EBAMRWatershed::deleteTemporaries");
    //number of active levels
    int numLevels = m_finestLevel+1;
    //clean up storage
    for (int ilev=0; ilev<numLevels; ilev++)
    {
	if (m_params.m_includeSUBSurfaceSolver)
	{
	    delete m_psi_new[ilev];
	    delete m_rhs_richards[ilev];
	    m_psi_new[ilev]           = NULL;
	    m_rhs_richards[ilev]    = NULL;

	}
	else
	{
	    delete m_psiIrr_new[ilev];
	    delete m_rhsIrr_richards[ilev];
	    m_psiIrr_new[ilev]           = NULL;
	    m_rhsIrr_richards[ilev]    = NULL;
	}
	
	if (m_params.m_includeSurfaceSolver)
	{
	    
	}
    }
}
/*****************/
void
EBAMRWatershed::
advance()
{

    CH_TIME("EBAMRWatershed::advance");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::advance, nstep " << m_curStep
	       << ", starting time = "  << m_time
	       << ", dt = " << m_dt << endl;
    }


    //allocate space for integration
    allocateTemporaries();

    allocateExtraTemporaries();

    predictor();
    corrector();
    
    pointsUpdated();
    
    //remove all unnecessary data
    deleteTemporaries();
    deleteExtraTemporaries();
    
}
/******************/
void
EBAMRWatershed::
predictor()
{
    updateModelParameters();
    if (m_params.m_includeSUBSurfaceSolver)
    { 
	coupledSurfaceSubsurfaceSolver();
    }
    else
    {
	standAloneSurfaceSolver();
    }
    
}
/******************/
void
EBAMRWatershed::
updateModelParameters()
{
    //setting rhs for newton-krylov solver as zero
    EBAMRDataOps::setToZero(m_rhs_richards);
    
    //update parameters for subsurface at next time step fori mplicity solver
    m_richardsCompGrid->updateParameters(m_time+m_dt);
    
    if (m_params.m_includeSurfaceSolver)
    {
	//adding rainsource as a constant for test now. 
	ParserFunc sourceFunc(m_params.m_scalarSurfaceSourceFunc);
	setParserFuncValue(m_EBSource, sourceFunc, 0, m_time+m_dt, 1.0);// or m_dt); ???
	m_richardsCompGrid->projectFluxOnXYPlane(m_EBSource);
	//m_richardsCompGrid->applyExtrapolation(m_richardsCompGrid->m_bndryPsiold, m_psi_current);
    }
}
/******************/
void
EBAMRWatershed::
coupledSurfaceSubsurfaceSolver()
{
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::watershedSolver" << endl;
    }
    //CH_TIMERS("EBAMRWatershed::watershedSolver");
    Vector<Vector<LevelData<EBCellFAB>*>*> SUBsurfaceCompositeVector;
    SUBsurfaceCompositeVector.push_back(&m_psi_new);
    SUBsurfaceCompositeVector.push_back(&m_rhs_richards);
    //SUBsurfaceCompositeVector.push_back(&m_psi_current);
    m_richardsCompGrid->setCurrentPressure(m_psi_current);

    //solving kinetic wave equation for the surface
    bool converged = false;
    // static bool resetDiverged = false;//only initial time step 
    
    while(!converged)
    {
	converged = m_watershedSolver->solve(SUBsurfaceCompositeVector, m_vectDx, m_dt);
	//pout()<<"converged="<<converged<<endl;
	if (!converged)
	{
	    if( m_params.m_useSmallerDt)
	    {
		pout()<<"EBAMRWatershed::Solver did not converged. Trying to decrease time step"<<endl;
		m_dt /= m_params.m_dtGrowFactor;
		updateModelParameters();
		if (m_params.m_verbosity > 3)
		{
		    pout() << "EBAMRWatershed::rerun advance, nstep " << m_curStep
			   << ", starting time = "  << m_time
			   << ", dt = " << m_dt << endl;
		}
		//if(resetDiverged)
		{
		    for (int ilev=0; ilev <= m_finestLevel; ilev++)
		    {
			m_psi_current[ilev]->copyTo(*m_psi_new[ilev]);
			m_psi_new[ilev]->exchange();
		    }
		}
	    }
	    else
	    {
		pout()<<"EBAMRWatershed::Solver did not converged, terminating the run. \n Set parameter rerun_smaller_dt if you want to continue with using a smaller time step."<<endl;
		std::exit(1);
	    }

	}
	//resetDiverged = true;
    }
}

void
EBAMRWatershed::
standAloneSurfaceSolver()
{
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::watershedSolver" << endl;
    }
    //CH_TIMERS("EBAMRWatershed::watershedSolver");
    Vector<Vector<LevelData<BaseIVFAB<Real> >*>*> surfaceCompositeVector;
    surfaceCompositeVector.push_back(&m_psiIrr_new);
    surfaceCompositeVector.push_back(&m_rhsIrr_richards);
    
    m_richardsCompGrid->setCurrentPressure(m_psi_current);

    //solving kinetic wave equation for the surface
    bool converged = false;
    //static bool resetDiverged = false;//only initial time step 
    
    while(!converged)
    {
	converged = m_watershedSolver->solveSurfaceFlow(surfaceCompositeVector, m_vectDx, m_dt);
	if (!converged)
	{
	    if (m_params.m_useSmallerDt)
	    {
		pout()<<"EBAMRWatershed::Solver did not converged. Trying to decrease time step"<<endl;
		m_dt /= m_params.m_dtGrowFactor;
		updateModelParameters();
		if (m_params.m_verbosity > 3)
		{
		    pout() << "EBAMRWatershed::rerun advance, nstep " << m_curStep
			   << ", starting time = "  << m_time
			   << ", dt = " << m_dt << endl;
		}
		//if(resetDiverged)
		{
		    copyEBCellFAB2BIVF(m_psi_current, m_psiIrr_new, m_eblg);
		}
	    }
	    else
	    {
		pout()<<"EBAMRWatershed::Solver did not converged, terminating the run. \n Set parameter rerun_smaller_dt if you want to continue with using a smaller time step."<<endl;
		PetscFinalize();
#ifdef CH_MPI
		MPI_Finalize();
#endif
		std::exit(1);
	    }
	}
    }
}

/******************/
void
EBAMRWatershed::corrector()
{
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::corrector()" << endl;
    }
    
    if (m_params.m_includeSUBSurfaceSolver)
    {
	
	//copying solution as a current
	int nlev = m_finestLevel + 1;
	for (int ilev = 0; ilev < nlev; ilev++)
	{
	    m_psi_new[ilev]->copyTo(*m_psi_current[ilev]);
	    m_psi_current[ilev]->exchange();
	}
    }
    else
    {
	copyBIVF2EBCellFAB(m_psiIrr_new, m_psi_current, m_eblg);
    }

}



void
EBAMRWatershed::postTimeStep()
{
    CH_TIME("EBAMRWatershed::postTimeStep");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::postTimeStep" << endl;
    }
    if (m_params.m_includeSurfaceSolver && m_params.m_printOutputRate)
    {
	m_richardsCompGrid->printOutputRate(m_psi_current, m_time, m_dt);
    }
}


void
EBAMRWatershed::pointsUpdated()
{
    CH_TIME("EBAMRWatershed::pointsUpdated");
    
    int numLevels = m_finestLevel + 1;
    for (int ilev = 0; ilev < numLevels; ilev++)
    {
	for (LayoutIterator lit = m_grids[ilev].layoutIterator(); lit.ok(); ++lit)
	{
	    m_pointsUpdated += m_grids[ilev][lit()].numPts();
	}
      
	if (m_params.m_includeSurfaceSolver)
	{
	}
    }
}


/*****************/
#ifdef CH_USE_HDF5
/*****************/
void
EBAMRWatershed::writePlotFile()
{
    CH_TIME("EBAMRWatershed::writePlotFile");

    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::writePlotFile" << endl;
    }
    int curNumLevels = m_finestLevel + 1;

    int nlev = m_finestLevel + 1;
    Vector<LevelData<EBCellFAB>* > outputData(nlev, NULL);

    Vector<string> hnames(1, string("pressureHead"));
    Vector<string> hdnames(1, string("hydraulicHead"));
    Vector<string> Swnames(1, string("saturation"));
    Vector<string> Fluxnames(3);
    for (int idir = 0; idir < 3; idir++)
    {
	char fluxchar[100];
	sprintf(fluxchar, "DarcyFlux-%d", idir);
	Fluxnames[idir] = string(fluxchar);
    }
    Vector<string> Kxnames(3);
    for (int idir = 0; idir < 3; idir++)
    {
	char Kxchar[100];
	sprintf(Kxchar, "AbsPermbl-%d", idir);
	Kxnames[idir] = string(Kxchar);
    }
    Vector<string> Ssnames(1, string("SpecStorage"));
    Vector<string> Phinames(1, string("porosity"));
    Vector<string> Slopenames(2);
    Slopenames[0] = "Slope-x";
    Slopenames[1] = "Slope-y";

    //  Vector<string> Rhonames(1, string("density"));
  
    Vector<string> names;
    names = hnames;
    names.append(hdnames);
    names.append(Swnames);
    names.append(Fluxnames);
    names.append(Kxnames);
    names.append(Ssnames);
    names.append(Phinames);
    names.append(Slopenames);
    
    // names.append(Rhonames);
  
    //For adding extra data
    //  names.append(extraNames());
  
    char fileChar[1000];
    int ncells = m_domain[0].size(0);
    sprintf(fileChar, "plotSubsurface.nx%d.step.%07d.hdf5", ncells, m_curStep);
  
    bool replaceCovered = false;
    Vector<Real> coveredValues;

    //get Darcy Flux
    Vector<LevelData<EBCellFAB>* >  darcyFlux(nlev, NULL);
    Vector<LevelData<EBCellFAB>* >  Slopes   (nlev,NULL);
    for (int ilev = 0; ilev < nlev; ilev++)
    {
	darcyFlux[ilev] = new LevelData<EBCellFAB>();
	Slopes   [ilev] = new LevelData<EBCellFAB>();
	
	EBCellFactory ebcellfact(m_ebisl[ilev]);
	darcyFlux[ilev]->define(m_grids[ilev], 3,  IntVect::Zero, ebcellfact);
	Slopes   [ilev]->define(m_grids[ilev], 2,  IntVect::Zero, ebcellfact);
    }
    EBAMRDataOps::setToZero(darcyFlux);
    m_richardsCompGrid->getCCDarcyFlux(darcyFlux, m_psi_current);
    EBAMRDataOps::setToZero(Slopes);
    if (m_params.m_includeSurfaceSolver) copyBIVF2EBCellFAB(m_richardsCompGrid->m_domainSlopes, Slopes, m_eblg);
   
    for (int ilev = 0; ilev < nlev; ilev++)
    {
	int nvar = 13; 
	int counter=0;      
	//For adding extra data
	//      int numExtraVars = getNumExtraVars();
	// nvar += numExtraVars;
      
	EBCellFactory ebcellfact(m_ebisl[ilev]);
	outputData[ilev] = new LevelData<EBCellFAB>(m_grids[ilev], nvar, IntVect::Zero, ebcellfact);
 	LevelData<EBCellFAB> tempData(m_grids[ilev], 1,  m_numGhostCells*IntVect::Unit, ebcellfact);   
     
	Interval srcInterv, dstInterv;
	srcInterv = Interval(0, 0);
	dstInterv = Interval(counter, counter);
	counter++;
	m_psi_current[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

	dstInterv = Interval(counter, counter);
	EBLevelDataOps::axby(tempData, *m_psi_current[ilev], *(m_richardsCompGrid->m_z[ilev]), 1.0, 1.0);
        tempData.copyTo(srcInterv, *outputData[ilev], dstInterv);
	counter++;
      
	
	m_richardsCompGrid->computeSaturation(tempData,*m_psi_current[ilev], m_domain[ilev]); 
	srcInterv = Interval(0, 0);
	dstInterv = Interval(counter, counter);
	counter++;
	tempData.copyTo(srcInterv, *outputData[ilev], dstInterv);
	
 
	srcInterv = Interval(0, 2);
	dstInterv = Interval(counter, counter+2);
	counter +=3;
	darcyFlux[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);
	

	srcInterv = Interval(0, 2);
	dstInterv = Interval(counter, counter+2);
	counter +=3;
	m_richardsCompGrid->m_Kx[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);
      
	srcInterv = Interval(0, 0);
	dstInterv = Interval(counter, counter);
	counter++;
	m_richardsCompGrid->m_Ss[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);
      
	srcInterv = Interval(0, 0);
	dstInterv = Interval(counter, counter);
	counter++;
	m_richardsCompGrid->m_Phi[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

	srcInterv = Interval(0, 1);
	dstInterv = Interval(counter, counter+1);
	counter +=2;
	Slopes[ilev]->copyTo(srcInterv, *outputData[ilev], dstInterv);

	setCoveredStuffToZero(*outputData[ilev]);
    }

    Vector<IntVect> vectRatios;
    vectRatios.resize(nlev);
 
    for (int ilev = 1; ilev < nlev; ilev++)
    {
	vectRatios[ilev-1] = m_params.m_refRatio[ilev-1]*IntVect::Unit;
    }

    string filename1(fileChar);
    writeEBHDF5(filename1,
		m_grids,
		outputData,
		names,
		m_domain[0].domainBox(),
		m_vectDx[0],
		m_dt,
		m_time,
		vectRatios,
		curNumLevels,
		replaceCovered,
		coveredValues);
    for (int ilev = 0; ilev <nlev; ilev++)
    {
	delete outputData[ilev];
	delete darcyFlux[ilev];
	delete Slopes[ilev];
    }
  
    if (m_params.m_includeSurfaceSolver)
    {
    }
  
    addExtraDatatoPlotFile(filename1);
}


void
EBAMRWatershed::writeCheckpointFile()
{
    CH_TIME("EBAMRWatershed::writeCheckpointFile");
    CH_assert(m_isSetup);
    //Setup the level header information
    HDF5HeaderData header;
    
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRNoSubcycle::writeCheckpointFile" << endl;
    }
    
    //all the stuff in m_params had to come in at define time
    
    //bool conversion to int
    
    header.m_real["time"]                   = m_time;
    header.m_int ["cur_step"]               = m_curStep;
    header.m_int ["finest_level"]           = m_finestLevel;

    char iter_str[1000];

    int ncells = m_domain[0].size(0);
    sprintf(iter_str, "check%d.nx%d.%dd.hdf5", m_curStep, ncells, SpaceDim);
    
    HDF5Handle handleOut(iter_str, HDF5Handle::CREATE);
    //Write the header for this level
    header.writeToFile(handleOut);
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	handleOut.setGroupToLevel(ilev);
	write(handleOut,m_grids[ilev]);
	write(handleOut,*m_psi_current[ilev],"psi");
    }
    handleOut.close();
}


void
EBAMRWatershed::readCheckpointFile(const string& a_restartFile)
{
    CH_TIME("EBAMRWatershed::readCheckpointFile");
    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRWatershed::readCheckpointFile" << endl;
    }

    HDF5Handle handleIn(a_restartFile, HDF5Handle::OPEN_RDONLY);
    HDF5HeaderData header;
    header.readFromFile(handleIn);

    //all the stuff in m_params had to come in at define time
    m_time          =   header.m_real["time"]                   ;
    m_curStep       =   header.m_int ["cur_step"]               ;
    m_finestLevel   =   header.m_int ["finest_level"]           ;

    int finestLevelFromParmParse   =   m_params.m_maxLevel;
    if (m_finestLevel > finestLevelFromParmParse)
    {
	m_finestLevel = finestLevelFromParmParse;
    }

    //get all the grids
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	handleIn.setGroupToLevel(ilev);
	//Get the grids
	Vector<Box> vboxGrids;
	const int gridStatus = read(handleIn, vboxGrids);
	if (gridStatus != 0)
        {
	    MayDay::Error("readCheckpointLevel: file has no grids");
        }

	Vector<int> proc_map;
	LoadBalance(proc_map,  vboxGrids);
	//EBEllipticLoadBalance(proc_map,vboxGrids, m_domain[ilev], false, m_ebisPtr );

	m_grids[ilev]= DisjointBoxLayout(vboxGrids,proc_map);
    }
    
    //define stuff using grids
    defineEBISLs();
    defineExtraEBISLs();
    defineVariables();
    defineExtraTerms();
    
    //now input the actual data
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	handleIn.setGroupToLevel(ilev);
	read<EBCellFAB>(handleIn, *m_psi_current[ilev], "psi", m_grids[ilev], Interval(), false);
    }
    handleIn.close();
    
}
/*****************/
#endif //CH_USE_HDF5
/*****************/
void
EBAMRWatershed::
setupForRestart(const string& a_restartFile)
{
    CH_TIME("EBAMRWatershed::setupForRestart");

    if (m_params.m_verbosity > 3)
    {
	pout() << "EBAMRNoSubcycle::setupForRestart" << endl;
    }
    m_isSetup = true;
    m_doRestart  = true;
#ifdef CH_USE_HDF5
    readCheckpointFile(a_restartFile);
#else
    MayDay::Error("cannot restart from checkpoint without hdf5");
#endif
    int finestLevelFromParmParse = m_params.m_maxLevel;
    int oldFinestLevel = m_finestLevel;
    Interval intervScalar(0, 0);

    if (oldFinestLevel < finestLevelFromParmParse)
    {
	//cache data on old levels
	Vector<LevelData<EBCellFAB>* > tempLDPtr;
	tempLDPtr.resize(oldFinestLevel+1);
	for (int ilev=0; ilev<= oldFinestLevel; ilev++)
        {
	    EBCellFactory ebcellfact(m_ebisl[ilev]);
	    tempLDPtr[ilev] = new LevelData<EBCellFAB>();
	    tempLDPtr[ilev]->define(m_grids[ilev], 1,  m_numGhostCells*IntVect::Unit, ebcellfact);
	    m_psi_current[ilev]->copyTo(intervScalar, *(tempLDPtr[ilev]), intervScalar);
        }
	
	
	m_finestLevel = finestLevelFromParmParse;
	
	bool moreLevels = (m_finestLevel > oldFinestLevel);
	
	BRMeshRefine meshrefine;
	if (moreLevels)
        {
	    meshrefine.define(m_domain[0],              m_params.m_refRatio,
			      m_params.m_fillRatio,     m_params.m_blockFactor,
			      m_params.m_nestingRadius, m_params.m_maxBoxSize);
        }

	Vector<Vector<Box> > new_boxes;
	Vector<Vector<Box> > old_boxes(m_finestLevel+1);
	
	for (int ilev=0; ilev<m_finestLevel+1; ilev++)
        {
	    if (ilev>oldFinestLevel)
            {
		old_boxes[ilev] = Vector<Box>(1, m_domain[ilev].domainBox());
            }
	    else
            {
		old_boxes[ilev] = m_grids[ilev].boxArray();
            }
        }
	
	int base_level = oldFinestLevel;
	int top_level = m_finestLevel-1;
	
	IntVectSet tagsVect;
	tagCellsLevel(tagsVect, 0);

	int new_finest = meshrefine.regrid(new_boxes, tagsVect,
					   base_level, top_level,
					   old_boxes);
	if (new_finest != m_finestLevel)
        {
	    MayDay::Error("EBAMRNoSubcycle::setupForRestart -- new_finest not equal m_finestLevel");
        }
	
	//from defineGrids
	for (int ilev=oldFinestLevel+1; ilev<=m_finestLevel; ilev++)
        {
	    Vector<int> procAssign;
	    mortonOrdering((Vector<Box>&)(new_boxes[ilev]));
	    LoadBalance(procAssign,  new_boxes[ilev]);
	    //EBEllipticLoadBalance(procAssign,  new_boxes[ilev], m_domain[ilev], false, m_ebisPtr );
	    m_grids[ilev] = DisjointBoxLayout();
	    m_grids[ilev].define(new_boxes[ilev], procAssign);
        }
	
	int startLevel = oldFinestLevel+1;
	defineEBISLs();
	defineExtraEBISLs();
	defineVariables(startLevel);
	defineExtraTerms(startLevel);
	
	for (int ilev=0; ilev<=oldFinestLevel; ilev++)
        {
	    tempLDPtr[ilev]->copyTo(intervScalar, *m_psi_current[ilev], intervScalar);
        }
	
	for (int ilev=oldFinestLevel+1; ilev<=m_finestLevel; ilev++)
        {
	    //interpolate everywhere
	    EBPWLFineInterp ebInterpScalar(m_grids[ ilev  ],
					   m_grids[ ilev-1],
					   m_ebisl[ ilev  ],
					   m_ebisl[ ilev-1],
					   m_domain[ilev-1],
					   m_params.m_refRatio[ilev-1],
					   1,
					   m_ebisPtr);
	    ebInterpScalar.interpolate(*m_psi_current[ilev  ],
				       *m_psi_current[ilev-1],
				       intervScalar);
	    
        }
	
	defineIrregularData();
	postInitialize();
    }
    else  if (oldFinestLevel > finestLevelFromParmParse)
    {
	m_finestLevel = finestLevelFromParmParse;
	defineIrregularData();
	postInitialize();
    }
    else//no change in finest level from inputs
    {
	defineIrregularData();
	postInitialize();
    }
    //setup parameters for start time
    m_richardsCompGrid->setInitialParameters(m_time);    
}

void
EBAMRWatershed::setCoveredStuffToZero(LevelData<EBCellFAB>& a_h)
{
    CH_TIME("setCoveredStuffToZero");
    for (DataIterator dit = a_h.dataIterator(); dit.ok(); ++dit)
    {
	EBCellFAB&  hFAB =a_h[dit()];
	Real covVal = 0.0;
	for (int icomp = 0; icomp < hFAB.nComp(); icomp++)
        {
	    hFAB.setCoveredCellVal(covVal, icomp);
        }
    }
}

void
EBAMRWatershed::printTotalGridPoints()
{
    long long totalPoints    = 0;
    long long totalBoxes     = 0;
    int numLevels = m_finestLevel + 1;
    for (int ilev = 0; ilev < numLevels; ilev++)
    {
	long long pointsThisLevel = 0;
	for (LayoutIterator lit = m_grids[ilev].layoutIterator(); lit.ok(); ++lit)
	{
	    pointsThisLevel    += m_grids[ilev][lit()].numPts();
	}
	totalBoxes    += m_grids[ilev].size();
	
	
	if (m_params.m_includeSurfaceSolver)
	{
	    
	}
	pout() << "getAllIrregRefineLayouts:level[" << ilev
	       << "],  number of 3D boxes = " << m_grids[ilev].size();
      
	pout()<< ",  number of 3D points = " << pointsThisLevel;

	pout()<< endl;
    }
    pout() << "getAllIrregRefineLayouts:"
	   <<  "   total 3D boxes = " << totalBoxes;
  
    pout()<<  ", total points = " << totalPoints <<  endl;
}

void
EBAMRWatershed::setParserFuncValue(Vector<LevelData<EBCellFAB>*>&   a_var, ParserFunc& a_parserfunc, const int a_comp, const Real a_time, const Real a_scale)
{
    Interval interv(a_comp,a_comp);
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	const DisjointBoxLayout& dbl    = m_grids[ilev];
	for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
        {
	    const Box& myBox = dbl.get(dit());
	    const EBISLayout& ebisl = m_eblg[ilev].getEBISL();

	    const EBISBox& ebisBox = ebisl[dit()];
	    IntVectSet ivsTot(myBox);
	    for (VoFIterator vofit(ivsTot, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
		
		const VolIndex& vof = vofit();
		const IntVect& iv = vof.gridIndex();
		RealVect point = 0.5*RealVect::Unit;
		point += iv;

		if (ebisBox.isIrregular(iv) && !m_params.m_irregular_value_cc)
		{
		    //point += ebisBox.centroid(vof);
		    point += ebisBox.bndryCentroid(vof);
		}
		point *= m_vectDx[ilev];

		Real depth;
		if (ebisBox.isIrregular(iv))
		{
		    // RealVect distVec = ebisBox.bndryCentroid(vof);
		    // distVec *=m_vectDx[ilev];

		    // depth  = -distVec.dotProduct(getAnisotropicNormal(ebisBox.normal(vof), m_vectDx[ilev]));
		    if (ebisBox.bndryArea(vof)>1.e-12)
			depth = 0.0;
		    else
			depth = -m_params.m_implicitBaseIF->value(point/m_vectDx[0])*m_vectDx[0][2];
		}
		else
		{
		    depth = -m_params.m_implicitBaseIF->value(point/m_vectDx[0])*m_vectDx[0][2];
		}
		if (depth < 0.0) depth = 0.0;
		//	if (ebisBox.isIrregular(iv)) pout()<<vof<<"\t"<<point<<" \t depth="<<depth<<"\t centroid="<<ebisBox.centroid(vof)<<endl;

		ParserInput parserinput;
		parserinput.point = point;
		parserinput.time  = a_time;
		parserinput.depth = depth;
		parserinput.psi   = (*m_psi_current[ilev])[dit()](vof,0);
		parserinput.phi   = (*m_richardsCompGrid->m_Phi[ilev])[dit()](vof,0);
		if (m_params.m_includeSurfaceSolver && ebisBox.isIrregular(iv))
		{    
		    parserinput.slope_x = (*m_richardsCompGrid->m_domainSlopes[ilev])[dit()](vof,0);
		    parserinput.slope_y = (*m_richardsCompGrid->m_domainSlopes[ilev])[dit()](vof,1);
		    parserinput.manncoef = (*m_richardsCompGrid->m_manningCoeff[ilev])[dit()](vof,0); 
		}
		(*a_var[ilev])[dit()](vof,a_comp) = a_scale*a_parserfunc.Eval(parserinput);
		
	    }
	}
	a_var[ilev]->exchange(interv);
    }
}

void
EBAMRWatershed::setParserFuncValue(Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&   a_var, 
				   ParserFunc& a_parserfunc, 
				   const int a_comp, 
				   const Real a_time, 
				   const Real a_scale)
{
    Interval interv(a_comp,a_comp);
    //const int nghostTemp = 2;
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	const EBISLayout& ebisl = m_eblg[ilev].getEBISL();

	for (DataIterator dit = m_grids[ilev].dataIterator(); dit.ok(); ++dit)
	{
	    const EBISBox& ebisBox = ebisl[dit()];
	    //has to correspond to number of ghost cells
	    Box grid = m_grids[ilev].get(dit());
	    IntVectSet ivsIrreg = ebisl[dit()].getIrregIVS(grid);
	    for (VoFIterator vofit(ivsIrreg, ebisl[dit()].getEBGraph()); vofit.ok(); ++vofit)
	    {
		const VolIndex& vof = vofit();
		const IntVect& iv = vof.gridIndex();
		RealVect point = 0.5*RealVect::Unit;
		point += iv;

		if (!m_params.m_irregular_value_cc) point += ebisBox.bndryCentroid(vof);
		point *= m_vectDx[ilev];
		ParserInput parserinput;
		parserinput.point = point;
		parserinput.time  = a_time;
		parserinput.psi   = (*m_psi_current[ilev])[dit()](vof,0);
		parserinput.phi   = (*m_richardsCompGrid->m_Phi[ilev])[dit()](vof,0); 
		parserinput.slope_x = (*m_richardsCompGrid->m_domainSlopes[ilev])[dit()](vof,0);
		parserinput.slope_y = (*m_richardsCompGrid->m_domainSlopes[ilev])[dit()](vof,1);
		parserinput.manncoef = (*m_richardsCompGrid->m_manningCoeff[ilev])[dit()](vof,0); 

		(*a_var[ilev])[dit()](vof,a_comp) = a_scale*a_parserfunc.Eval(parserinput);
	    }
	}
	a_var[ilev]->exchange(interv);
    }
}
