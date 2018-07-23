#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRWATERSHED_H_
#define _EBAMRWATERSHED_H_

class EBAMRWatershedSolver;
class WatershedIBCFactory;
class WatershedIBC;

#include "EBWatershedUtils.H"
#include "EBPetscCompGridRichards.H"
#include "ParserFunc.H"
#include "LevelData.H"
#include "EBCellFAB.H"
#include "EBISLayout.H"
#include "EBLevelGrid.H"


#include "UsingNamespace.H"

/// Class to manage advance for EBAMRWatershed.
/**
 */
class EBAMRWatershed
{
    friend class EBPetscCompGridKinWave;
    friend class EBPetscCompGridRichards;
    friend class EBAMRWatershedSolver;
    friend class WatershedIBC;

public:
    // weak construction is bad
    EBAMRWatershed()
	{
	    MayDay::Error("invalid constructor");
	}

    ///
    /**
       Note that the EBAMRWatershed class defines its own grid hierarchy
       Initial and boundary conditions come through the WaterIBCFactory.
       If fixedGrids != NULL, regridding gets turned off and the input
       grids are used
    */
    EBAMRWatershed(const WatershedParameters&      a_params,
		   const WatershedIBCFactory&       a_IBC,
		   const ProblemDomain&      a_coarsestDomain,
		   const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance());



    //fixes timestep to input
    void setInitialDt(Real a_dt);

    /// destructor.
    virtual ~EBAMRWatershed();

    /// set up grids, etc for AMR run
    void setupForAMRRun();

    /// set up for restart
    void setupForRestart(const string& a_restartFile);

    /// advance solution.  returns time of final solution
    Real run(Real a_maxTime, int a_maxStep);

    ///for convergence tests
    Vector<DisjointBoxLayout>& getGrids()
	{
	    return m_grids;
	}
    ///for convergence tests
    Vector<EBISLayout>& getEBISLayouts()
	{
	    return m_ebisl;
	}

    Vector<RealVect>& getVectDx()
	{
	    return m_vectDx;
	}

    Vector<ProblemDomain> & getDomain()
	{
	    return m_domain;
	}
    /// write final checkpoint and plot files
    void conclude();

    /// accumulate number of points updated
    void pointsUpdated();

    /*!
     * \brief caches surface and subsurface solver parameters for regriding
     *
     *
     */
    virtual void cacheExtraTerms(const int a_endLevel = 9999999);
 
    /*!
     * \brief interpolates extra terms after regriding 
     *
     *   it is assumed that method destroys cached terms allocated in cacheExtraTerms once interpolation is done.
     *
     */
    virtual void interpolateExtraTerms(const int a_startLevel = -1, const int a_endLevel = 9999999);

    virtual void interpolateExtraTermsForRestart(const int a_startLevel = -1, const int a_endLevel = 9999999)
	{
	}

#ifdef CH_USE_HDF5
    ///
    void writeCheckpointFile();

    ///
    virtual void writeExtraDataToCheckpoint(HDF5Handle& a_handleOut, const int ilev)
	{
	}

    ///
    virtual void readExtraDataFromCheckpoint(HDF5Handle& a_handleIn, const int ilev)
	{
	}

    ///
    virtual void readCheckpointFile(const string& a_restartFile);

    ///
    virtual void writePlotFile();

    ///
    virtual void addExtraOutputData(LevelData<EBCellFAB>& a_outputData, const int a_startDstInterv, const int a_ilev)
	{
	}

    ///
    virtual void addExtraDatatoPlotFile(const string& a_filename)
	{
	}

    ///
    virtual int getNumExtraVars()
	{
	    return 0;
	}

    ///
    virtual Vector<string> extraNames()
	{
	    return 0;
	}

#endif

    virtual void advance();
    virtual void predictor();
    virtual void corrector();

    /// things to do after a timestep
    virtual void postTimeStep();

    /// create tags
    void tagCells(Vector<IntVectSet>& a_tags);

    /// create tags
    virtual void tagCellsLevel(IntVectSet& a_tags, const int a_level);

    /// regrid
    void regrid();

    ///
    /**
       Perform any post-regridding ops.   This includes the final smoothing and
       projecting the velocity and reinitializing the pressure.
    */
    void postRegrid();

    ///
    /**
       Perform any pre-regrdding ops.   This is vel = (I-mu*lapl)vel.
    */
    void preRegrid();

    ///
    void initialGrid(const Vector<DisjointBoxLayout>& a_new_grids);

    /// initialize data
    void initialData();

    /// initialize extra data
    virtual void initialExtraData();

    /// things do do after initialization
    void postInitialize();

    virtual void defineGrids(const Vector<Vector<Box> >& a_vectBoxes);
    // virtual void defineEBISLs(const int a_endLevel = 9999999);
    // virtual void defineExtraEBISLs(const int a_endLevel = 9999999){}
    virtual void defineEBISLs();
    virtual void defineExtraEBISLs();
  
    void defineVariables(const int a_startLevel=-1);

    ///
    virtual void defineExtraTerms(const int a_startLevel=-1);

    virtual void allocateDataHolders();
    /// compute dt
    virtual Real computeDt();


    virtual void updateModelParameters();
    virtual void coupledSurfaceSubsurfaceSolver();
    virtual void standAloneSurfaceSolver();

    /// set up the AMR hierarchy, etc for a_vectGrids.
    void setup();

    /// regrid  once new boxes have been generated
    void regrid(const Vector<Vector<Box> >& a_new_grids);

    void initialGrid(const Vector<Vector<Box> >& a_vectBoxes);

    //fills temporary stuff.  need to do this for some convergence tests.  do not do this unless you are sure
    //you know the context
    void allocateTemporaries();
    virtual void allocateExtraTemporaries()
	{
	}
    //delete temporary stuff.  need to do this for some convergence tests.  do not do this unless you are sure
    //you know the context
    void deleteTemporaries();
    virtual void deleteExtraTemporaries()
	{
	}
    void computeGradient(LevelData<EBCellFAB>& a_grad, int a_level);

protected:

    void setCoveredStuffToZero(LevelData<EBCellFAB>& a_h);
    void setParserFuncValue(Vector<LevelData<EBCellFAB>* >& a_var, 
			    ParserFunc& a_parserfunc, 
			    const int a_comp, 
			    const Real a_time,
			    const Real a_scale = 1.0);
    void setParserFuncValue(Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >&   a_var, 
			    ParserFunc& a_parserfunc, 
			    const int a_comp, 
			    const Real a_time, 
			    const Real a_scale = 1.0);
    //parameters
    WatershedParameters          m_params;
    WatershedIBC*                 m_ibc;
    //3D grid description
    Vector<ProblemDomain>     m_domain;
    Vector<RealVect>              m_vectDx;
    Vector<DisjointBoxLayout> m_grids;
    Vector<EBISLayout>        m_ebisl;
    Vector<EBLevelGrid>       m_eblg;
    const EBIndexSpace*    m_ebisPtr;

    //Vector<ProblemDomain>     m_domain_surface;
    //Vector<DisjointBoxLayout> m_grids_surface;
 
    EBPetscCompGridRichards    *m_richardsCompGrid; 
    EBAMRWatershedSolver *m_watershedSolver;

    //  finest extant level
    int m_finestLevel;
    int m_curStep;
    int m_numGhostCells;

    //  time stepping information
    Real m_time;
    Real m_dt;
    Real m_scaleFactorDt;
    Real m_maxDt;
 
    long long m_pointsUpdated;

    //switches that i would rather lose but cannot
    bool m_isSetup;
    bool m_doRestart;

    //data held in memory
    Vector<LevelData<EBCellFAB>* > m_psi_current; // pressure head
    //temporary data for time advance
    Vector<LevelData<EBCellFAB>* > m_psi_new; 
    Vector<LevelData<EBCellFAB>* > m_rhs_richards; // is cc rhs
    
    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_EBSource;

    //holders for surface flow only case
    Vector<LevelData<BaseIVFAB<Real> >* > m_psiIrr_new;
    Vector<LevelData<BaseIVFAB<Real> >* > m_rhsIrr_richards;

    //defines the above stuff
    void defineIrregularData();

private:
    /*!
     * \brief prints total points for each refinement level
     *
     */
    void printTotalGridPoints();

    //we disallow copy construction and assignment of objects
    //which contain pointered data
    EBAMRWatershed(const EBAMRWatershed& a_input)
	{
	    MayDay::Error("invalid constructor");
	}

    //we disallow copy construction and assignment of objects
    //which contain pointered data
    void operator=(const EBAMRWatershed& a_input)
	{
	    MayDay::Error("invalid constructor");
	}

};

//#include "NamespaceFooter.H"

#endif
