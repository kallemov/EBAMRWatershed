#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBPETSCCOMPOSITEGRIDRICHARDS_H_
#define _EBPETSCCOMPOSITEGRIDRICHARDS_H_

#include <stack>
#include "BaseIF.H"
#include "EBPetscMap.H"
#include "EBQuadCFInterp.H"
#include "EBCoarseAverage.H"
#include "EBCoarseAverage.H"
#include "RichardsDomainBC.H"
#include "RichardsEBBC.H"
#include "EBDarcyOp.H"
#include "EBQuadCFStencilColor.H"
#include "IFData2.H"

#include "NamespaceHeader.H"

typedef struct ParserExpr ParserExpr;
struct ParserExpr{
    std::string SpecStorage;
    std::vector<std::string> SatConductivity;
    std::string Density;
    std::string Porosity;
    std::string Viscosity;
    std::vector<std::string> surfaceBC; 
    std::string ManningCoeff;
    std::vector<std::string> frictionSlopes;
};

// This struct contains all geometry information for
// 2D face of unstructured surface grid
//since we use 3D vector to store 2D vector 
//the third component is used to pair data
//faceData = outer normal[2] + length
typedef struct UGData UGData;
struct UGData{
    VolIndex vof=VolIndex();//connected vof
    RealVect faceData=RealVect::Zero;    
    Vector<Real> slopes;    
    int bndryDir=6;
    bool operator==(const UGData& data) const
    {
	return (vof==data.vof
		&& faceData  == data.faceData
		&& bndryDir == data.bndryDir
		//   && numFaces==data.numFaces
	    );
    }
//    RealVect face2dCentroid = RealVect::Zero;
//    Real interpWeight    = 0.;
};

typedef struct PairPoints PairPoints;
struct PairPoints
{
    RealVect first;
    RealVect second;
//    IntVect iv;
    Vector<IntVect> alternIVs;
};

    enum surfaceSolverType
{
    KinematicWaveUnstructured,
    DiffusionWaveUnstructured,
};

/*!
 * \brief constructor of Kinetic Wave Grid helper class
 *
 */
class EBPetscCompGridRichards : public EBPetscMap
{
public:

    void  fillRandom(Vector<LevelData<EBCellFAB>* >& a_psi);

    /*!
     * \brief base class constructor
     *
     */
    EBPetscCompGridRichards(const int a_maxLevel, const bool a_includeSurfaceSolver);

    virtual ~EBPetscCompGridRichards();

    /*!
     * \brief allocates model's parameters such as domain slopes, Manning's coefficients, etc.
     *
     */
    void allocateParameters();

    /*!
     * \brief define model's parameters for current grid
     *
     */
    void defineParameters(const int a_startLevel, 
			  const int a_finestLevel,
			  const int a_numGhostCells);

    /*!
     * \brief saves anisotropic dx on all levels
     *
     */
    void setVectDXs(const Vector<RealVect>& a_vectDxs)
	{
	    m_vectDxs=a_vectDxs;
	}


    /*!
     * \brief sets values of model's parameters at initialization
     *
     */
    void setInitialParameters(const Real a_time = 0.0);

    /*!
     * \brief update values of model's parameters at specific time
     *
     */
    void updateParameters(const Real a_time);

    /*!
     * \brief saves references for the current pressure on all levels
     *
     */
    void setCurrentPressure(const Vector<LevelData<EBCellFAB>* >& a_psi_current);

    
    /*!
     * \brief saves references for the implicit BaseIF function from geometryshop
     *
     */
    void setImplicitBaseIFPtr(BaseIF* a_implicitBaseIF);
    
    /*!
     * \brief caches parameter values for regriding
     *
     */
    void cacheParameters();

    /*!
     * \brief interpolates parameters after regriding and destroys cached terms
     *
     */
    void interpolateParameters();

    /*!
     * \brief defines model's interpolation and averaging operators
     *
     */
    void defineOperators(const int a_startLevel, 
			 const int a_finestLevel,
			 const int a_numGhostCells,
			 RefCountedPtr<RichardsDomainBCFactory> a_bcRichards,
			 RefCountedPtr<RichardsEBBCFactory>     a_ebbcRichards,
			 const bool includeOverlandFlow);
  
    /*!
     * \brief averages EBFluxFAB data down from fine to coarse 
     *
     */
    void averageDown(Vector<LevelData<EBFluxFAB>* >&  a_data);

    /*!
     * \brief averages EBCellFAB data down from fine to coarse 
     *
     */
    void averageDown(Vector<LevelData<EBCellFAB>* >&  a_data);

    void averageDownHydrolic(Vector<LevelData<EBCellFAB>* >&  a_data);

    /*!
     * \brief computes vector for Newton-Krylov step 
     *
     */
    void getNewtonKrylovVectorOp(Vector<LevelData<EBCellFAB>* >& a_F_mfree,
				 Vector<LevelData<EBCellFAB>* >& a_psi_mfree, 
				 const Vector<RealVect>& a_vectDx,
				 const Real a_dt);
    
    /*!
     * \brief computes Jacobinan-vector product for Newton-Krylov step of Richards equation
     *
     */
    void getNewtonKrylovJacobianOp(Vector<LevelData<EBCellFAB>* >& a_Jac_out,
				   Vector<LevelData<EBCellFAB>* >& a_Jac_in,
				   const Vector<LevelData<EBCellFAB>* >& a_sur_mfree, 
				   const Vector<RealVect>& a_vectDx,
				   const Real a_dt);

    void getNewtonKrylovVectorOp(Vector<LevelData<BaseIVFAB<Real> >* >& a_F_mfree,
				 Vector<LevelData<BaseIVFAB<Real> >* >& a_psi_mfree, 
				 const Vector<RealVect>& a_vectDx,
				 const Real a_dt);

    /*!
     * \brief computes SUBsuface part of the global Jacobinan matrix
     *
     */
    void createSubsurfaceJacobianMatrix(Mat a_Jacobian,
					const Vector<LevelData<EBCellFAB>* >& a_F_mfree, 
					const Vector<RealVect>& a_vectDx,
					const Real a_dt);
    /*!
     * \brief generate vectors for coloring a Jacobian matrix
     *
     */
    void getJacobianColoringVectors(std::vector<PetscInt>&  a_rows,
				    std::vector<std::vector<PetscInt> >& a_collumns,
				    std::vector<PetscInt>& a_num_diag,
				    std::vector<PetscInt>& a_num_offdiag,
				    const bool a_include_SUBsurface,
				    const int a_num_comps=1);

    /*!
     * \brief generate vectors for a Jacobian matrix
     *
     */
    void getJacobianVectors(std::vector<PetscInt>&                  a_rows,
			    std::vector<std::vector<PetscInt> >&    a_collumns,
			    std::vector<std::vector<PetscScalar> >& a_data,
			    Vector<LevelData<EBCellFAB>* >&         a_psi, 
			    const Real                              a_dt);
    
    void computeSaturation(LevelData<EBCellFAB>& a_Sw, 
			   const LevelData<EBCellFAB>& a_psi,
			   const ProblemDomain& a_domain);

    void  getDivDarcyFlux(Vector<LevelData<EBCellFAB>* >& a_Div,
			  const Vector<LevelData<EBCellFAB>* >& a_psi);

    void getCCDarcyFlux(Vector<LevelData<EBCellFAB>* >& a_flux,
			const Vector<LevelData<EBCellFAB>* >& a_psi);

    void computeGradientPsi(const Vector<LevelData<EBCellFAB>* >& a_psi);

    void addAMRIrregExchangeFlux(const Vector<LevelData<EBCellFAB>* >& a_psinew,
				 const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_source,
				 const Real a_dt);

    void  addAMRIrregExchangeFlux(const Vector<LevelData<BaseIVFAB<Real> >*>& a_psinew,
				  const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_source,
				  const Real a_dt);

    void setDomainSlopes();

    void defineUnstructuredSurfaceGrid();
    

    void printOutputRate(const Vector<LevelData<EBCellFAB>* >& a_psi, 
			 const Real a_time,
			 const Real a_dt=0.0);
    void  projectFluxOnXYPlane(Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_val);

    void addExchangeFlux(BaseIVFAB<Real>&       a_exflux,
			 const EBCellFAB& a_psinew,
			 const EBCellFAB& a_psiold,
			 const BaseIVFAB<Real>& a_domSl,
			 const BaseIVFAB<Real>& a_mannC,
			 const BaseIVFAB<Vector<UGData> >& a_UGData,
			 const BaseIVFAB<Real>& a_bndryArea,
			 const BaseIVFAB<Real>& a_2DGrad,
			 const ProblemDomain& a_domain,
			 const Real a_dt,
			 const RealVect& a_vectDx,
			 const Box&       a_box,
			 const EBISBox& a_ebisBox,
			 const bool isDiffusionWave,
			 const bool a_useRegularization = false,
			 const Real a_regConstant = 1.0);
    
    
    void  floorIrregular(Vector<LevelData<BaseIVFAB<Real> >* >& a_psi);

    int  findIntersectionPoints(IFData2<3>& a_ifdata,
				Vector<RealVect>& a_points,
				Vector<IntVect>&  a_sides,
				const VolIndex& a_vof,
				const int a_lev,
				const Vector<int>& a_refRatios,
				const bool a_skipDuplicates=true);
    
    int findPairPoints(Vector<PairPoints>& a_pairPoints,
		       const Vector<RealVect>& a_points,
		       const Vector<IntVect>&  a_sides,
		       const IntVect& a_iv);
    
    //protected:
    //Richards equation
    Real m_gconst;// gravitational constant (assuming acting along negative z direction)
    //cc data 
    Vector<LevelData<EBCellFAB>* > m_Ss;  // Specific Storage
    Vector<LevelData<EBCellFAB>* > m_Kx;  // Absolute permeability of the medium
    Vector<LevelData<EBCellFAB>* > m_Phi; // Porosity of the medium
    Vector<LevelData<EBCellFAB>* >           m_z; //depth from the surface
    Vector<RefCountedPtr<LevelData<EBFluxFAB> > >  m_DarcyFlux;


    Vector<LevelData<BaseIVFAB<Real> >* >    m_domainSlopes;
    Vector<LevelData<BaseIVFAB<Real> >* >    m_manningCoeff;
    Vector<LevelData<BaseIVFAB<Real> >* >    m_2DGrad;

    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_spreadFlux;
    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_bndryPsinew;
    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > > m_bndryPsiold;
    Vector<LevelData<BaseIVFAB<Real> >* >                m_bndryArea;
    Vector<LayoutData<BaseIVFAB<VoFStencil> >* >         m_2DGradStencil;
    Vector<LayoutData<BaseIVFAB<Vector<UGData> > >* >    m_surfaceGridData;

    //define these parameters as constants for now
    Real m_Ssat;  // relative saturated water content(degree of saturation)
    Real m_Sres;  // relative residual saturation
    Real m_alpha; // alpha soil parameter
    Real m_nsoil; // n soil parameter
    Real m_visc; // water viscosity
    Real m_lambda; // alpha soil parameter
    Real m_rho_zero; // scaling paprameter for density
    Real m_regConstant;
    
    bool m_includeSurfaceSolver;
    bool m_includeSUBSurfaceSolver;
    bool m_irregular_value_cc;
    bool m_useBCRegularization;

    surfaceSolverType m_surface_solver_type;
    
    BaseIF* m_implicitBaseIF; //implicit function from GeometryShop
private:
  
    /*!
     * \brief compute coefficients for Darcy fluxes(conductivity)
     *
     */
    void  computeKxKrFaceValues(const Vector<LevelData<EBCellFAB>* >& a_psi);

    /*!
     * \brief define stencil for irregular cells(surface solver)
     *
     */

    void  setUGStencil(Vector<UGData>&           a_vofStencil2D,
		       Vector<int>&              a_savedIndex,
		       IFData2<3>&               a_ifdata,
		       const Vector<PairPoints>& a_pairPoints,
		       const VolIndex&           a_vof,
		       const int                 a_lev,
		       const EBISBox&            a_ebisBox,
		       const ProblemDomain&      a_domain,
		       const ProblemDomain&      a_domainZero,
		       const RealVect&           a_vectDx);
    
    void setUnderresolvedUGStencil(Vector<UGData>&           a_vofStencil2D,
				   Vector<PairPoints>&       a_pairPoints,
				   BaseIVFAB<Real>&          a_bndryAres,
				   IFData2<3>&               a_ifdata,
				   const Vector<int>&        a_savedIndex,
				   const VolIndex&           a_vof,
				   const int                 a_lev,
				   const EBISBox&            a_ebisBox,
				   const ProblemDomain&      a_domain,
				   const ProblemDomain&      a_domainZero,
				   const RealVect&           a_vectDx);
    bool  computeClustersData(Vector<RealVect>&         a_normals,
			      Vector<Real>&             a_areas,
			      Vector<int>&              a_indexCluster, 
			      const Vector<PairPoints>& a_pairs,
			      const RealVect&           a_vectDx);
    
    void resolveCell(Vector<PairPoints>&       a_resolvedPairs,
		     IFData2<3>&               a_ifdata,
		     const IntVect&            a_vofIV,
		     const int                 a_lev,
		     const Vector<int>&        a_ratios);

    void  setUnderresolvedPairs(Vector<PairPoints>&       a_pairPoints,
				std::stack<RealVect>&     a_allStacks,
				Vector<int>&              a_stenIndex,
				IFData2<3>&               a_ifdata,
				const Vector<UGData>&     a_vofStencil2D,
				const Vector<int>&        a_savedIndex,
				const VolIndex&           a_vof,
				const int                 a_lev);
    
    void define2DGradStencil(const int a_startLevel, 
			     const int a_finestLevel);

    bool  setUGData(UGData&           a_ugdata,
		    const PairPoints& a_pair,
		    const VolIndex&   a_vofCenter,
		    const VolIndex&   a_vofConnect,
		    const EBISBox&    a_ebisBox,
		    const RealVect&   a_vectDx,
		    const int         a_bndryDir);

    int isPointOnCoarseEdge(IntVect& a_edgeInds,
			    const RealVect& a_point,
			    const RealVect& a_cellCorner,
			    const Real a_eps);

    
    bool isUGDataInStencil(const UGData&         a_ugdata,
			   const Vector<UGData>& a_vofStencil);

    void computeGradImplicitFunction(RealVect& a_grad,
				     const RealVect& a_point,
				     const int a_dim=2,
				     const Real delta=1.e-5);

    
    bool  isPointsConnected(std::stack<RealVect>& a_stack,
			    Vector<RealVect>& a_visited,
			    const RealVect& a_secondPoint,
			    const Vector<PairPoints>& a_resolvedPairs,
			    const RealVect& a_cellCorner);

    
    bool checkVoFPairs(const Vector<PairPoints>& a_checkPairs,
		       const PairPoints&         a_vofPair,
		       const IntVect&            a_vofIV);
    
    Real computeShoeLaceArea(const Vector<RealVect>&   a_points,
			     const RealVect&           a_vectDx);
    
    /*!
     * \brief sets elevation
     *
     */
    void  setZ(Vector<LevelData<EBCellFAB>* >& a_z);


    void getCFInterfaceColoringVectors(std::vector<PetscInt>&                         a_collumns,
				       int&                                           a_inside,
				       int&                                           a_outside,
				       const VolIndex&                                a_vof,
				       const DataIndex&                               a_dit,
				       const int                                      a_lev,
				       const EBISBox&                                 a_ebisBox,
				       const Vector<LevelData<BaseEBCellFAB<int> >*>& a_procIDGlobal,
				       const Vector<EBQuadCFStencilColor*>&           a_quadCFStencilColor,
				       const Vector<LevelData<BaseEBCellFAB<int64_t> >*>& a_bufferVofGlobalMap,
				       const Vector<LevelData<BaseEBCellFAB<int> >*>& a_bufferprocIDGlobal);
    
    void getSurfaceColoringVectors(std::vector<PetscInt>&                         a_collumns,
				   int&                                           a_inside,
				   int&                                           a_outside,
				   const VolIndex&                                a_vof,
				   const DataIndex&                               a_dit,
				   const int                                      a_lev,
				   const EBISBox&                                 a_ebisBox,
				   const Vector<LevelData<BaseEBCellFAB<int> >*>& a_procIDGlobal,
				   const Vector<EBQuadCFStencilColor*>&           a_quadCFStencilColor,
				   const Vector<LevelData<BaseEBCellFAB<int64_t> >*>& a_bufferVofGlobalMap,
				   const Vector<LevelData<BaseEBCellFAB<int> >*>& a_bufferprocIDGlobal);
    
    
    int m_meanKr_method; //identifier of averaging method for Kr 
    //fc data for Darcy fluxes
    Vector<RefCountedPtr<LevelData<EBFluxFAB> > >           m_KxKrRho;
    Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >    m_KxKrIVF;

    Vector<LevelData<EBCellFAB>* >           m_tempCC1;
    Vector<LevelData<EBFluxFAB>* >           m_tempKxRho;
    Vector<LevelData<EBFluxFAB>* >           m_tempFlux1;
    
    //this is reference to watershed solver
    Vector<LevelData<EBCellFAB>* > m_psi_saved;  // current step h

    int m_maxLevel;// max level
    int m_finestLevel;// max level

    //  Vector<RefCountedPtr<EBConstantCFInterp> > m_constantCFI;
    Vector<RefCountedPtr<EBQuadCFInterp> >    m_quadCFI;
    Vector<RefCountedPtr<EBCoarseAverage> >   m_aveOper;
    Vector<RefCountedPtr<EBCoarseAverage> >   m_aveSpac;
    Vector<RefCountedPtr<EBDarcyOp> >  m_DarcyOp;

    //storage to cache parameter value for interpolation
    Vector<Vector<LevelData<EBCellFAB>* > > tempLDVec;

    ParserExpr m_paramExpr;

    RefCountedPtr<RichardsDomainBC> m_DomainBC;

    Real m_time;
    Vector<RealVect> m_vectDxs;
};

#include "NamespaceFooter.H"

#endif
