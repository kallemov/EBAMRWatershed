#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
//#include "EBDebugOut.H"
#include "EBAMRIO.H"
#include "ParmParse.H"
#include "Notation.H"
#include "IndexTM.H"
#include "EBAMRDataOps.H"
#include "EBLevelDataOps.H"
#include "EBQuadCFInterp.H"
#include "EBCellFactory.H"
#include "EBPWLFineInterp.H"
#include "EBSUBSurfaceOps.H"
#include "EBPetscCompGridRichards.H"
#include "ReductionOps.H"
#include "Stencil.H"
#include "ParserFunc.H"
#include "EBQuadCFStencilColor.H"
#include "BaseEBCellFactory.H"
#include "ScattBilinearInterp.H"

#include "NamespaceHeader.H"
// extern void init_scatt_bilinear(const double*, const double*, const double*);
// extern const double interp_scatt_bilinear(const double, const double);
// extern const double scatt_bilinear_derivative(const double, const double, const int);
extern void scatt_bilinear_derivative_stencil(double*, const double*, const double*, const int);
static const Real Tol = 1e-15;

typedef IndexTM<int,3>                      EdgeIndex;
typedef map<EdgeIndex,Real > EdgeIntersections;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

EBPetscCompGridRichards::EBPetscCompGridRichards(const int a_maxLevel, const bool a_includeSurfaceSolver)
  :EBPetscMap(),
   m_gconst(1.0),
   m_visc(1.0),
   m_lambda(0.5),
   m_rho_zero(1.0),
   m_includeSurfaceSolver(a_includeSurfaceSolver),
   m_irregular_value_cc(true),
   m_meanKr_method(0),
   m_maxLevel(a_maxLevel),
   m_time(0.0)
{
    CH_TIME("EBPetscCompGridRichards::EBPetscCompGridRichards");
  
  const int nlevels = m_maxLevel + 1;

  m_Ss.resize(nlevels, NULL);  
  m_Kx.resize(nlevels, NULL);  
  m_Phi.resize(nlevels, NULL); 
  //m_Rho.resize(nlevels, NULL); 

  //ref pointers
  m_DarcyFlux.resize(nlevels);
  m_KxKrRho.resize(nlevels);
  m_KxKrIVF.resize(nlevels);

  if  (m_includeSurfaceSolver)
  {
      m_spreadFlux.resize(nlevels);
      m_bndryPsinew.resize(nlevels);
      m_bndryPsiold.resize(nlevels);
      m_domainSlopes.resize(nlevels, NULL);
      m_manningCoeff.resize(nlevels, NULL);
      m_bndryArea.resize(nlevels, NULL);
      m_surfaceGridData.resize(nlevels, NULL);
  }
  m_z.resize(nlevels, NULL);
  m_tempCC1.resize(nlevels, NULL);
  m_tempKxRho.resize(nlevels, NULL);


  m_quadCFI.resize(nlevels);
  m_aveOper.resize(nlevels);
  m_aveSpac.resize(nlevels);
  m_DarcyOp.resize(nlevels);

  ParmParse pp;
  //Specific storage
  pp.get("Ss", m_paramExpr.SpecStorage);
  //Saturated hydraulic conductivity
  pp.getarr("Kx", m_paramExpr.SatConductivity, 0, 3);
  //Porosity
  pp.get("porosity", m_paramExpr.Porosity);
  //Density
  pp.get("rho", m_paramExpr.Density);

  if  (m_includeSurfaceSolver)
  {
      string tempString;
      pp.get("surfaceBC_Lo_X_func",tempString);
      m_paramExpr.surfaceBC.push_back(tempString);

      tempString.clear();
      pp.get("surfaceBC_Lo_Y_func",tempString);
      m_paramExpr.surfaceBC.push_back(tempString);

      tempString.clear();
      pp.get("surfaceBC_Hi_X_func",tempString);
      m_paramExpr.surfaceBC.push_back(tempString);

      tempString.clear();
      pp.get("surfaceBC_Hi_Y_func",tempString);
      m_paramExpr.surfaceBC.push_back(tempString);

      
      pp.get("ManningCoeff", m_paramExpr.ManningCoeff);
      pp.getarr("friction_slopes", m_paramExpr.frictionSlopes, 0, 2);
      pp.get("use_bc_regularization", m_useBCRegularization);
      if (m_useBCRegularization) pp.get("regularization_parameter", m_regConstant);
	  
      
      int surface_solver_type;
      pp.get("surface_solver_type", surface_solver_type);
      switch (surface_solver_type)
      {
      case 0: m_surface_solver_type=KinematicWaveUnstructured; break;
      case 1: m_surface_solver_type=DiffusionWaveUnstructured;
	  m_2DGradStencil.resize(nlevels, NULL);
	  m_2DGrad.resize(nlevels, NULL);
	  break;
      default: m_surface_solver_type=KinematicWaveUnstructured;
      }
  }
  //Viscosity
  pp.query("viscosity", m_visc);
 
  //vanGenuchten parameters
  //relative residual saturation
  pp.get("Sres", m_Sres);
  //relative saturated water content
  pp.get("Ssat", m_Ssat);
  //soil parameters
  pp.get("alpha_soil", m_alpha);
  pp.get("n_soil", m_nsoil);
  pp.query("lambda_soil", m_lambda);

 //gravitational constant
  pp.get("gravity", m_gconst);

  pp.get("permeability_averaging_method", m_meanKr_method);

  pp.query("irregular_value_cc", m_irregular_value_cc);

  allocateParameters();
}


EBPetscCompGridRichards::~EBPetscCompGridRichards()
{   
    CH_TIME("EBPetscCompGridRichards::~");
  EBPetscMap::clean();
  for (int ilev = 0; ilev <= m_maxLevel; ilev++)
  {
      delete m_Ss[ilev];  
      delete m_Kx[ilev];  
      delete m_Phi[ilev]; 

      delete m_z[ilev];
      delete m_tempCC1[ilev];
      delete m_tempKxRho[ilev];
      //delete m_tempFlux1[ilev];

      if  (m_includeSurfaceSolver)
      {
	  delete m_domainSlopes[ilev];
	  delete m_manningCoeff[ilev];
	  delete m_bndryArea[ilev];
	  delete m_surfaceGridData[ilev];
	  if (m_surface_solver_type==DiffusionWaveUnstructured)
	  {
	      delete m_2DGradStencil[ilev];
	      delete m_2DGrad[ilev];
	  }
      }
  }
  if (m_implicitBaseIF) delete m_implicitBaseIF;
}

void 
EBPetscCompGridRichards::allocateParameters()
{
    CH_TIME("EBPetscCompGridRichards::allocateParameters");

  for (int ilev = 0; ilev <= m_maxLevel; ilev++)
  {
      m_Ss[ilev]     = new LevelData<EBCellFAB>();
      m_Kx[ilev]     = new LevelData<EBCellFAB>();
      m_Phi[ilev]    = new LevelData<EBCellFAB>();
      //m_Rho[ilev]    = new LevelData<EBCellFAB>();

      //additional variables used in computations
      m_DarcyFlux[ilev] = RefCountedPtr<LevelData<EBFluxFAB> >(new LevelData<EBFluxFAB>());
      m_KxKrRho[ilev]   = RefCountedPtr<LevelData<EBFluxFAB> >(new LevelData<EBFluxFAB>());
      m_KxKrIVF[ilev]   = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >());

      m_z[ilev]         = new LevelData<EBCellFAB>();

      m_tempCC1[ilev]   = new LevelData<EBCellFAB>();
      m_tempKxRho[ilev] = new LevelData<EBFluxFAB>();
      //m_tempFlux1[ilev] = new LevelData<EBFluxFAB>();

      if  (m_includeSurfaceSolver)
      {
	  m_domainSlopes[ilev]      = new LevelData<BaseIVFAB<Real> >();
	  m_manningCoeff[ilev]      = new LevelData<BaseIVFAB<Real> >();
	  m_bndryArea[ilev]         = new LevelData<BaseIVFAB<Real> >();
	  m_surfaceGridData[ilev]   = new LayoutData<BaseIVFAB<Vector<UGData> > >();
	  m_spreadFlux[ilev]     = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >());
	  m_bndryPsinew[ilev]     = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >());
	  m_bndryPsiold[ilev]     = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(new LevelData<BaseIVFAB<Real> >());
	  if (m_surface_solver_type==DiffusionWaveUnstructured)
	  {
	      m_2DGradStencil[ilev]     = new LayoutData<BaseIVFAB<VoFStencil> >();
	      m_2DGrad[ilev]            = new LevelData<BaseIVFAB<Real> >();
	  }
      }
  }
}

void
EBPetscCompGridRichards::defineParameters(const int a_startLevel, 
					 const int a_finestLevel,
					 const int a_numGhostCells)
{
    CH_TIME("EBPetscCompGridRichards::defineParameters");

    m_finestLevel = a_finestLevel;
    int startLevel = 0;
    // if (a_startLevel > startLevel)
    // {
    //   startLevel=a_startLevel;
    // }
    for (int ilev = startLevel; ilev <= a_finestLevel; ilev++)
    {
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	EBCellFactory ebcellfact(m_eblgsPtr[ilev]->getEBISL());
      m_Ss[ilev]->define(dbl, 1,  a_numGhostCells*IntVect::Unit, ebcellfact);
      m_Kx[ilev]->define(dbl, SpaceDim,  a_numGhostCells*IntVect::Unit, ebcellfact);
      m_Phi[ilev]->define(dbl, 1,  a_numGhostCells*IntVect::Unit, ebcellfact);
      //m_Rho[ilev]->define(dbl, 1,  a_numGhostCells*IntVect::Unit, ebcellfact);

      //we need to define flux terms before setting conductivity operator
      int nghostTemp = a_numGhostCells; 
      m_z[ilev]->define(dbl,       1,  a_numGhostCells*IntVect::Unit, ebcellfact);
      m_tempCC1[ilev]->define(dbl, 1,  nghostTemp*IntVect::Unit, ebcellfact);
      
      int nghostFlux = a_numGhostCells;
      EBFluxFactory ebfluxfact(m_eblgsPtr[ilev]->getEBISL());
      m_KxKrRho[ilev]->define(dbl,  1,  nghostFlux*IntVect::Unit, ebfluxfact);
      m_tempKxRho[ilev]->define(dbl,  1,  nghostFlux*IntVect::Unit, ebfluxfact);
      m_DarcyFlux[ilev]->define(dbl,  1,  nghostFlux*IntVect::Unit, ebfluxfact);

      if  (m_includeSurfaceSolver)  
      {
	  m_surfaceGridData[ilev]    ->define(dbl);
	  if (m_surface_solver_type==DiffusionWaveUnstructured)
	  {
	      m_2DGradStencil[ilev]    ->define(dbl);
	  }
      }

      LayoutData<IntVectSet> irregSets(dbl);
      
      DataIterator dit = m_eblgsPtr[ilev]->getDBL().dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
          //has to correspond to number of ghost cells
	  Box box = grow(dbl.get(datInd), nghostTemp);
	  box &= m_eblgsPtr[ilev]->getDomain();
	  //const Box& box = dbl.get(datInd);
	  const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	  irregSets[datInd] = ebisl[datInd].getIrregIVS(box);
	  if (m_includeSurfaceSolver)
	  {
	      (*m_surfaceGridData[ilev])[datInd].define(irregSets[datInd],ebisl[datInd].getEBGraph(),1);
	      if (m_surface_solver_type==DiffusionWaveUnstructured)
	      {
		  (*m_2DGradStencil[ilev])[datInd].define(irregSets[datInd],ebisl[datInd].getEBGraph(),2);
	      }
	  }
      }
      BaseIVFactory<Real>  baseivfact(m_eblgsPtr[ilev]->getEBISL(), irregSets);
      m_KxKrIVF[ilev]->define(dbl, 1,  nghostTemp*IntVect::Unit, baseivfact);
      if  (m_includeSurfaceSolver)
      {
	  m_domainSlopes[ilev]->define(dbl, 2,  nghostTemp*IntVect::Unit, baseivfact);
	  m_manningCoeff[ilev]->define(dbl, 1,  nghostTemp*IntVect::Unit, baseivfact);
	  m_bndryArea[ilev]   ->define(dbl, 1,  nghostTemp*IntVect::Unit, baseivfact);
	  m_spreadFlux[ilev]  ->define(dbl, 1,  nghostTemp*IntVect::Unit, baseivfact);
	  m_bndryPsinew[ilev] ->define(dbl, 1,  nghostTemp*IntVect::Unit, baseivfact);
	  m_bndryPsiold[ilev] ->define(dbl, 1,  nghostTemp*IntVect::Unit, baseivfact);
	  if (m_surface_solver_type==DiffusionWaveUnstructured)
	  {
	      m_2DGrad[ilev]->define(dbl, 2,  nghostTemp*IntVect::Unit, baseivfact);
	  }
      }
    }
  // if  (m_includeSurfaceSolver)
  // {
  //     //setting m_domainSlopes  
  //     setDomainSlopes();
  //     //define stencil for surface kin wave operator
  //     defineUnstructuredSurfaceGrid(a_startLevel, a_finestLevel);
  //     if (m_surface_solver_type==DiffusionWaveUnstructured)
  //     {
  // 	  define2DGradStencil(a_startLevel, a_finestLevel);
  //     }
  // }

  
  //setting KxKrIVF to unity since ebbc fluxes are set explicitely
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      EBLevelDataOps::setToZero(*m_DarcyFlux[ilev]);
//      EBLevelDataOps::setVal(*m_DarcyFlux[ilev],5.0);
      EBLevelDataOps::setToZero(*m_KxKrRho[ilev]);      
      DataIterator dit = m_eblgsPtr[ilev]->getDBL().dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  (*m_KxKrIVF[ilev])[datInd].setVal(1.0);
	  if  (m_includeSurfaceSolver)
	  {
	      // (*m_GradPsi[ilev])[datInd].setVal(0.0);
	      //(*m_manningCoeff[ilev])[datInd].setVal(m_manningCoeffVal);
	      (*m_spreadFlux[ilev])[datInd].setVal(0.0);
	  }
      }
      
      m_KxKrIVF[ilev]->exchange();
      if  (m_includeSurfaceSolver)
      {
	  //m_GradPsi[ilev]->exchange();
	  //m_manningCoeff[ilev]->exchange();
	  m_spreadFlux[ilev]->exchange();
      }
  }

  //set z, a depth from the surface
  //we need it to be initialized here since we use it in initial pressure
  setZ(m_z);
  EBAMRDataOps::setToZero(m_tempKxRho);
  // writeEBFluxLDname(m_tempKxRho[0], 0,"dump2.hdf5");
}



void
EBPetscCompGridRichards::setInitialParameters(const Real a_time)
{
    //setZ(m_z);
    //setting constants
    updateParameters(a_time);
}

void
EBPetscCompGridRichards::setCurrentPressure(const Vector<LevelData<EBCellFAB>* >& a_psi_current)
{
  m_psi_saved = a_psi_current;
  return;
}

void
EBPetscCompGridRichards::setImplicitBaseIFPtr(BaseIF* a_implicitBaseIF)
{
    m_implicitBaseIF = a_implicitBaseIF;
}

void 
EBPetscCompGridRichards::cacheParameters()
{
  CH_TIME("EBPetscCompGridRichards::cacheParameters");
  //interval for scalar and vector
  Interval intervScalar(0, 0);
  Interval intervVec(0, SpaceDim-1);

  Vector<LevelData<EBCellFAB>* > tempLD1, tempLD2, tempLD3, tempLD4;
  tempLD1.resize(m_finestLevel+1);
  tempLD2.resize(m_finestLevel+1);
  tempLD3.resize(m_finestLevel+1);
  tempLD4.resize(m_finestLevel+1);
  for (int ilev=0; ilev<= m_finestLevel; ilev++)
  {
      const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
      EBCellFactory ebcellfact(m_eblgsPtr[ilev]->getEBISL());

      //cache vars, using only one ghost cell(?)
      const int ghostnum = 2;

      tempLD1[ilev] = new LevelData<EBCellFAB>();
      tempLD1[ilev]->define(dbl, 1,  ghostnum*IntVect::Unit, ebcellfact);
      m_Ss[ilev]->copyTo(intervScalar, *(tempLD1[ilev]), intervScalar);

      tempLD2[ilev] = new LevelData<EBCellFAB>();
      tempLD2[ilev]->define(dbl, SpaceDim,  ghostnum*IntVect::Unit, ebcellfact);
      m_Kx[ilev]->copyTo(intervVec, *(tempLD2[ilev]), intervVec);

      tempLD3[ilev] = new LevelData<EBCellFAB>();
      tempLD3[ilev]->define(dbl, 1,  ghostnum*IntVect::Unit, ebcellfact);
      m_Phi[ilev]->copyTo(intervScalar, *(tempLD3[ilev]), intervScalar);

      tempLD4[ilev] = new LevelData<EBCellFAB>();
      tempLD4[ilev]->define(dbl, 1,  ghostnum*IntVect::Unit, ebcellfact);
      //m_Rho[ilev]->copyTo(intervScalar, *(tempLD4[ilev]), intervScalar);
    }

  //save cached leveldata
  tempLDVec.push_back(tempLD1);
  tempLDVec.push_back(tempLD2);
  tempLDVec.push_back(tempLD3);
  tempLDVec.push_back(tempLD4);
}


void 
EBPetscCompGridRichards::interpolateParameters()
{
    CH_TIME("EBPetscCompGridRichards::interpolateParameters");

  //number of cached leveldata terms
  const int numCached = tempLDVec.size();

  //interval for scalar and vector
  Interval intervScalar(0, 0);
  Interval intervVec(0, SpaceDim-1);

  //fill zero level with cached values, since all parameter are already defined on a new grid
  tempLDVec[0][0]->copyTo(intervScalar, *m_Ss[0], intervScalar);
  tempLDVec[1][0]->copyTo(intervVec, *m_Kx[0], intervVec);
  tempLDVec[2][0]->copyTo(intervScalar, *m_Phi[0], intervScalar);
  //tempLDVec[3][0]->copyTo(intervScalar, *m_Rho[0], intervScalar);

  //update ghost cells
  m_Ss[0]->exchange();
  m_Kx[0]->exchange();
  m_Phi[0]->exchange();
  //m_Rho[0]->exchange();

  for (int ilev=1; ilev<= m_finestLevel; ilev++)
  {
      const DisjointBoxLayout& dbl     = m_eblgsPtr[ilev]->getDBL();
      const DisjointBoxLayout& dblCoar = m_eblgsPtr[ilev-1]->getDBL();

      //interpolate everywhere
      EBPWLFineInterp ebInterpScalar( dbl,
				      dblCoar,
				      m_eblgsPtr[ilev]->getEBISL(),
				      m_eblgsPtr[ilev-1]->getEBISL(),
				      m_eblgsPtr[ilev-1]->getDomain(),
				      m_refRatios[ilev-1],
				      1,
				      m_eblgsPtr[ilev]->getEBIS());

      EBPWLFineInterp ebInterpVec( dbl,
				   dblCoar,
				   m_eblgsPtr[ilev]->getEBISL(),
				   m_eblgsPtr[ilev-1]->getEBISL(),
				   m_eblgsPtr[ilev-1]->getDomain(),
				   m_refRatios[ilev-1],
				   SpaceDim,
				   m_eblgsPtr[ilev]->getEBIS());

      //iterpolate all parameters and copy cached data from current level
      ebInterpScalar.interpolate(*m_Ss[ilev  ],
				 *m_Ss[ilev-1],
				 intervScalar);
      tempLDVec[0][ilev]->copyTo(intervScalar, *m_Ss[ilev], intervScalar);


      ebInterpVec.interpolate(*m_Kx[ilev  ],
                              *m_Kx[ilev-1],
                              intervVec);
      tempLDVec[1][ilev]->copyTo(intervVec, *m_Kx[ilev], intervVec);
      

      ebInterpScalar.interpolate(*m_Phi[ilev  ],
				 *m_Phi[ilev-1],
				 intervScalar);
      tempLDVec[2][ilev]->copyTo(intervScalar, *m_Phi[ilev], intervScalar);


      // ebInterpScalar.interpolate(*m_Rho[ilev  ],
      // 				 *m_Rho[ilev-1],
      // 				 intervScalar);
      // tempLDVec[3][ilev]->copyTo(intervScalar, *m_Rho[ilev], intervScalar);

      //update ghost cells
      m_Ss[ilev]->exchange();
      m_Kx[ilev]->exchange();
      m_Phi[ilev]->exchange();
      // m_Rho[ilev]->exchange();
    }
  //set z
  //  setZ(m_z);

  //clean cached terms
  for (int ipart=0; ipart< numCached; ++ipart)
    {
      for (int ilev=0; ilev<= m_finestLevel; ilev++)
	{
	  delete tempLDVec[ipart][ilev];
	}
    }
  tempLDVec.clear();
}

void
EBPetscCompGridRichards::defineOperators(const int a_startLevel, 
					 const int a_finestLevel,
					 const int a_numGhostCells,
					 RefCountedPtr<RichardsDomainBCFactory> a_bcRichards,
					 RefCountedPtr<RichardsEBBCFactory>     a_ebbcRichards,
					 const bool a_includeOverlandFlow)
{
    CH_TIME("EBPetscCompGridRichards::defineOperators");

  int startLevel = 1;
  // if (a_startLevel > startLevel)
  //   {
  //     startLevel=a_startLevel;
  //   }
  for (int ilev = startLevel; ilev <= a_finestLevel; ilev++)
    {
      //always one component for quadcfi---only way to get reuse
      int nvarquad = 1;
      const DisjointBoxLayout& dbl     = m_eblgsPtr[ilev]->getDBL();
      const DisjointBoxLayout& dblCoar = m_eblgsPtr[ilev-1]->getDBL();

      m_quadCFI[ilev]  = RefCountedPtr<EBQuadCFInterp>(new  EBQuadCFInterp(dbl, dblCoar,
									   m_eblgsPtr[ilev]->getEBISL(), 
									   m_eblgsPtr[ilev-1]->getEBISL(),
									   m_eblgsPtr[ilev-1]->getDomain(),
									   m_refRatios[ilev-1], 
									   nvarquad,
									   *m_eblgsPtr[ilev]->getCFIVS(),
									   m_eblgsPtr[ilev]->getEBIS()));
      
      m_aveOper[ilev]  = RefCountedPtr<EBCoarseAverage>(new EBCoarseAverage(dbl, dblCoar,
									    m_eblgsPtr[ilev]->getEBISL(), 
									    m_eblgsPtr[ilev-1]->getEBISL(),
									    m_eblgsPtr[ilev-1]->getDomain(),
									    m_refRatios[ilev-1], 
									    nvarquad,
									    m_eblgsPtr[ilev]->getEBIS()));
      
      m_aveSpac[ilev]  = RefCountedPtr<EBCoarseAverage>(new EBCoarseAverage(dbl, dblCoar,
									    m_eblgsPtr[ilev]->getEBISL(), 
									    m_eblgsPtr[ilev-1]->getEBISL(),
									    m_eblgsPtr[ilev-1]->getDomain(),
									    m_refRatios[ilev-1], 
									    SpaceDim,
									    m_eblgsPtr[ilev]->getEBIS()));
    }

  //Now define conductivity operator for computing Darcy fluxes
  IntVect ghostVect = a_numGhostCells*IntVect::Unit;
  //  IntVect ghostVect = IntVect::Unit;

  Real beta = 1.0; // rescale it later for the whole sum //m_scalingFactor;//scales all terms in richards equation

  for (int ilev = 0; ilev <= a_finestLevel; ilev++)
  {

      RefCountedPtr<RichardsEBBC>      ebbcPtr(a_ebbcRichards->create(m_eblgsPtr[ilev]->getDomain(),
								      m_eblgsPtr[ilev]->getEBISL(), 
								      m_vectDxs[ilev],
								      &ghostVect, 
								      &ghostVect));
      if (a_includeOverlandFlow)
      {
	  ebbcPtr->setData(m_spreadFlux[ilev]);
      }
      
      RefCountedPtr<RichardsDomainBC> dombcPtr(a_bcRichards->create(m_eblgsPtr[ilev]->getDomain(),
								    m_eblgsPtr[ilev]->getEBISL(), 
								    m_vectDxs[ilev]));
      
      RefCountedPtr<LevelData<EBFluxFAB> >            darcyFlux  = RefCountedPtr<LevelData<EBFluxFAB       > >(m_DarcyFlux[ilev]);
      RefCountedPtr<LevelData<EBFluxFAB> >            bcoef      = RefCountedPtr<LevelData<EBFluxFAB       > >(m_KxKrRho[ilev]);
      RefCountedPtr<LevelData<BaseIVFAB<Real> > >     bcoefIrreg = RefCountedPtr<LevelData<BaseIVFAB<Real> > >(m_KxKrIVF[ilev]);
      
      int refToFiner   = -1;
      int refToCoarser = -1;
      if (ilev > 0)
      {
	refToCoarser= m_refRatios[ilev-1];
      }
      if (ilev < a_finestLevel)
      {
	refToFiner = m_refRatios[ilev];
      }
      EBLevelGrid dummyEBLG;
      m_DarcyOp[ilev]  = RefCountedPtr<EBDarcyOp> (
	  new EBDarcyOp( ilev < a_finestLevel ? *m_eblgsPtr[ilev+1] :dummyEBLG,
			 *m_eblgsPtr[ilev],
			 ilev > 0 ? *m_eblgsPtr[ilev-1] : dummyEBLG,
			 dombcPtr, ebbcPtr,  
			 m_quadCFI[ilev], 
			 refToFiner, refToCoarser,
			 m_vectDxs[ilev], beta, bcoef, bcoefIrreg, darcyFlux,
			 ghostVect, ghostVect));

    }

  //we need this to get access to muParser functions for computing boundary coefficients
  m_DomainBC = RefCountedPtr<RichardsDomainBC>(a_bcRichards->create(m_eblgsPtr[0]->getDomain(), m_eblgsPtr[0]->getEBISL(), m_vectDxs[0])); 
}

void 
EBPetscCompGridRichards::averageDown(Vector<LevelData<EBFluxFAB>* >&  a_data)
{
  CH_TIME("EBPetscCompGridRichards::averageDown(fluxData)");
  //do average down here
  for (int ilev = m_finestLevel; ilev > 0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      if (ilev==m_finestLevel) a_data[m_finestLevel]->exchange(interval);

      RefCountedPtr<EBCoarseAverage> avePtr = m_aveOper[ilev];
      if (ncomp == SpaceDim)
        {
          avePtr = m_aveSpac[ilev];
        }

      avePtr->average(*a_data[ilev-1],
		      *a_data[ilev  ],
		      interval);

      a_data[ilev-1]->exchange(interval);
    } //end loop over levels
}


void
EBPetscCompGridRichards::averageDownHydrolic(Vector<LevelData<EBCellFAB>* >&  a_data)
{
  CH_TIME("EBPetscCompGridRichards::averageDown(celldata)");
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      a_data[ilev]->exchange();
      const LevelData<EBCellFAB>& psi_LD  = *a_data[ilev];
      LevelData<EBCellFAB>& h_LD = *m_tempCC1[ilev];
      
      EBLevelDataOps::axby(h_LD, psi_LD, *m_z[ilev], 1.0, 1.0); //m_gconst/m_rho_zero);
      h_LD.exchange();
  }
  //do average down here
  for (int ilev = m_finestLevel; ilev > 0; ilev--)
  {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      if (ilev==m_finestLevel) a_data[m_finestLevel]->exchange(interval);
      RefCountedPtr<EBCoarseAverage> avePtr = m_aveOper[ilev];
      if (ncomp == SpaceDim)
      {
          avePtr = m_aveSpac[ilev];
      }
      avePtr->average(*m_tempCC1[ilev-1],
		      *m_tempCC1[ilev  ],
		      interval);
      m_tempCC1[ilev-1]->exchange(interval);
      EBLevelDataOps::axby(*a_data[ilev-1], *m_tempCC1[ilev-1], *m_z[ilev-1], 1.0, -1.0); //m_gconst/m_rho_zero);
      a_data[ilev-1]->exchange();
    } //end loop over levels

}
  

void
EBPetscCompGridRichards::averageDown(Vector<LevelData<EBCellFAB>* >&  a_data)
{
  CH_TIME("EBPetscCompGridRichards::averageDown(celldata)");
  //do average down here
  for (int ilev = m_finestLevel; ilev > 0; ilev--)
    {
      int ncomp = a_data[ilev]->nComp();
      Interval interval(0, ncomp-1);
      if (ilev==m_finestLevel) a_data[m_finestLevel]->exchange(interval);
      RefCountedPtr<EBCoarseAverage> avePtr = m_aveOper[ilev];
      if (ncomp == SpaceDim)
        {
          avePtr = m_aveSpac[ilev];
        }

      avePtr->average(*a_data[ilev-1],
		      *a_data[ilev  ],
		      interval);
      a_data[ilev-1]->exchange(interval);

    } //end loop over levels
}
void 
EBPetscCompGridRichards::getNewtonKrylovVectorOp(Vector<LevelData<EBCellFAB>* >& a_F_mfree,
						 Vector<LevelData<EBCellFAB>* >& a_psi_mfree, 
						 const Vector<RealVect>& a_vectDx,
						 const Real a_dt)
{
  //update fine-corse interface
  a_psi_mfree[0]->exchange();
  CH_TIME("EBPetscCompGridRichards::getNewtonKrylovVecotrOp");
  //EBAMRDataOps::setToZero(a_F_mfree);
  Interval interval(0, 0);
/*
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      const LevelData<EBCellFAB>& psi_LD  = *a_psi_mfree[ilev];
      LevelData<EBCellFAB>& h_LD = *m_tempCC1[ilev];
      
      EBLevelDataOps::axby(h_LD, psi_LD, *m_z[ilev], 1.0, 1.0); //m_gconst/m_rho_zero);
      h_LD.exchange();
  }

 for (int ilev = 1; ilev <= m_finestLevel; ilev++)
  {
      EBQuadCFInterp& ebinterpolate = *m_quadCFI[ilev];
      
      // ebinterpolate.interpolate(*a_psi_mfree[ilev  ],
      // 				*a_psi_mfree[ilev-1],
      // 				interval);
      
      // a_psi_mfree[ilev]->exchange();

      ebinterpolate.interpolate(*m_tempCC1[ilev  ],
      				*m_tempCC1[ilev-1],
      				interval);
      

      EBLevelDataOps::axby(*a_psi_mfree[ilev], *m_tempCC1[ilev], *m_z[ilev], 1.0, -1.0);
      a_psi_mfree[ilev]->exchange();

  }
*/
  
  //for (int ilev = 0; ilev <= m_finestLevel; ilev++)  m_DarcyFlux[ilev]->exchange();
  //compute Darcy Flux
  getDivDarcyFlux(a_F_mfree, a_psi_mfree);
  // static int count=0;
  // char fileChar[1000];
  // sprintf(fileChar, "vectorQQ.step.%d.hdf5", count);
  // writeEBAMRname(&a_F_mfree, fileChar);

      
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	//m_DarcyFlux[ilev]->exchange();

      LevelData<EBCellFAB>& psi_LD  = *a_psi_mfree[ilev];
      LevelData<EBCellFAB>& F_LD  = *a_F_mfree[ilev];
      
      //now add saturation terms and partial derivates
      const DisjointBoxLayout &dbl = a_psi_mfree[ilev]->getBoxes();
      const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();

      DataIterator dit = psi_LD.dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  EBCellFAB&         F_FAB     =  F_LD[datInd];
	  const EBCellFAB&   psi_FAB     =  psi_LD[datInd];
	  const EBCellFAB&   psi_oldFAB  =  (*m_psi_saved[ilev])[datInd];
	  const EBCellFAB&   Ss        = (*m_Ss[ilev])[datInd];
	  const EBCellFAB&   Phi       = (*m_Phi[ilev])[datInd];
	  // const EBCellFAB&   rho       = (*m_Rho[ilev])[datInd];

	  //F_FAB.mult(-a_dt/m_visc);
	  F_FAB.mult(-a_dt);

	  //const Box& grid = surFAB.getRegion();
	  const Box& grid = dbl.get(datInd);

	  EBCellFAB& Sw_FAB     =  (*m_tempCC1[ilev])[datInd];
	  Sw_FAB.setVal(0.0);

	  addSaturationTerms(Sw_FAB, psi_FAB, psi_oldFAB, Ss, Phi, 1.0, //rho, 
			     m_Ssat, m_Sres, m_alpha, m_nsoil, 
			     m_eblgsPtr[ilev]->getDomain(), grid);

	  const EBISBox& ebisBox = ebisl[datInd];
	  IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
	  for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      // const RealVect normal = getAnisotropicNormal(ebisBox.normal(vof), a_vectDx[ilev]);
	      // const Real VolFrac = getAnisotropicVolFrac(ebisBox.bndryCentroid(vof),
	      // 						 normal, 
	      // 						 a_vectDx[ilev]);
	      Sw_FAB(vof,0) *=ebisBox.volFrac(vof);// VolFrac;
	  }

	  F_FAB +=Sw_FAB;
      }//mybox
      a_F_mfree[ilev]->exchange();
    }
  // sprintf(fileChar, "vectorFF.step.%d.hdf5", count);
  // writeEBAMRname(&a_F_mfree, fileChar);

  //update fine-corse interface
  //averageDown(a_F_mfree);
  //averageDown((Vector<LevelData<EBFluxFAB>*>) m_DarcyFlux);

}

void 
EBPetscCompGridRichards::getNewtonKrylovVectorOp(Vector<LevelData<BaseIVFAB<Real> >* >& a_F_mfree,
						 Vector<LevelData<BaseIVFAB<Real> >* >& a_psi_mfree, 
						 const Vector<RealVect>& a_vectDx,
						 const Real a_dt)
{
    CH_TIME("EBPetscCompGridRichards::getNewtonKrylovVecotrOp");
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	a_psi_mfree[ilev]->exchange();
    }
    
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      LevelData<BaseIVFAB<Real> >& psi_LD  = *a_psi_mfree[ilev];
      LevelData<BaseIVFAB<Real> >& F_LD  = *a_F_mfree[ilev];
      
      //now add saturation terms and partial derivates
      const DisjointBoxLayout &dbl = a_psi_mfree[ilev]->getBoxes();
      const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();

      DataIterator dit = psi_LD.dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  BaseIVFAB<Real>&         F_FAB     =  F_LD[datInd];
	  const BaseIVFAB<Real>&   psi_FAB     =  psi_LD[datInd];
	  const EBCellFAB&   psi_oldFAB  =  (*m_psi_saved[ilev])[datInd];
	  const EBCellFAB&   Ss        = (*m_Ss[ilev])[datInd];
	  const EBCellFAB&   Phi       = (*m_Phi[ilev])[datInd];

	  const Box& grid = dbl.get(datInd);

	  F_FAB.setVal(0.0);


	  const EBISBox& ebisBox = ebisl[datInd];

	  IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);

	  irregSaturationTerms(F_FAB, psi_FAB, psi_oldFAB, Ss, Phi, 1.0, //rho, 
			       m_Ssat, m_Sres, m_alpha, m_nsoil, 
			       m_eblgsPtr[ilev]->getDomain(), grid);
	  
	  for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      F_FAB(vof,0) *=ebisBox.volFrac(vof); //VolFrac;
	      
	      Real flux = (*m_spreadFlux[ilev])[datInd](vof, 0);
	      const Real& areaFrac = ebisBox.bndryArea(vof);
	      if (areaFrac > 0.0)
	      {
		  //const RealVect anisotNormal = getAnisotropicNormal(ebisBox.normal(vof), a_vectDx[ilev]);
		  flux *= areaFrac;//*anisotNormal[2];
		  
		  Real anisotr_corr=0.0;
		  const RealVect&   normalRef = ebisBox.normal(vof);
		  
		  for (int idir = 0; idir < SpaceDim; idir++)
		  {
		      anisotr_corr += normalRef[idir]*normalRef[idir] 
			  /(a_vectDx[ilev][idir]*a_vectDx[ilev][idir]);
		  }
		  F_FAB(vof,0) -= a_dt*flux*sqrt(anisotr_corr);
	      }
	  }
	}
      a_F_mfree[ilev]->exchange();
    }
}

void 
EBPetscCompGridRichards::getJacobianColoringVectors(std::vector<PetscInt>&  a_rows,
						    std::vector<std::vector<PetscInt> >& a_collumns,
						    std::vector<PetscInt>& a_num_diag,
						    std::vector<PetscInt>& a_num_offdiag,
						    const bool a_include_SUBsurface,
						    const int a_num_comps)
{

    CH_TIME("EBPetscCompGridRichards::getJacobianColoringVector");
    
    // int rank;
    // MPI_Comm_rank(Chombo_MPI::comm, &rank);
    
    const int numGhostCells = 4; //must be consistent with EBPetscMap
    const Interval interval(0, 0);
    const int nvarquad = 1;

    Vector<EBQuadCFStencilColor*>    quadCFStencilColor;
    Vector<LevelData<BaseEBCellFAB<int> >* > procIDGlobal;
    procIDGlobal.resize(m_finestLevel+1,NULL);

    //containers for coarser grid
    Vector<LevelData<BaseEBCellFAB<int64_t> >* > bufferVofGlobalMap;
    Vector<LevelData<BaseEBCellFAB<int> >* > bufferprocIDGlobal;
    if (m_finestLevel > 0)
    {
	quadCFStencilColor.resize(m_finestLevel+1,NULL);
	bufferVofGlobalMap.resize(m_finestLevel,NULL);
	bufferprocIDGlobal.resize(m_finestLevel,NULL);
    }
    
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
        procIDGlobal[ilev] = new LevelData<BaseEBCellFAB<int> >();
	BaseEBCellFactory<int> baseEBFact(m_eblgsPtr[ilev]->getEBISL());
	procIDGlobal[ilev]->define(dbl, 1, numGhostCells*IntVect::Unit, baseEBFact);
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
#pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& datInd = dit[mybox];
	    (*procIDGlobal[ilev])[datInd].setVal(procID());
	}
	procIDGlobal[ilev]->exchange(interval);
	
	if(ilev > 0)
	{
	    //build coloring stencil for CF interpolation 
	    //This is used for coloring jacobian so there are no correct weights in the stencils
	    const DisjointBoxLayout& dblCoar = m_eblgsPtr[ilev]->getDBL();
	    quadCFStencilColor[ilev]  = new  EBQuadCFStencilColor(dbl, dblCoar,
								  m_eblgsPtr[ilev]->getEBISL(), 
								  m_eblgsPtr[ilev-1]->getEBISL(),
								  m_eblgsPtr[ilev-1]->getDomain(),
								  m_refRatios[ilev-1], 
								  nvarquad,
								  *m_eblgsPtr[ilev]->getCFIVS(),
								  m_eblgsPtr[ilev]->getEBIS());
	    
	    const RefCountedPtr<EBCFData>& ebcfData = 	quadCFStencilColor[ilev]->getEBCFData();
    
	    //define coarse ebis layouts buffers
	    bufferVofGlobalMap[ilev-1] = new LevelData<BaseEBCellFAB<int64_t> >();
	    bufferprocIDGlobal[ilev-1] = new LevelData<BaseEBCellFAB<int> >();
	    BaseEBCellFactory<int> coarsenedFact(ebcfData->m_ebislCoarsenedFine);
	    BaseEBCellFactory<int64_t> coarsenedFact64(ebcfData->m_ebislCoarsenedFine);
		
	    bufferVofGlobalMap[ilev-1]->define(ebcfData->m_gridsCoarsenedFine, 1, numGhostCells*IntVect::Unit, coarsenedFact64);
	    bufferprocIDGlobal[ilev-1]->define(ebcfData->m_gridsCoarsenedFine, 1, numGhostCells*IntVect::Unit, coarsenedFact);
	    DataIterator dit = dbl.dataIterator();
	    int nbox=dit.size();
#pragma omp parallel for
	    for (int mybox=0;mybox<nbox; mybox++)
	    {
		const DataIndex& datInd = dit[mybox];
		(*bufferVofGlobalMap[ilev-1])[datInd].setVal(-1);
		(*bufferprocIDGlobal[ilev-1])[datInd].setVal(-1);
	    }
		
	    m_VofGlobalMap[ilev-1]->copyTo(interval, *bufferVofGlobalMap[ilev-1], interval);
	    procIDGlobal[ilev-1]  ->copyTo(interval, *bufferprocIDGlobal[ilev-1], interval);
	    bufferVofGlobalMap[ilev-1]->exchange(interval);
	    bufferprocIDGlobal[ilev-1]->exchange(interval); 
	}
    }
    

    
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	
	if (a_include_SUBsurface)
	{
	    m_DarcyOp[ilev]->defineStencils();
	    
	    //now add saturation terms and partial derivates
	    const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	    DataIterator dit = dbl.dataIterator();
	    int nbox=dit.size();
// #pragma omp parallel for
	    for (int mybox=0;mybox<nbox; mybox++)
	    {
		const DataIndex& datInd = dit[mybox];
		const Box& myBox = dbl.get(datInd);
		// Box grownBox = grow(myBox, 2);
		// grownBox &= m_eblgsPtr[ilev]->getDomain();
		
		//has to correspond to number of ghost cells
		const EBISBox& ebisBox = ebisl[datInd];
		
		VoFIterator& vofit = (*m_vofIterator[ilev])[datInd];
		for(vofit.reset(); vofit.ok(); ++vofit)
		{
		    const VolIndex& vof = vofit();
		    const IntVect& iv = vof.gridIndex();
		    //check flags whether include irregular and covered cells into pets vector
		    if (ebisBox.isIrregular(iv) && !m_includeIrregular) continue;
		    if (ebisBox.isCovered(iv) && !m_includeCovered) continue;
		    const PetscInt row_idx = (*m_VofGlobalMap[ilev])[datInd](vof,0);
		    //pout()<<vof<<"\t indx= "<<row_idx<<endl;
		    if(row_idx<0) continue;
		    
		    //here we try to exclude outside irregular cells with zero area
		    bool skipVoF=false;
		    if (ebisBox.isIrregular(iv) && !ebisBox.bndryArea(vof))
		    {
			skipVoF=true;
			for (int idir=0; idir<3; ++idir)
			{
			    for (SideIterator sit; sit.ok(); ++sit)
			    {
				Vector<FaceIndex> faces = ebisBox.getFaces(vof, idir, sit());
				for (int iface = 0; iface < faces.size(); iface++)
				{
				    const FaceIndex& face = faces[iface];
				    //pout()<<"check "<<vof<<"\t face= "<<face<<"\t side="<<sit()<<"\t bndryArea="<<ebisBox.areaFrac(face)<<endl;
				    if (ebisBox.areaFrac(face)) skipVoF = false;
				    if(!skipVoF) break;
				}
				if(!skipVoF) break;
			    }
			    if(!skipVoF) break;
			}
			// if (skipVoF) pout()<<vof<<"\t indx= "<<row_idx<<"\t bndry="<<ebisBox.bndryArea(vof)<<endl;
		    }
		    if (skipVoF)
		    {
			continue;
		    }
		    
		    VoFStencil  vofStencil;
		    m_DarcyOp[ilev]->getJacobianStencil(vofStencil, vof, datInd);
		    
		    const int numCol = vofStencil.size();
		    std::vector<PetscInt> collumns;
		    int inside=0;
		    int outside=0;
		    
		    if (numCol)
		    {
			//const PetscInt row_idx = m_PetscOffsetIdx+m_VofPetscMap[ilev].at(vof)*a_num_comps+icomp;
			a_rows.push_back(row_idx);
			//inside++;		      
			for (int i = 0; i < numCol; i++)
			{
			    const VolIndex& curVoF = vofStencil.vof(i);
			    
			    const int procid = (*procIDGlobal[ilev])[datInd](curVoF,0);
			    if (procid==procID()) 
				inside++;
			    else 
				outside++;
			    
			    //const Real weight = vofStencil.weight(i);
			    //pout()<<rank<<"-check-0"<<endl;
			    //const PetscInt col_idx = m_PetscOffsetIdx+m_VofPetscMap[ilev].at(curVoF)*a_num_comps+icomp;
			    const PetscInt col_idx =(*m_VofGlobalMap[ilev])[datInd](curVoF,0);
			    if(col_idx>=0) collumns.push_back(col_idx);
			    else
			    {
				//check if it's a ghost cell for finer levels
				if (ilev > 0)
				{
				    if (!myBox.contains(curVoF.gridIndex()))//&&grownBox.contains(curVoF.gridIndex()))
				    {
					getCFInterfaceColoringVectors(collumns, inside, outside,
								      curVoF, datInd, ilev, ebisBox,  procIDGlobal,
								      quadCFStencilColor, bufferVofGlobalMap, bufferprocIDGlobal);
					
				    }
				}
			    }
			    // const IntVect& ivtest = curVoF.gridIndex();
			    // pout()<<numCol<<"\t"<<vof<<"\t"<<curVoF<<"\t"<<row_idx<<"\t"<<col_idx<<"\t"<<ebisBox.isIrregular(ivtest)<<"\t"<<ebisBox.isCovered(iv)<<"\t"
			    // 	    <<weight<<endl;
			}
		    }
		    else
		    {
			// pout()<<"OffStencil \t"<<vof<<"\t indx= "<<row_idx<<endl;
			//const PetscInt row_idx = (*m_VofGlobalMap[ilev])[datInd](vof,0);
			a_rows.push_back(row_idx);		      
			collumns.push_back(row_idx);
			inside++;
		    }
		    
		    
		    // //now flux part 
		    if  (m_includeSurfaceSolver && ebisBox.isIrregular(iv))
		    {
			if (ebisBox.bndryArea(vof)>=Tol)
			{

			    getSurfaceColoringVectors(collumns, inside, outside,
						      vof, datInd, ilev, ebisBox,  procIDGlobal,
						      quadCFStencilColor, bufferVofGlobalMap, bufferprocIDGlobal);
			}
		    }
		    
		    
		    //now EBCFInterpolation part
		    
		    
		    a_collumns.push_back(collumns);
		    a_num_diag.push_back(inside);
		    a_num_offdiag.push_back(outside);
		    //pout()<<"diag="<<inside<<"  offdiag="<<outside<<endl;
		}//vof
	    }//dit
	}
	else
	{
	    const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	    DataIterator dit = dbl.dataIterator();
	    int nbox=dit.size();
// #pragma omp parallel for
	    for (int mybox=0;mybox<nbox; mybox++)
	    {
		const DataIndex& datInd = dit[mybox];
		//has to correspond to number of ghost cells
		const EBISBox& ebisBox = ebisl[datInd];
		
		
		const EBGraph& ebgraph = ebisBox.getEBGraph();
		const Box& grid = dbl.get(datInd);
		IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
		
		//subtract finer level from ivs
		if(ilev  < m_finestLevel)
		{
		    const DisjointBoxLayout& dblFine = m_eblgsPtr[ilev+1]->getDBL();
		    //subtract off finer level regions from ivs
		    const DisjointBoxLayout& finer =  dblFine;
		    for(LayoutIterator lit = finer.layoutIterator(); lit.ok(); ++lit)
		    {
			Box coarsenedFine = finer[lit()];
			coarsenedFine.coarsen(m_refRatios[ilev]);
			ivsIrreg -= coarsenedFine;
		    }
		}
	      
		for(VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
		{
		    const VolIndex& vof = vofit();
		    //const IntVect& iv = vof.gridIndex();
		    
		    const PetscInt row_idx = (*m_VofGlobalMap[ilev])[datInd](vof,0);
		    
		    if(row_idx<0) continue;
		    if (ebisBox.bndryArea(vof)<Tol) continue;
		    std::vector<PetscInt> collumns;
		    int inside=1;
		    int outside=0;
		    
		    a_rows.push_back(row_idx);		      
		    collumns.push_back(row_idx);
		    
		    getSurfaceColoringVectors(collumns, inside, outside,
					      vof, datInd, ilev, ebisBox, procIDGlobal,
					      quadCFStencilColor, bufferVofGlobalMap, bufferprocIDGlobal);
		  
		    a_collumns.push_back(collumns);
		    a_num_diag.push_back(inside);
		    a_num_offdiag.push_back(outside);
		}//vof
	    }//dit
	}
    }//ilev
    
    {
	for (int ilev = 0; ilev <=m_finestLevel; ilev++)
	{
	    if (ilev > 0)
	    {
		delete quadCFStencilColor[ilev];
		delete bufferVofGlobalMap[ilev-1];
		delete bufferprocIDGlobal[ilev-1];
	    }
	    delete procIDGlobal[ilev];
	}
	
    }
}


void 
EBPetscCompGridRichards::getCFInterfaceColoringVectors(std::vector<PetscInt>&                         a_collumns,
						       int&                                           a_inside,
						       int&                                           a_outside,
						       const VolIndex&                                a_vof,
						       const DataIndex&                               a_dit,
						       const int                                      a_lev,
						       const EBISBox&                                 a_ebisBox,
						       const Vector<LevelData<BaseEBCellFAB<int> >*>& a_procIDGlobal,
						       const Vector<EBQuadCFStencilColor*>&           a_quadCFStencilColor,
						       const Vector<LevelData<BaseEBCellFAB<int64_t> >*>& a_bufferVofGlobalMap,
						       const Vector<LevelData<BaseEBCellFAB<int> >*>& a_bufferprocIDGlobal) 
{
    CH_TIME("EBPetscCompGridRichards::getCFInterfaceColoringVectors");    
    if (a_lev > 0)
    {
	VoFStencil fineStencil, coarStencil;
	a_quadCFStencilColor[a_lev]->getCFVoFStencils(a_vof, a_dit, fineStencil, coarStencil);
	const int numFineCol =  fineStencil.size();
	const int numCoarCol =  coarStencil.size();
	if (numFineCol || numCoarCol)
	{
	    for (int i = 0; i < numFineCol; i++)
	    {
		const VolIndex& curVoF = fineStencil.vof(i);
		if (!curVoF.isDefined()) continue;
		const PetscInt col_idx =(*m_VofGlobalMap[a_lev])[a_dit](curVoF,0);
		if(col_idx>=0) 
		{
		    a_collumns.push_back(col_idx);
		    //pout()<<"finer "<<col_idx<<endl;
		    const int procid = (*a_procIDGlobal[a_lev])[a_dit](curVoF,0);
		    if (procid==procID()) 
			a_inside++;
		    else 
			a_outside++;
		}
	    }
	    for (int i = 0; i < numCoarCol; i++)
	    {
		const VolIndex& curVoF = coarStencil.vof(i);
		if (!curVoF.isDefined()) continue;
		const PetscInt col_idx =(*a_bufferVofGlobalMap[a_lev-1])[a_dit](curVoF,0);
		if(col_idx>=0) 
		{
		    a_collumns.push_back(col_idx);
		    //pout()<<"coarser "<<col_idx<<endl;
		    const int procid = (*a_bufferprocIDGlobal[a_lev-1])[a_dit](curVoF,0);
		    if (procid==procID()) 
			a_inside++;
		    else 
			a_outside++;
		}
	    }
	}
    }
}

void 
EBPetscCompGridRichards::getSurfaceColoringVectors(std::vector<PetscInt>&                         a_collumns,
						   int&                                           a_inside,
						   int&                                           a_outside,
						   const VolIndex&                                a_vof,
						   const DataIndex&                               a_dit,
						   const int                                      a_lev,
						   const EBISBox&                                 a_ebisBox,
						   const Vector<LevelData<BaseEBCellFAB<int> >*>& a_procIDGlobal,
						   const Vector<EBQuadCFStencilColor*>&           a_quadCFStencilColor,
						   const Vector<LevelData<BaseEBCellFAB<int64_t> >*>& a_bufferVofGlobalMap,
						   const Vector<LevelData<BaseEBCellFAB<int> >*>& a_bufferprocIDGlobal)
{
    CH_TIME("EBPetscCompGridRichards::getSurfaceColoringVectors");
    Vector<VolIndex> grad2dVoFs;
    
    const DisjointBoxLayout& dbl = m_eblgsPtr[a_lev]->getDBL();
    const Box& myBox = dbl.get(a_dit);
    // Box grownBox = grow(myBox, 2);
    // grownBox &= m_eblgsPtr[a_lev]->getDomain();
    const Vector<UGData>& ugDataVec = (*m_surfaceGridData[a_lev])[a_dit](a_vof,0);
    const int numCol =  ugDataVec.size();
    if (numCol)
    {
	for (int i = -1; i < numCol; i++)
	{
	    VolIndex curVoF;
	    if (i<0) 
	    {
		curVoF = a_vof;
	    }else
	    {
		curVoF = ugDataVec[i].vof;
	    }
	    if (!curVoF.isDefined()) continue;
	    
	    if (m_surface_solver_type==DiffusionWaveUnstructured) grad2dVoFs.push_back(curVoF);
	    
	    const PetscInt col_idx =(*m_VofGlobalMap[a_lev])[a_dit](curVoF,0);
	    if(col_idx>=0) 
	    {
		a_collumns.push_back(col_idx);
		
		const int procid = (*a_procIDGlobal[a_lev])[a_dit](curVoF,0);
		if (procid==procID()) 
		    a_inside++;
		else 
		    a_outside++;
	    }
	    else
	    {
		//check if it's a ghost cell for finer levels
		if (a_lev > 0)
		{
		    if (!myBox.contains(curVoF.gridIndex()))//&&grownBox.contains(curVoF.gridIndex()))
		    {
			getCFInterfaceColoringVectors(a_collumns, a_inside, a_outside,
						      curVoF, a_dit, a_lev, a_ebisBox,  a_procIDGlobal,
						      a_quadCFStencilColor, a_bufferVofGlobalMap, a_bufferprocIDGlobal);
		    }
		}
	    }
	}
    }
    
    //now add 2D grad stencil
    if (m_surface_solver_type==DiffusionWaveUnstructured)
    {
	//grad2dVoFs.push_back(a_vof);
	for (int idir = 0; idir < 2; idir++)
	{

	    VoFStencil& IrrStencil = (*m_2DGradStencil[a_lev])[a_dit](a_vof,idir);
 	    const int numIrrCol =  IrrStencil.size();
	    if (numIrrCol)
	    {
		
		for (int i = 0; i < numIrrCol; i++)
		{
		    const VolIndex& curVoF = IrrStencil.vof(i);
		    grad2dVoFs.push_back(curVoF);
		}
	    }
	}

	const int num2DGrad =  grad2dVoFs.size();
	if (num2DGrad)
	{
	    for (int i = 0; i < num2DGrad; i++)
	    {
		const VolIndex& curVoF = grad2dVoFs[i];
		if (curVoF.isDefined())
		{
		    const PetscInt col_idx =(*m_VofGlobalMap[a_lev])[a_dit](curVoF,0);
		    if(col_idx>=0) 
		    {
			a_collumns.push_back(col_idx);
			
			const int procid = (*a_procIDGlobal[a_lev])[a_dit](curVoF,0);
			if (procid==procID()) 
			    a_inside++;
			else 
			    a_outside++;
		    }
		    else
		    {
			//check if it's a ghost cell for finer levels
			if (a_lev > 0)
			{
			    if (!myBox.contains(curVoF.gridIndex()))
			    {
				getCFInterfaceColoringVectors(a_collumns, a_inside, a_outside,
							      curVoF, a_dit, a_lev, a_ebisBox,  a_procIDGlobal,
							      a_quadCFStencilColor, a_bufferVofGlobalMap, a_bufferprocIDGlobal);
				
			    }
			}
		    }
		}
	    }
	}
    }
}


void 
EBPetscCompGridRichards::getJacobianVectors(std::vector<PetscInt>&  a_rows,
					    std::vector<std::vector<PetscInt> >& a_collumns,
					    std::vector<std::vector<PetscScalar> >& a_data,
					    Vector<LevelData<EBCellFAB>* >& a_psi, 
					    const Real a_dt)
{
    MayDay::Error("Function is not implemented yet");
}



void  
EBPetscCompGridRichards::computeKxKrFaceValues(const Vector<LevelData<EBCellFAB>* >& a_psi)
{
  CH_TIME("EBPetscCompGridRichards::computeKxKrFaceValues");
  //const int nghostTemp = 2;
  Interval interval(0, 0);
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
      const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
      
      //reuse tempKxRho each iteration
      DataIterator dit = dbl.dataIterator();
      int nbox=dit.size();

#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  const Box grid = dbl.get(datInd);
      	  // Box grid = grow(dbl.get(datInd), nghostTemp);
      	  // grid &= m_eblgsPtr[ilev]->getDomain();

	  (*m_KxKrRho[ilev])[datInd].copy(grid, interval, grid, (*m_tempKxRho[ilev])[datInd], interval);
      	  IntVectSet ivs(grid);
      	  const EBGraph& ebgraph = ebisl[datInd].getEBGraph();
      	  for(int idir = 0; idir < SpaceDim; idir++)
      	  {
      	      const EBFaceFAB& kxrhoFAB   = (*m_tempKxRho[ilev])[datInd][idir];
      	      EBFaceFAB& kxkrrhoFAB       = (*m_KxKrRho[ilev])[datInd][idir];
      	      //set face values on domain bc
      	      //iterate over boundary faces
      	      for (SideIterator sit; sit.ok(); ++sit)
      	  	{
      	  	  IntVect ivSideGrid =   grid.sideEnd(sit());
      	  	  IntVect ivSideDom  = m_eblgsPtr[ilev]->getDomain().domainBox().sideEnd(sit());
		  
      	  	  //check to see if we actually are at domain boundary
      	  	  if (ivSideGrid[idir] == ivSideDom[idir])
      	  	    {
      	  	      //create box and ivs adjacent to domain boundary
      	  	      //and interior to domain
      	  	      Box sideBox = adjCellBox(grid, idir, sit(), 1);
      	  	      int ishift = -sign(sit());
      	  	      sideBox.shift(idir, ishift);
      	  	      IntVectSet sideIVS(sideBox);
      	  	      FaceStop::WhichFaces stopCritSide;
      	  	      if (sit() == Side::Lo)
      	  		{
      	  		  stopCritSide = FaceStop::LoBoundaryOnly;
      	  		}
      	  	      else
      	  		{
      	  		  stopCritSide = FaceStop::HiBoundaryOnly;
      	  		}
		      
      	  	      for (FaceIterator faceit(sideIVS, ebgraph, idir, stopCritSide);  faceit.ok(); ++faceit)
      	  	      {
      	  		  const FaceIndex & bndryFace = faceit();
      	  		  kxkrrhoFAB(bndryFace,0) = kxrhoFAB(bndryFace,0);
      	  	      }
      	  	    }
      	  	}//side
      	  }//idir
      }//dit
      //m_KxKrRho[ilev]->exchange();

      
      scaleLevelFluxByKr(*m_KxKrRho[ilev], *a_psi[ilev], *m_DarcyFlux[ilev],
			 m_alpha, m_nsoil, m_lambda,
			 dbl, m_eblgsPtr[ilev]->getEBISL(), m_eblgsPtr[ilev]->getDomain(),
			 m_meanKr_method,
			 ilev);
      
      
      //set bc values
      // DataIterator dit = m_eblgsPtr[ilev]->getDBL().dataIterator();
      // int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  // Box grid = grow(dbl.get(datInd), nghostTemp);
	  // grid &= m_eblgsPtr[ilev]->getDomain();
	  const Box grid = dbl.get(datInd);
	  IntVectSet ivs(grid);
	  const EBGraph& ebgraph = ebisl[datInd].getEBGraph();
	  for(int idir = 0; idir < SpaceDim; idir++)
	    {
		EBFaceFAB& bcoefFAB       = (*m_KxKrRho[ilev])[datInd][idir];
		const EBFaceFAB& kxrhoFAB = (*m_tempKxRho[ilev])[datInd][idir];
		
		const EBCellFAB& psi = (*a_psi[ilev])[datInd];
		//set face values on domain bc

		//iterate over boundary faces
		for (SideIterator sit; sit.ok(); ++sit)
		{
		    IntVect ivSideGrid =   grid.sideEnd(sit());
		    IntVect ivSideDom  = m_eblgsPtr[ilev]->getDomain().domainBox().sideEnd(sit());
		    
		    //check to see if we actually are at domain boundary
		    if (ivSideGrid[idir] == ivSideDom[idir])
		    {
			//create box and ivs adjacent to domain boundary
			//and interior to domain
			Box sideBox = adjCellBox(grid, idir, sit(), 1);
			int ishift = -sign(sit());
			sideBox.shift(idir, ishift);
			IntVectSet sideIVS(sideBox);
			FaceStop::WhichFaces stopCritSide;
		      if (sit() == Side::Lo)
		      {
			  stopCritSide = FaceStop::LoBoundaryOnly;
		      }
		      else
		      {
			  stopCritSide = FaceStop::HiBoundaryOnly;
		      }
		      
		      //do simple extrapolation to boundary faces from neighboring faces in same direction
		      for (FaceIterator faceit(sideIVS, ebgraph, idir, stopCritSide);  faceit.ok(); ++faceit)
		      {
			  const FaceIndex & bndryFace = faceit();
			  
			  const Real localPsi = m_DomainBC->getBCValue(bndryFace,psi, bcoefFAB, sit(),idir,m_vectDxs[ilev], m_time,0);
			  const Real KrBndry  = getRelativePermeability(localPsi, m_alpha, m_nsoil, m_lambda);
			  
			  //Real extrapValue = EBArith::extrapFaceValueToDomain(bndryFace,sit(),a_idir,a_ebGraph,bcoefFAB, 0);
			  bcoefFAB(bndryFace,0) = kxrhoFAB(bndryFace,0)*KrBndry;
			  //pout()<<idir<<"-"<<bndryFace<<"\t"<<localPsi<<"\t"<<bcoefFAB(bndryFace,0)<<"\t"<<(*m_KxKrRho[ilev])[datInd][idir](bndryFace,0)<<endl;
		      }
		    }
		}//side
	    }//idir
      }//dit
      //m_KxKrRho[ilev]->exchange();
  }//ilev
}


void
EBPetscCompGridRichards::computeGradientPsi(const Vector<LevelData<EBCellFAB>* >& a_psi)
{
    CH_TIME("EBPetscCompGridRichards::computeGradientPsi");

    //const int nghostTemp = 1;
    const Interval interval(0,1);
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	//const ProblemDomain& domain  = m_eblgsPtr[ilev]->getDomain();
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
#pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& di = dit[mybox];
	    //has to correspond to number of ghost cells
	    // Box myBox;
	    // if (ilev>0)
	    // {
	    // 	myBox = grow(dbl.get(di), nghostTemp);
	    // 	myBox &= m_eblgsPtr[ilev]->getDomain();
	    // }
	    // else
	    // {
	    // 	myBox= dbl.get(di);
	    // }

	    const Box& myBox= dbl.get(di);

	    const EBISBox&   ebisBox = ebisl[di];
	    IntVectSet ivsIrreg = ebisl[di].getIrregIVS(myBox);
	    const EBCellFAB& psi = (*a_psi[ilev])[di];
	    BaseIVFAB<Real>& grad2DIVF = (*m_2DGrad[ilev])[di];
	    grad2DIVF.setVal(0.);
	    
	    for (VoFIterator vofit(ivsIrreg, ebisl[di].getEBGraph()); vofit.ok(); ++vofit)
	    {
		const VolIndex& vof = vofit();
		const RealVect normal = ebisBox.normal(vof);
		for (int idir=0; idir<2; idir++)
		{
		    Real& grad2D = grad2DIVF(vof,idir);
		    VoFStencil& gradStencil = (*m_2DGradStencil[ilev])[di](vof,idir);
		    const int numIrrCol = gradStencil.size();
		    for (int i = 0; i < numIrrCol; i++)
		    {
			const VolIndex& curVoF = gradStencil.vof(i);
			Real val  = std::max(psi(curVoF, 0), 0.0);
			grad2D -=val*gradStencil.weight(i);
		    }
		}
	    }
	}
	m_2DGrad[ilev]->exchange(interval);
    }
}

void  
EBPetscCompGridRichards::addAMRIrregExchangeFlux(const Vector<LevelData<EBCellFAB>*>& a_psinew,
						 const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_source,
						 const Real a_dt)
{
    
    a_psinew[0]->exchange();
    CH_TIME("EBPetscCompGridRichards::addAMRIrregExchangeFlux");
//    const int nghostTemp = 2;

    Interval interval(0, 0);
    // for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    // {
    // 	const LevelData<EBCellFAB>& psi_LD  = *a_psinew[ilev];
    // 	LevelData<EBCellFAB>& h_LD = *m_tempCC1[ilev];
	
    // 	EBLevelDataOps::axby(h_LD, psi_LD, *m_z[ilev], 1.0, 1.0); //m_gconst/m_rho_zero);
    // 	h_LD.exchange();
    // }

    for (int ilev = 1; ilev <= m_finestLevel; ilev++)
    {
	EBQuadCFInterp& ebinterpolate = *m_quadCFI[ilev];
	
	ebinterpolate.interpolate(*a_psinew[ilev],
				  *a_psinew[ilev-1],
				  interval);
	
	a_psinew[ilev]->exchange();
	
	// ebinterpolate.interpolate(*m_tempCC1[ilev  ],
	// 			  *m_tempCC1[ilev-1],
	// 			  interval);
	
	
	// EBLevelDataOps::axby(*a_psinew[ilev], *m_tempCC1[ilev], *m_z[ilev], 1.0, -1.0);
//	a_psinew[ilev]->exchange();
	
    }

    // for (int ilev = 1; ilev <= m_finestLevel; ilev++)
    // {
    // 	EBQuadCFInterp& ebinterpolate = *m_quadCFI[ilev];
	
    // 	ebinterpolate.interpolate(*a_psinew[ilev  ],
    // 				  *a_psinew[ilev-1],
    // 				  interval);
	
    // 	a_psinew[ilev]->exchange();
    // }
    
    if (m_surface_solver_type == DiffusionWaveUnstructured)
    {
	computeGradientPsi(a_psinew);
    }
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	a_source[ilev]->copyTo(*m_spreadFlux[ilev]);
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
#pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& datInd = dit[mybox];
	    // Box grid;
	    // if (ilev)
	    // {
	    // 	grid = grow(dbl.get(datInd), nghostTemp);
	    // 	grid &= m_eblgsPtr[ilev]->getDomain();

	    // }
	    // else
	    // {
	    // 	grid = dbl.get(datInd);
	    // }
	    const Box& 	grid = dbl.get(datInd);
	    const EBISBox& ebisBox = ebisl[datInd];
	    if(m_surface_solver_type == KinematicWaveUnstructured)
	    {
	    	addExchangeFlux((*m_spreadFlux[ilev])[datInd],
				(*a_psinew[ilev])[datInd],
				(*m_psi_saved[ilev])[datInd],
				(*m_domainSlopes[ilev])[datInd],
				(*m_manningCoeff[ilev])[datInd],
				(*m_surfaceGridData[ilev])[datInd],
				(*m_bndryArea[ilev])[datInd],
				BaseIVFAB<Real>(),
				m_eblgsPtr[ilev]->getDomain(),
				a_dt,
				m_vectDxs[ilev],
				grid,
				ebisBox,
				false,
				m_useBCRegularization,
				m_regConstant);
	    } else if(m_surface_solver_type == DiffusionWaveUnstructured)
	    {
	    	addExchangeFlux((*m_spreadFlux[ilev])[datInd],
				(*a_psinew[ilev])[datInd],
				(*m_psi_saved[ilev])[datInd],
				(*m_domainSlopes[ilev])[datInd],
				(*m_manningCoeff[ilev])[datInd],
				(*m_surfaceGridData[ilev])[datInd],
				(*m_bndryArea[ilev])[datInd],
				(*m_2DGrad[ilev])[datInd],
				m_eblgsPtr[ilev]->getDomain(),
				a_dt,
				m_vectDxs[ilev],
				grid,
				ebisBox,
				true,
				m_useBCRegularization,
				m_regConstant);
	    	
	    } else
	    {
	    	MayDay::Error("Unknown surface solver type");
	    }
		
	}
	m_spreadFlux[ilev]->exchange();
    }
}

void  
EBPetscCompGridRichards::addAMRIrregExchangeFlux(const Vector<LevelData<BaseIVFAB<Real> >*>& a_psinew,
						 const Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_source,
						 const Real a_dt)
{
    CH_TIME("EBPetscCompGridRichards::addAMRIrrefExchangeFlux");
//    const int nghostTemp = 2;

    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	a_psinew[ilev]->exchange();
    }

    copyBIVF2EBCellFAB(a_psinew, m_tempCC1, m_eblgsPtr);

    if (m_surface_solver_type == DiffusionWaveUnstructured)
    {
	computeGradientPsi(m_tempCC1);
    }
    
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	a_source[ilev]->copyTo(*m_spreadFlux[ilev]);
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
#pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& datInd = dit[mybox];
	    // Box grid;
	    // if (ilev)
	    // {
	    // 	grid = grow(dbl.get(datInd), nghostTemp);
	    // 	grid &= m_eblgsPtr[ilev]->getDomain();

	    // }
	    // else
	    // {
	    // 	grid = dbl.get(datInd);
	    // }
	    const Box& 	grid = dbl.get(datInd);

	    const EBISBox& ebisBox = ebisl[datInd];
	    if(m_surface_solver_type == KinematicWaveUnstructured)
	    {
		addExchangeFlux((*m_spreadFlux[ilev])[datInd],
				(*m_tempCC1[ilev])[datInd],
				(*m_psi_saved[ilev])[datInd],
				(*m_domainSlopes[ilev])[datInd],
				(*m_manningCoeff[ilev])[datInd],
				(*m_surfaceGridData[ilev])[datInd],
				(*m_bndryArea[ilev])[datInd],
				BaseIVFAB<Real>(),
				m_eblgsPtr[ilev]->getDomain(),
				a_dt,
				m_vectDxs[ilev],
				grid,
				ebisBox,
				false,
				m_useBCRegularization,
				m_regConstant);

	    }
	    else if(m_surface_solver_type == DiffusionWaveUnstructured)
	    {
		addExchangeFlux((*m_spreadFlux[ilev])[datInd],
				(*m_tempCC1[ilev])[datInd],
				(*m_psi_saved[ilev])[datInd],
				(*m_domainSlopes[ilev])[datInd],
				(*m_manningCoeff[ilev])[datInd],
				(*m_surfaceGridData[ilev])[datInd],
				(*m_bndryArea[ilev])[datInd],
				(*m_2DGrad[ilev])[datInd],
				m_eblgsPtr[ilev]->getDomain(),
				a_dt,
				m_vectDxs[ilev],
				grid,
				ebisBox,
				true,
				m_useBCRegularization,
				m_regConstant);
	    }else
	    {
		MayDay::Error("Unknown surface solver type");
	    }
	}
	m_spreadFlux[ilev]->exchange();
    }
}
void  
EBPetscCompGridRichards::setZ(Vector<LevelData<EBCellFAB>* >& a_z)
{
    CH_TIME("EBPetscCompGridRichards::setZ");
    const int nghostTemp = 2;
    Interval interv(0, 0);
    
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
#pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& datInd = dit[mybox];
	    //has to correspond to number of ghost cells
	    Box grownBox;
	    if (ilev>0)
	    {
	    	grownBox = grow(dbl.get(datInd), nghostTemp);
	    	grownBox &= m_eblgsPtr[ilev]->getDomain();
	    }
	    else
	    {
	    	grownBox = dbl.get(datInd);
	    }
	    //const Box& 	grownBox = dbl.get(datInd);
	    
	    const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	    const EBISBox& ebisBox = ebisl[datInd];
	    IntVectSet ivsTot(grownBox);
	    //IntVectSet ivsIrreg = ebisl[datInd].getIrregIVS(grownBox);
	    for (VoFIterator vofit(ivsTot, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
		const VolIndex& vof = vofit();
		
		Real Z = 0.5;
		const IntVect& iv = vof.gridIndex();
		Z += iv[2];
		if (ebisBox.isIrregular(iv) && !m_irregular_value_cc)
		{
		    Z += ebisBox.centroid(vof)[2];
		    // Z += ebisBox.bndryCentroid(vof)[2];
		}
		Z *= m_vectDxs[ilev][2];
		
		//Z += a_probLo;
		(*a_z[ilev])[datInd](vof,0) = Z;
	    }
	}
      a_z[ilev]->exchange(interv);
    }
}

void  
EBPetscCompGridRichards::computeSaturation(LevelData<EBCellFAB>& a_Sw, 
					   const LevelData<EBCellFAB>& a_psi,
					   const ProblemDomain& a_domain)
{
    CH_TIME("EBPetscCompGridRichards::computeSaturation");
    //const int nghostTemp = 1;

    EBLevelDataOps::setToZero(a_Sw);
    const DisjointBoxLayout &dbl = a_psi.getBoxes();
    DataIterator dit = a_psi.dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
    for (int mybox=0;mybox<nbox; mybox++)
    {
	const DataIndex& datInd = dit[mybox];
	EBCellFAB&         Sw_FAB    =  a_Sw[datInd];
	const EBCellFAB&   psi_FAB     =  a_psi[datInd];

	// Box  grid = grow(dbl.get(datInd), nghostTemp);
	// grid &= a_domain;
	const Box& grid = dbl.get(datInd);
	
	getSaturation(Sw_FAB, psi_FAB,  
		      m_Ssat, m_Sres, m_alpha, m_nsoil, 
		      a_domain, grid);
    }
    a_Sw.exchange();
}

void  
EBPetscCompGridRichards::fillRandom(Vector<LevelData<EBCellFAB>* >& a_psi)
{
    CH_TIME("EBPetscCompGridRichards::fillRandom");
  //const int nghostTemp = 2;
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
      const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
      DataIterator dit = m_eblgsPtr[ilev]->getDBL().dataIterator();
      int nbox=dit.size();
//#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  //has to correspond to number of ghost cells
	  const Box& grid = (*a_psi[ilev])[datInd].getRegion();
	  // Box grownBox = grow(dbl.get(datInd), nghostTemp);
	  // grownBox &= m_eblgsPtr[ilev]->getDomain();
	  IntVectSet ivsTot(grid);
   	  for(VoFIterator vofit(ivsTot,  ebisl[datInd].getEBGraph()); vofit.ok(); ++vofit)
	    {
	      const VolIndex& vof = vofit();
	      (*a_psi[ilev])[datInd](vof, 0) = rand();
	    }
	}
      a_psi[ilev]->exchange();
    }
}

void  
EBPetscCompGridRichards::getDivDarcyFlux(Vector<LevelData<EBCellFAB>* >& a_Div,
					 const Vector<LevelData<EBCellFAB>* >& a_psi)
{
    CH_TIME("EBPetscCompGridRichards::getDivDarcyFlux");
  //set harmonic averages KxKr on faces for conductive operator as b_coeff 
  computeKxKrFaceValues(a_psi);

  EBAMRDataOps::setToZero(m_tempCC1);

 //const int nghostTemp = 2;
  
  // no need since we set explicitely ebbc flux //set irregular Kxkr 

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      const LevelData<EBCellFAB>& psi_LD  = *a_psi[ilev];
      LevelData<EBCellFAB>& h_LD = *m_tempCC1[ilev];
      
      EBLevelDataOps::axby(h_LD, psi_LD, *m_z[ilev], 1.0, 1.0);
      h_LD.exchange();
  }
  //EBAMRDataOps::setVal(m_tempCC1, 1.0);
  // m_tempCC1[0]->exchange();
  // Interval interval(0, 0);
  // for (int ilev = 1; ilev <= m_finestLevel; ilev++)
  // {
  //     EBQuadCFInterp& ebinterpolate = *m_quadCFI[ilev];
      
      
  //     ebinterpolate.interpolate(*m_tempCC1[ilev  ],
  // 				*m_tempCC1[ilev-1],
  // 				interval);
  //     m_tempCC1[ilev]->exchange();
  // }

  LevelData<EBCellFAB> dummyLDdata;

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      LevelData<EBCellFAB>& F_LD  = *a_Div[ilev];
      EBLevelDataOps::setToZero(F_LD);

      const LevelData<EBCellFAB>& h_LD = *m_tempCC1[ilev];
      LevelData<EBCellFAB> *hfineLDPtr = &dummyLDdata;
      LevelData<EBCellFAB> *hcoarLDPtr = &dummyLDdata;
      EBDarcyOp* finerOp=NULL;
      if (ilev < m_finestLevel) 
      {
	 hfineLDPtr = m_tempCC1[ilev+1];
	 finerOp = m_DarcyOp[ilev+1];
      }
      
      if (ilev > 0) 
      {
	 hcoarLDPtr = m_tempCC1[ilev-1];
      }

      m_DarcyOp[ilev]->defineStencils();
      m_DarcyOp[ilev]->applyOp(F_LD, h_LD, false);
      // m_DarcyOp[ilev]->applyAMROp(F_LD, *hfineLDPtr, h_LD, *hcoarLDPtr, false, finerOp);


      F_LD.exchange();
    }
}

//depricated function for computing Darcy fluxes 
void  
EBPetscCompGridRichards::getCCDarcyFlux(Vector<LevelData<EBCellFAB>* >& a_CCflux,
					const Vector<LevelData<EBCellFAB>* >& a_psi)
{
    CH_TIME("EBPetscCompGridRichards::getCCDarcyFlux");
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();

      //use saved darcy fluxes
      ccpLevelccpAverageFaceToCells(*a_CCflux[ilev],
				    *m_DarcyFlux[ilev],
				    dbl,
				    m_eblgsPtr[ilev]->getEBISL(),
				    m_eblgsPtr[ilev]->getDomain(),
				    m_vectDxs[ilev],
				    ilev);
      
      //add eb flux    
      m_DarcyOp[ilev]->defineStencils();
      DataIterator dit = m_eblgsPtr[ilev]->getDBL().dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  m_DarcyOp[ilev]->getFluxIrregular((*a_CCflux[ilev])[datInd],
					    *a_psi[ilev],
					    datInd);
      }
  }
}


void
EBPetscCompGridRichards::updateParameters(const Real a_time)
{
    CH_TIME("EBPetscCompGridRichards::updateParameters");
  m_time = a_time;

  ParserFunc SsFunc(m_paramExpr.SpecStorage);
  ParserFunc KxFunc_0(m_paramExpr.SatConductivity[0]);
  ParserFunc KxFunc_1(m_paramExpr.SatConductivity[1]);
  ParserFunc KxFunc_2(m_paramExpr.SatConductivity[2]);
  ParserFunc PhiFunc(m_paramExpr.Porosity);
  ParserFunc MannCoeffFunc;
  if (m_includeSurfaceSolver) MannCoeffFunc.define(m_paramExpr.ManningCoeff);
  //  ParserFunc RhoFunc(m_paramExpr.Density);
 
  //const int nghostTemp = 2;
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
      const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
      DataIterator dit = dbl.dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  const EBISBox& ebisBox = ebisl[datInd];
	  //has to correspond to number of ghost cells
	  // Box grownBox;
	  // if (ilev>0)
	  // {
	  //     grownBox = grow(dbl.get(datInd), nghostTemp);
	  //     grownBox &= m_eblgsPtr[ilev]->getDomain();
	  // }
	  // else
	  // {
	  //     grownBox = dbl.get(datInd);
	  // }
	  const Box& grownBox = dbl.get(datInd);

	  IntVectSet ivsTot(grownBox);
	  // VoFIterator& vofit = (*m_vofIterator[ilev])[datInd];
	  // for(vofit.reset(); vofit.ok(); ++vofit)
   	  for(VoFIterator vofit(ivsTot,  ebisl[datInd].getEBGraph()); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      const IntVect& iv = vof.gridIndex();
	      RealVect point = 0.5*RealVect::Unit;
	      point += iv;
	      if (!m_irregular_value_cc && ebisBox.isIrregular(iv)) point += ebisBox.bndryCentroid(vof);//point += ebisBox.centroid(vof);
	      for (int idir = 0; idir < SpaceDim; idir++)
	      {
		  point[idir] *= m_vectDxs[ilev][idir];
	      }
	      Real depth;
	      if (ebisBox.isIrregular(iv))
	      {
		  // RealVect distVec = ebisBox.bndryCentroid(vof);
		  // distVec *=m_vectDx[ilev];
		  
		  // depth  = -distVec.dotProduct(getAnisotropicNormal(ebisBox.normal(vof), m_vectDx[ilev]));
		  if (ebisBox.bndryArea(vof)>1.e-12)
		      depth = 0.0;
		  else
		      depth = -m_implicitBaseIF->value(point/m_vectDxs[0])*m_vectDxs[0][2];
	      }
	      else
	      {
		  depth = -m_implicitBaseIF->value(point/m_vectDxs[0])*m_vectDxs[0][2];
	      }
	      if (depth < 0.0) depth = 0.0;
	      
	      ParserInput parserinput;
	      parserinput.point = point;
	      parserinput.time  = m_time;
	      parserinput.depth = depth;
	      //parserinput.psi   = (*m_psi_current[ilev])[datInd](vof,0);
	      //parserinput.slope = ;
	      
	      //(*m_Rho[ilev])[datInd](vof,0) = RhoFunc.Eval(point, a_time);
	      //parserinput.rho   = m_richardsCompGrid->m_Pho[ilev][datInd](vof,0);
	      
	      //update phi before parser other functions
	      (*m_Phi[ilev])[datInd](vof,0) = PhiFunc.Eval(parserinput);
	      parserinput.phi   = (*m_Phi[ilev])[datInd](vof,0);
	      
	      if (m_includeSurfaceSolver && ebisBox.isIrregular(iv))
	      {
		  //get default slopes (not sure how to deal here with underresolved cells)
		  parserinput.slope_x = (*m_domainSlopes[ilev])[datInd](vof,0);
		  parserinput.slope_y = (*m_domainSlopes[ilev])[datInd](vof,1);
		  //update Manning Coeff before parser other functions
		  (*m_manningCoeff[ilev])[datInd](vof,0) = MannCoeffFunc.Eval(parserinput);
		  parserinput.manncoef =(*m_manningCoeff[ilev])[datInd](vof,0);
	      }
		
	      (*m_Ss[ilev])[datInd](vof,0) = SsFunc.Eval(parserinput);
	      (*m_Kx[ilev])[datInd](vof,0) = KxFunc_0.Eval(parserinput);
	      (*m_Kx[ilev])[datInd](vof,1) = KxFunc_1.Eval(parserinput);
	      (*m_Kx[ilev])[datInd](vof,2) = KxFunc_2.Eval(parserinput);

	    }
	}
      
      m_Ss[ilev]->exchange();
      m_Kx[ilev]->exchange();
      m_Phi[ilev]->exchange();
      if (m_includeSurfaceSolver)  m_manningCoeff[ilev]->exchange();
      
      ccpHarmonicAverageToFaces(*m_tempKxRho[ilev],*m_Kx[ilev],  
      				dbl, m_eblgsPtr[ilev]->getEBISL(), m_eblgsPtr[ilev]->getDomain(),
				ilev);

      //set bc values for Kx by direct computation
      // DataIterator dit = m_eblgsPtr[ilev]->getDBL().dataIterator();
      // int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  // Box grownBox;
	  // if (ilev>0)
	  // {
	  //     grownBox = grow(dbl.get(datInd), nghostTemp);
	  //     grownBox &= m_eblgsPtr[ilev]->getDomain();
	  // }
	  // else
	  // {
	  //     grownBox = dbl.get(datInd);
	  // }
	  const Box& grownBox = dbl.get(datInd);

	  IntVectSet ivs(grownBox);
	  const EBISBox& ebisBox = ebisl[datInd];
	  const EBGraph& ebgraph = ebisBox.getEBGraph();
	  for(int idir = 0; idir < SpaceDim; idir++)
	    {
	      EBFaceFAB& kxrhoFAB = (*m_tempKxRho[ilev])[datInd][idir];
	      //set face values on domain bc
	      //iterate over boundary faces
	      for (SideIterator sit; sit.ok(); ++sit)
	      {
		  IntVect ivSideGrid =   grownBox.sideEnd(sit());
		  IntVect ivSideDom  = m_eblgsPtr[ilev]->getDomain().domainBox().sideEnd(sit());
		  
		  //check to see if we actually are at domain boundary
		  if (ivSideGrid[idir] == ivSideDom[idir])
		  {
		      //create box and ivs adjacent to domain boundary
		      //and interior to domain
		      Box sideBox = adjCellBox(grownBox, idir, sit(), 1);
		      int ishift = -sign(sit());
		      sideBox.shift(idir, ishift);
		      IntVectSet sideIVS(sideBox);
		      FaceStop::WhichFaces stopCritSide;
		      if (sit() == Side::Lo)
		      {
			  stopCritSide = FaceStop::LoBoundaryOnly;
		      }
		      else
		      {
			  stopCritSide = FaceStop::HiBoundaryOnly;
		      }
		      
		      //do simple extrapolation to boundary faces from neighboring faces in same direction
		      for (FaceIterator faceit(sideIVS, ebgraph, idir, stopCritSide);  faceit.ok(); ++faceit)
		      {
			  const FaceIndex & bndryFace = faceit();
			  RealVect point = EBArith::getFaceLocation(bndryFace, m_vectDxs[ilev], RealVect::Zero);
			  
			  ParserFunc KxFunc(m_paramExpr.SatConductivity[idir]);
			  
			  Real depth;
			  depth = -m_implicitBaseIF->value(point/m_vectDxs[0])*m_vectDxs[0][2];
			  if (depth < 0.0) depth = 0.0;
			  
			  //ParserFunc RhoFunc(m_paramExpr.Density);
			  ParserInput parserinput;
			  parserinput.point = point;
			  parserinput.time  = m_time;
			  parserinput.depth = depth;
			  parserinput.phi   = PhiFunc.Eval(parserinput);
			  
			  Real kxrhoBndry = KxFunc.Eval(parserinput);
			  
			  //pout()<<bndryFace<<"\t"<<point<<"\t"<<localPsi<<"\t"<<kxrhoBndry<<endl;
			  //Real extrapValue = EBArith::extrapFaceValueToDomain(bndryFace,sit(),a_idir,a_ebGraph,a_faceData,0);
			  kxrhoFAB(bndryFace,0) = kxrhoBndry;
		      }
		  }
	      }//side
	    }//idir
      }//dit
      //m_tempKxRho[ilev]->exchange();
  }//ilev
}


void
EBPetscCompGridRichards::setDomainSlopes()
{
    CH_TIME("EBPetscCompGridRichards::setDomainSlopes");
    //const int nghostTemp = 2;
    const Interval interval(0,1);
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	//const ProblemDomain& domain  = m_eblgsPtr[ilev]->getDomain();
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	const RealVect& vectDx = m_vectDxs[ilev];
	
	const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
#pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& di = dit[mybox];
	    //has to correspond to number of ghost cells
	    // Box myBox;
	    // if (ilev>0)
	    // {
	    // 	myBox = grow(dbl.get(di), nghostTemp);
	    // 	myBox &= m_eblgsPtr[ilev]->getDomain();
	    // }
	    // else
	    // {
	    // 	myBox= dbl.get(di);
	    // }
	    const Box& 	myBox= dbl.get(di);
	    
	    const EBISBox&   ebisBox = ebisl[di];
	    IntVectSet ivsIrreg = ebisl[di].getIrregIVS(myBox);
	    BaseIVFAB<Real>& slope = (*m_domainSlopes[ilev])[di];
	    
	    for (VoFIterator vofit(ivsIrreg, ebisl[di].getEBGraph()); vofit.ok(); ++vofit)
	    {
		const VolIndex& vof = vofit();
		const RealVect normal = ebisBox.normal(vof);
		for (int idir=0; idir<2; idir++)
		{
		    if (m_paramExpr.frictionSlopes[idir] == "auto")
		    {
			if (fabs(normal[2])>1.e-12)
			{
			    slope(vof, idir) = -normal[idir]*vectDx[2]/(normal[2]*vectDx[idir]);
			}
			else
			{
			    slope(vof, 0) = 0.0;
			    slope(vof, 1) = 0.0;
			}
		    }else
		    {
			ParserFunc slopeFunc(m_paramExpr.frictionSlopes[idir]);
			const IntVect& iv = vof.gridIndex();
			RealVect point = 0.5*RealVect::Unit;
			point += iv;
			point += ebisBox.bndryCentroid(vof);
			//point += ebisBox.centroid(vof);
			for (int iidir = 0; iidir < SpaceDim; iidir++)
			{
			    point[iidir] *= m_vectDxs[ilev][iidir];
			}

			
			ParserInput parserinput;
			parserinput.point = point;
			parserinput.time  = m_time;
			parserinput.depth = 0.0; //surface cells have zero depth
			slope(vof, idir) = slopeFunc.Eval(parserinput);
		    }
		}
	    }
	    
	}
	m_domainSlopes[ilev]->exchange(interval);
    }
}

//#define UGDEBUGOUT
//std::ofstream a_dumpFile;
void
EBPetscCompGridRichards::defineUnstructuredSurfaceGrid()
{
    CH_TIME("EBPetscCompGridRichards::defineUnstructuredSurfaceGrid");
    //a_dumpFile.open("dataEB.dat", std::ofstream::out);
   
//    const int nghostTemp = 1;

    //create IFData2 class to reuse it in findIntersection methods 
    IFData2<3> ifdata(*m_implicitBaseIF);
     //gets VofStencil for Irr surface
    const ProblemDomain& domain0  = m_eblgsPtr[0]->getDomain();
//    const RealVect& vectDx0 = m_vectDxs[0];

    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	int numUnderresolved=0;
	const ProblemDomain& domain  = m_eblgsPtr[ilev]->getDomain();
	//const IntVect domainSize = domain.size();
	const RealVect& vectDx = m_vectDxs[ilev];
	const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
//// #pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& datInd = dit[mybox];
	    const EBISBox&   ebisBox = ebisl[datInd];
	    (*m_bndryArea[ilev])[datInd].setVal(0.0);

	    // Box grid;
	    // if (ilev)
	    // {
	    // 	grid = grow(dbl.get(datInd), nghostTemp);
	    // 	grid &= domain;
	    // }
	    // else
	    // {
	    // 	grid = dbl.get(datInd);
	    // 
	    
	    const Box& grid = dbl.get(datInd);
	    
	    IntVectSet irregSet =  ebisBox.getIrregIVS(grid);
	    IntVectSet ivsTot(grid);

	    //pout()<<" containsTot_before="<<ivsTot.contains(IntVect(52,95,74))<<" containsIrr_before="<<irregSet.contains(IntVect(52,95,74));

	    //subtract off finer level regions from ivs
	    if(ilev  < m_maxLevel)
	    {
	    	//subtract off finer level regions from ivs
	    	const DisjointBoxLayout& finer =  m_eblgsPtr[ilev+1]->getDBL();
	    	for(LayoutIterator lit = finer.layoutIterator(); lit.ok(); ++lit)
	    	{
	    	    Box coarsenedFine = finer[lit()];
	    	    coarsenedFine.coarsen(m_refRatios[ilev]);
	    	    irregSet -= coarsenedFine;
	    	}
	    }

	    //pout()<<" containsCoarsened="<<containsCoarsened<<" containsIrr_after="<<irregSet.contains(IntVect(52,95,74))<<endl;

	    for (VoFIterator vofit(irregSet,ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
		const VolIndex& vof = vofit();
		
		if (ebisBox.bndryArea(vof)<Tol) continue;
		const IntVect& vofIV = vof.gridIndex();
		//pout()<<vof<<" ";
		// if(vofIV==IntVect(52,95,74)) 
		// {
		//     pout()<<" contains_after="<<irregSet.contains(vofIV)<<endl;
		//     continue; //this is a weird thing, some vofs are in the iterator
		// }

		const RealVect cellCorner = RealVect(vofIV);
		Vector<UGData> vofStencil2D;

		Vector<RealVect> points;
		Vector<IntVect>  sides;
		Vector<PairPoints> pairPoints;
		Vector<int> savedIndex;
#ifdef UGDEBUGOUT
		pout()<< "*** New VoF "<< vof<<" ***************************************************** "<<endl;
		pout()<< "Level = "<<ilev<<"\t number of vofs_in_cell= "<< ebisBox.numVoFs(vofIV)<<endl;
		pout()<< "isIrregular="<<ebisBox.isIrregular(vofIV)
		      <<"\t isCovered="<<ebisBox.isCovered(vofIV)
		      <<"\t isMultiValued="<<ebisBox.isMultiValued (vofIV)
		      <<endl;
#endif		    
		
		const int numPoints=findIntersectionPoints(ifdata, points, sides,vof, ilev, m_refRatios);
		
		// if (vof.gridIndex()==IntVect(40,43,70))
		// {
		    
		//     int size=points.size();
		//     for (int ipoint=0;ipoint<size;ipoint++)
		//     {
		// 	pout()<< "Edges:  "<<points[ipoint]<<"\t"<<sides[ipoint]<<endl;
		//     }
		// }
#ifdef UGDEBUGOUT
		int size=points.size();
		for (int ipoint=0;ipoint<size;ipoint++)
		{
			pout()<< "Edges:  "<<points[ipoint]<<"\t"<<sides[ipoint]<<endl;
		}
#endif		    
		const int numPairs=findPairPoints(pairPoints,points,sides,vofIV);

#ifdef UGDEBUGOUT
		const int sizePair = pairPoints.size();
		pout()<<"numPoints="<<numPoints<<"\t numPairs="<<numPairs<<endl;
		for (int ipoint=0;ipoint<sizePair;ipoint++)
		{
		    pout()<< "Pair:  first= "<<pairPoints[ipoint].first<<"\t second= "<<pairPoints[ipoint].second;
		    const int size2 = pairPoints[ipoint].alternIVs.size();
		    for (int ii=0;ii<size2;ii++)
		    {
			pout()<< "\t alternIV="<<pairPoints[ipoint].alternIVs[ii];
		    }
		    pout()<<endl;
		}
#endif
		//set UGData stencil 
		setUGStencil(vofStencil2D,  savedIndex, ifdata, pairPoints, vof, ilev, ebisBox, domain, domain0,vectDx);

#ifdef UGDEBUGOUT
		pout()<<"numPoints="<<numPoints
		      <<"\t numPairs="<<numPairs
		      <<"\t size="<<vofStencil2D.size()				
		      <<"\t bndryAreaEbis="<<ebisBox.bndryArea(vof)
		      <<"\t normalEbis="<<ebisBox.normal(vof)
		      <<"\t volFrac="<<ebisBox.volFrac(vof)
		      <<endl;
#endif
		CH_assert(numPairs==vofStencil2D.size());

		//check for underresolved eb cells
		if (numPairs!=numPoints)
		{
		    // if (vofIV==IntVect(21,31,274))
		    // {
		    // 	pout()<<"=> vof="<<vof
		    // 	      <<"\t bndryArea="<<(*m_bndryArea[ilev])[datInd](vof,0)
		    // 	      <<endl;
		    // 	for (int is=0;is<vofStencil2D.size();is++)
		    // 	{
		    // 	    //pout()<<scientific;
		    // 	    pout()<<"\t conVoF="<<vofStencil2D[is].vof
		    // 		  <<"\t normal="<<vofStencil2D[is].faceData
		    // 		  <<"\t bndryDir="<<vofStencil2D[is].bndryDir
		    // 		  <<endl;
		    // 	}
		    // 	const int sizePair = pairPoints.size();
		    // 	pout()<<"numPoints="<<numPoints<<"\t numPairs="<<numPairs<<endl;
		    // 	for (int ipoint=0;ipoint<sizePair;ipoint++)
		    // 	{
		    // 	    pout()<< "Pair:  first= "<<pairPoints[ipoint].first<<"\t second= "<<pairPoints[ipoint].second;
		    // 	    const int size2 = pairPoints[ipoint].alternIVs.size();
		    // 	    for (int ii=0;ii<size2;ii++)
		    // 	    {
		    // 		pout()<< "\t alternIV="<<pairPoints[ipoint].alternIVs[ii];
		    // 	    }
		    // 	    pout()<<endl;
		    // 	}
		    // 	pout()<<endl;
		    // }

		    setUnderresolvedUGStencil(vofStencil2D,
					      pairPoints,
					      (*m_bndryArea[ilev])[datInd],
					      ifdata,
					      savedIndex,
					      vof,
					      ilev,
					      ebisBox, domain, domain0, vectDx);
		    (*m_surfaceGridData[ilev])[datInd](vof,0)=vofStencil2D;
		    pout()<<"Underresolved cell "<<vof<<endl;
		    numUnderresolved++;
 		}//underresolved if
		else
		{
		    (*m_surfaceGridData[ilev])[datInd](vof,0)=vofStencil2D;
		    const Real ebArea = computeShoeLaceArea(points, m_vectDxs[0]);//vectDx);
		    (*m_bndryArea[ilev])[datInd](vof,0)=ebArea;

		    if (!numPoints)
		    {
			pout()<< "This is a big problem- zero number of points for the VoF "<< vof<<" ! "<<endl;
			pout()<< "Level = "<<ilev<<"\t number of vofs_in_cell= "<< ebisBox.numVoFs(vofIV)<<endl;
			pout()<< "isIrregular="<<ebisBox.isIrregular(vofIV)
			      <<"\t isCovered="<<ebisBox.isCovered(vofIV)
			      <<"\t isMultiValued="<<ebisBox.isMultiValued (vofIV)
			      <<"\t bndryAreaEbis="<<ebisBox.bndryArea(vof)
			      <<"\t normalEbis="<<ebisBox.normal(vof)
			      <<"\t volFrac="<<ebisBox.volFrac(vof)
			      <<endl;
			MayDay::Abort("numPairs is zero!");
		    }
		    
#ifdef UGDEBUGOUT
		    pout()<<"=> vof="<<vof
			  <<"\t bndryArea="<<(*m_bndryArea[ilev])[datInd](vof,0)
			  <<endl;
		    for (int is=0;is<vofStencil2D.size();is++)
		    {
		        //pout()<<scientific;
			pout()<<"\t conVoF="<<vofStencil2D[is].vof
			      <<"\t normal="<<vofStencil2D[is].faceData
			      <<"\t bndryDir="<<vofStencil2D[is].bndryDir
			      <<endl;
		    }
		     pout()<<endl;
#endif
		}
	    }//vofIterator
	}//boxIterator
	m_bndryArea[ilev]->exchange();
	int totNumUnderresolved;
	MPI_Reduce(&numUnderresolved, &totNumUnderresolved, 1, MPI_INT, MPI_SUM, 0, Chombo_MPI::comm);
	if(!procID()) pout()<<"Total number of underresolved cells on level "<< ilev<<" is "<< totNumUnderresolved<<endl;
    }

    //a_dumpFile.close();

    if (m_surface_solver_type==DiffusionWaveUnstructured)
    {
	define2DGradStencil(0, m_finestLevel);
    }

#if 0
  if (m_finestLevel==m_maxLevel)
  {
      pout()<<endl<<endl;
      int count=0;
      int ilev = m_finestLevel;
      //const ProblemDomain& domain  = m_eblgsPtr[ilev]->getDomain();
      const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
      const RealVect& vectDx = m_vectDxs[ilev];
      
      const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
      DataIterator dit = dbl.dataIterator();
      int nbox=dit.size();
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& di = dit[mybox];
	  //has to correspond to number of ghost cells
	  // Box myBox;
	  // if (ilev>0)
	  // {
	  //     myBox = grow(dbl.get(di),  a_numGhostCells);
	  //     myBox &= m_eblgsPtr[ilev]->getDomain();
	  // }
	  // else
	  // {
	  //     myBox= dbl.get(di);
	  // }
	  const EBISBox&   ebisBox = ebisl[di];
	  IntVectSet ivsIrreg = ebisl[di].getIrregIVS(myBox);
	  const BaseIVFAB<Vector<UGData> >&  gridData =  (*m_surfaceGridData[ilev])[di];
	  
	  for (VoFIterator vofit(ivsIrreg, ebisl[di].getEBGraph()); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      const Vector<UGData>&  data = gridData(vof,0);
	      count++;
	      for (int i=0; i<data.size(); i++)
	      {
		  pout()<<vof<<"\t"
			<<data[i].vof<<"\t"
			<<data[i].faceData<<"\t"
			<<data[i].slopes<<"\t"
			<<data[i].bndryDir<<"\t"
			<<(*m_bndryArea[ilev])[di](vof,0)
			<<endl;
	      }
	  }
      }
      pout()<<"total count"<<count<<endl;
  }
#endif
}

void  
EBPetscCompGridRichards:: setUGStencil(Vector<UGData>&           a_vofStencil2D,
				       Vector<int>&              a_savedIndex,
				       IFData2<3>&               a_ifdata,
				       const Vector<PairPoints>& a_pairPoints,
				       const VolIndex&           a_vof,
				       const int                 a_lev,
				       const EBISBox&            a_ebisBox,
				       const ProblemDomain&      a_domain,
				       const ProblemDomain&      a_domainZero,
				       const RealVect&           a_vectDx)
{
    CH_TIME("EBPetscCompGridRichards::setUGStencil");
    const int sizePair = a_pairPoints.size();
    const IntVect domainSize = a_domainZero.size();
    const IntVect& vofIV = a_vof.gridIndex();

    //IFData2<3> ifdata(*m_implicitBaseIF);

    for (int ipair=0;ipair<sizePair;ipair++)
    {
	
#ifdef UGDEBUGOUT
	pout()<<"   checking "<<a_pairPoints[ipair].alternIVs[0]<<endl;
#endif		    
	bool isBoundary = false;
//		    if (pairPoints[ipair].iv[0]<0 || pairPoints[ipair].iv[0] > domainSize[0]-1)
	if ((a_pairPoints[ipair].first[0]<=0 && a_pairPoints[ipair].second[0] <=0)
	    || (a_pairPoints[ipair].first[0] >= domainSize[0] && a_pairPoints[ipair].second[0] >= domainSize[0]))
	{
	    isBoundary = true;
	    if (a_pairPoints[ipair].first[1]-a_pairPoints[ipair].second[1])
	    {
#ifdef UGDEBUGOUT
		pout()<<"\t connected Pair \t"<<a_vof<<"\t bndryDir= "<< ((a_pairPoints[ipair].alternIVs[0][0]<0)?0:1) <<endl;
#endif
		UGData face_2d;
		const bool correctFace = setUGData(face_2d,
						   a_pairPoints[ipair],
						   a_vof,
						   VolIndex(),
						   a_ebisBox,
						   a_vectDx,
						   ((a_pairPoints[ipair].alternIVs[0][0]<0)?0:1));
		
		if (correctFace && !isUGDataInStencil(face_2d, a_vofStencil2D))
		{
		    a_vofStencil2D.push_back(face_2d);
		    a_savedIndex.push_back(ipair);
		}
	    }
	}
	if ((a_pairPoints[ipair].first[1]<=0 && a_pairPoints[ipair].second[1]<=0)
	    ||(a_pairPoints[ipair].first[1] >= domainSize[1] && a_pairPoints[ipair].second[1] >= domainSize[1]))
	{
	    isBoundary = true;
	    if (a_pairPoints[ipair].first[0]-a_pairPoints[ipair].second[0])
	    {
		
#ifdef UGDEBUGOUT
		pout()<<"\t connected Pair \t"<<a_vof<<"\t bndryDir= "<<((a_pairPoints[ipair].alternIVs[0][1]<0)?2:3)<<endl;
#endif
		UGData face_2d;
		const bool correctFace=setUGData(face_2d,
						 a_pairPoints[ipair],
						 a_vof,
						 VolIndex(),
						 a_ebisBox,
						 a_vectDx,
						 ((a_pairPoints[ipair].alternIVs[0][1]<0)?2:3));
		
		if (correctFace && !isUGDataInStencil(face_2d, a_vofStencil2D))
		{
		    a_vofStencil2D.push_back(face_2d);
		    a_savedIndex.push_back(ipair);
		}
		
	    }
	}
		    
	//special case of crossing boundary in z direction
	//no boundry BC - it will take care in Richards solver
	if (a_pairPoints[ipair].alternIVs[0][2]<0 || a_pairPoints[ipair].alternIVs[0][2] > a_domain.size()[2]-1)
	{
#ifdef UGDEBUGOUT
	    pout()<<"\t connected Pair in z domain boundary \t"<<a_vof<<"\t bndryDir= "<<4<<endl;
#endif
	    UGData face_2d;
	    const bool correctFace=setUGData(face_2d,
					     a_pairPoints[ipair],
					     a_vof,
					     VolIndex(),
					     a_ebisBox,
					     a_vectDx,
					     4);
	    if (correctFace && !isUGDataInStencil(face_2d, a_vofStencil2D))
	    {
		a_vofStencil2D.push_back(face_2d);
		a_savedIndex.push_back(ipair);
	    }
	    isBoundary = true;
	}
	
	if (!isBoundary)
	{
	    Vector<VolIndex> vofList;
	    const int size2 = a_pairPoints[ipair].alternIVs.size();
	    for (int ii=0;ii<size2;ii++)
	    {
		if (a_domain.contains(a_pairPoints[ipair].alternIVs[ii]))
		{
		    
		    vofList.append(a_ebisBox.getVoFs(a_pairPoints[ipair].alternIVs[ii]));
		}
	    }
	    Vector<VolIndex> vofList2;
	    EBArith::getAllVoFsWithinRadius(vofList2, a_vof, a_ebisBox, 1);
	    vofList.append(vofList2);
	    
	    int sizeVoFs = vofList.size();
	    for (int ivof=0;ivof<sizeVoFs;ivof++)
	    {
		Vector<RealVect> checkPoints;
		Vector<IntVect>  checkSides;
		Vector<PairPoints> checkPairs;
		
		const VolIndex& pairVoF = vofList[ivof];
#ifdef UGDEBUGOUT
		pout()<< "\t checking pairVoF:  "<<pairVoF<<"\t bndryArea= "<<a_ebisBox.bndryArea(pairVoF)<<endl;
#endif
		if (a_ebisBox.bndryArea(pairVoF)<Tol) continue;
		
		findIntersectionPoints(a_ifdata, checkPoints,checkSides,pairVoF,a_lev, m_refRatios);
		findPairPoints(checkPairs,checkPoints,checkSides,pairVoF.gridIndex());
#ifdef UGDEBUGOUT
		int size=checkPoints.size();
		for (int ipoint=0;ipoint<size;ipoint++)
		{
		    pout()<< "\t Edges2:  "<<checkPoints[ipoint]<<"\t"<<checkSides[ipoint]<<endl;
		}
		const int sizePair = checkPairs.size();
		pout()<<"\t -2 numPoints="<<size<<"\t numPairs="<<sizePair<<endl;
		for (int ipoint=0;ipoint<sizePair;ipoint++)
		{
		    pout()<< "\t Pair-2:  first= "<<checkPairs[ipoint].first<<"\t second= "<<checkPairs[ipoint].second;
		    const int size2 = checkPairs[ipoint].alternIVs.size();
		    for (int ii=0;ii<size2;ii++)
		    {
			pout()<< "\t \t alternIV-2="<<checkPairs[ipoint].alternIVs[ii];
		    }
		    pout()<<endl;
		}
#endif
		const bool found = checkVoFPairs(checkPairs, a_pairPoints[ipair], vofIV);   
		if (found)
		{
		    UGData face_2d;
#ifdef UGDEBUGOUT
		    pout()<<"connected Pair \t"<<a_vof<<"\t"<<vofList[ivof]<<"\t bndryDir= "<<-1<<endl;
#endif
		    const bool correctFace=setUGData(face_2d,
						     a_pairPoints[ipair],
						     a_vof,
						     vofList[ivof],
						     a_ebisBox,
						     a_vectDx,
						     -1);
		    if (correctFace && !isUGDataInStencil(face_2d, a_vofStencil2D))
		    {
			a_vofStencil2D.push_back(face_2d);
			a_savedIndex.push_back(ipair);
			break;
		    }
		}
	    }
	}
    }
}


void  
EBPetscCompGridRichards:: setUnderresolvedUGStencil(Vector<UGData>&           a_vofStencil2D,
						    Vector<PairPoints>&       a_pairPoints,
						    BaseIVFAB<Real>&          a_bndryArea,
						    IFData2<3>&               a_ifdata,
						    const Vector<int>&        a_savedIndex,
						    const VolIndex&           a_vof,
						    const int                 a_lev,
						    const EBISBox&            a_ebisBox,
						    const ProblemDomain&      a_domain,
						    const ProblemDomain&      a_domainZero,
						    const RealVect&           a_vectDx)
    
{
    CH_TIME("EBPetscCompGridRichards::setUnderresolvedUGStencil");
    Vector<PairPoints>      pairPoints;
    Vector<UGData>          newvofStencil2D;
    // const IntVect domainSize = a_domain.size();
    std::stack<RealVect> allStacks; // declare here to use in outpu
    const IntVect& vofIV = a_vof.gridIndex();
    Vector<int> stenIndex;
    setUnderresolvedPairs(a_pairPoints,
			  allStacks,
			  stenIndex,
			  a_ifdata,
			  a_vofStencil2D,
			  a_savedIndex,
			  a_vof,
			  a_lev);
    
    const int numPairs = a_pairPoints.size();

    // if (vofIV==IntVect(40,43,70))
    // {
    // 	const int sizePairs = a_pairPoints.size();
    // 	pout()<<"Underresolved case numPairs="<<sizePairs<<endl;
    // 	for (int ipoint=0;ipoint<sizePairs;ipoint++)
    // 	{
    // 	    pout()<< "Pair:  first= "<<a_pairPoints[ipoint].first<<"\t second= "<<a_pairPoints[ipoint].second;
    // 	    const int size2 = a_pairPoints[ipoint].alternIVs.size();
    // 	    for (int ii=0;ii<size2;ii++)
    // 	    {
    // 		pout()<< "\t alternIV="<<a_pairPoints[ipoint].alternIVs[ii];
    // 	    }
    // 	    pout()<<endl;
    // 	}
    // 	pout()<<endl;
    // }

    //Underresolved case
    // a_dumpFile<<"=> vof="<<a_vof
    // 	      <<"numPoints="<<numPairs
    // 	      <<"\t numPairs="<<numPairs//+allStacks.size()-1
    // 	      <<"\t size="<<0//newvofStencil2D.size() 
    // 	// <<"\t bndryAreaEB="<<ebisBox.bndryArea(vof)
    // 	// <<"\t normalEB="<<ebisBox.normal(vof)
    // 	      <<endl;
    for (int is=0;is<numPairs;is++)
    {
	const PairPoints& pair = a_pairPoints[is];
	const RealVect& first  = pair.first;
	const RealVect& second = pair.second;

	int included=-1;
	const VolIndex& conVoF = a_vofStencil2D[stenIndex[is]].vof;
	if (first[2]!=second[2])
	{
	    UGData& face2d = a_vofStencil2D[stenIndex[is]];
	    newvofStencil2D.push_back(face2d);

	    pairPoints.push_back(pair);
	    included=1;
	}
	else
	{
	    Vector<RealVect> conPoints;
	    Vector<IntVect>  conSides;
	    Vector<PairPoints>      conPairPoints;
	    std::stack<RealVect>    conAllStacks; // declare here to use in outpu
	    Vector<int> conStenIndex;
	    
	    const int numConPoints = findIntersectionPoints(a_ifdata, conPoints, conSides, conVoF, a_lev, m_refRatios);
	    const int numConPairs = findPairPoints(conPairPoints,conPoints,conSides,conVoF.gridIndex());
	    Vector<UGData> conVoFStencil2D;
	    Vector<int>        conSavedIndex;
	    setUGStencil(conVoFStencil2D,  conSavedIndex, a_ifdata, conPairPoints, conVoF, a_lev, a_ebisBox, a_domain, a_domainZero, a_vectDx);

	    if (numConPairs!=numConPoints)
	    {
		setUnderresolvedPairs(conPairPoints,
				      conAllStacks,
				      conStenIndex,
				      a_ifdata,
				      conVoFStencil2D,
				      conSavedIndex,
				      conVoF,
				      a_lev);
	    }
	    const bool found = checkVoFPairs(conPairPoints, pair, vofIV);   
	    if (found)
	    {
		UGData& face2d = a_vofStencil2D[stenIndex[is]];
		newvofStencil2D.push_back(face2d);
		
		pairPoints.push_back(pair);
		included=1;
	    }
	}

	// RealVect grad1 = RealVect::Zero;
	// RealVect grad2 = RealVect::Zero;
	// computeGradImplicitFunction(grad1, first);
	// computeGradImplicitFunction(grad2, second);
	// a_dumpFile<<"\t first="<<pair.first
	// 	  <<"\t second="<<pair.second
	// 	  <<"\t conVoF="<< conVoF
	// 	  <<"\t normal="<<a_vofStencil2D[stenIndex[is]].faceData
	// 	  <<"\t bndryDir="<<a_vofStencil2D[stenIndex[is]].bndryDir
	// 	  <<"\t gradFirst=" <<grad1
	// 	  <<"\t gradSecond="<<grad2
	// 	  <<"\t dotProduct="<<included
	// 	  <<endl;
    }

    a_vofStencil2D.clear();
    a_vofStencil2D = newvofStencil2D;



    if(a_vofStencil2D.size()!=pairPoints.size()) MayDay::Abort("numPairs is not equal vofStencil size!");


    // if (a_vof.gridIndex()==IntVect(40,43,70))
    // {
	 
    // 	a_dumpFile<<"=> vof="<<a_vof
    // 		  <<"numPoints="<<pairPoints.size()
    // 		  <<"\t numPairs="<<pairPoints.size()+allStacks.size()-1
    // 		  <<"\t size="<<a_vofStencil2D.size() 
    // 		  <<endl;
    // 	for (int is=0;is<a_vofStencil2D.size();is++)
    // 	{
    // 	    const PairPoints& pair = pairPoints[is];
    // 	    const RealVect& first  = pair.first;
    // 	    const RealVect& second = pair.second;
	    
    // 	    RealVect grad1 = RealVect::Zero;
    // 	    RealVect grad2 = RealVect::Zero;
    // 	    computeGradImplicitFunction(grad1, first);
    // 	    computeGradImplicitFunction(grad2, second);
	    
    // 	    a_dumpFile<<"\t first="<<pair.first
    // 		      <<"\t second="<<pair.second
    // 		      <<"\t conVoF="<< a_vofStencil2D[is].vof
    // 		      <<"\t normal="<< a_vofStencil2D[is].faceData
    // 		      <<"\t bndryDir="<<a_vofStencil2D[is].bndryDir
    // 		      <<"\t gradFirst=" <<grad1
    // 		      <<"\t gradSecond="<<grad2
    // 		      <<"\t dotProduct="<<1
    // 		      <<endl;
    // 	}
    // 	const int allstacksize=allStacks.size();
    // 	for (int is=0;is<allstacksize-1;is++)
    // 	{
    // 	    const RealVect first = allStacks.top();
    // 	    allStacks.pop();
    // 	    const RealVect second = allStacks.top();
	    
    // 	    a_dumpFile<<"\t first="<<first
    // 		      <<"\t second="<<second
    // 		      <<"\t conVoF="<<VolIndex()
    // 		      <<"\t normal="<<RealVect::Zero//newvofStencil2D[0].faceData[0]
    // 		      <<"\t bndryDir="<<-1//newvofStencil2D[0].bndryDirs[0]
    // 		      <<"\t gradFirst="<<RealVect::Zero
    // 		      <<"\t gradSecond="<<RealVect::Zero
    // 		      <<"\t dotProduct="<<0
    // 		      <<endl;
    // 	}
    // 	a_dumpFile<<endl;
    // }
    Vector<RealVect> normals;
    Vector<Real>     areas;
    Vector<int>      indexCluster;

    // pout()<<a_vof
    // 	  <<"\t stencilSize="<<a_vofStencil2D.size()
    // 	  <<"\t bndryAreaEB="<<a_ebisBox.bndryArea(a_vof)
    // 	  <<"\t normalEB="<<a_ebisBox.normal(a_vof);

    const bool goodClust=computeClustersData(normals,
					     areas,
					     indexCluster, 
					     pairPoints,
					     a_vectDx);
    if (!goodClust)
    {
	cerr<<"WARNING: Nonmatching number of points and pairs in the cluster! "<<a_vof<<"\t numPairs="<<pairPoints.size()<<endl;
	cerr<<"Try to use more resolved grid! "<<endl;

	// a_dumpFile<<"=> vof="<<a_vof
	// 	  <<"numPoints="<<pairPoints.size()
	// 	  <<"\t numPairs="<<pairPoints.size()+allStacks.size()-1
	// 	  <<"\t size="<<a_vofStencil2D.size() 
	// 	  <<endl;
	// for (int is=0;is<a_vofStencil2D.size();is++)
	// {
	//     const PairPoints& pair = pairPoints[is];
	//     const RealVect& first  = pair.first;
	//     const RealVect& second = pair.second;
	    
	//     RealVect grad1 = RealVect::Zero;
	//     RealVect grad2 = RealVect::Zero;
	//     computeGradImplicitFunction(grad1, first);
	//     computeGradImplicitFunction(grad2, second);
	    
	//     a_dumpFile<<"\t first="<<pair.first
	// 	      <<"\t second="<<pair.second
	// 	      <<"\t conVoF="<< a_vofStencil2D[is].vof
	// 	      <<"\t normal="<< a_vofStencil2D[is].faceData
	// 	      <<"\t bndryDir="<<a_vofStencil2D[is].bndryDir
	// 	      <<"\t gradFirst=" <<grad1
	// 	      <<"\t gradSecond="<<grad2
	// 	      <<"\t dotProduct="<<1
	// 	      <<endl;
	// }
	// const int allstacksize=allStacks.size();
	// for (int is=0;is<allstacksize-1;is++)
	// {
	//     const RealVect first = allStacks.top();
	//     allStacks.pop();
	//     const RealVect second = allStacks.top();
	    
	//     a_dumpFile<<"\t first="<<first
	// 	      <<"\t second="<<second
	// 	      <<"\t conVoF="<<VolIndex()
	// 	      <<"\t normal="<<RealVect::Zero//newvofStencil2D[0].faceData[0]
	// 	      <<"\t bndryDir="<<-1//newvofStencil2D[0].bndryDirs[0]
	// 	      <<"\t gradFirst="<<RealVect::Zero
	// 	      <<"\t gradSecond="<<RealVect::Zero
	// 	      <<"\t dotProduct="<<0
	// 	      <<endl;
	// }
	// a_dumpFile<<endl;
	// exit(1);
    }
    
    Real cell_area=0.;
    Vector<Vector<Real> > slopesCell;

    //RealVect normal = a_ebisBox.normal(a_vof);
    const int clustSize = normals.size();
    for (int iclust=0;iclust<clustSize;iclust++)
    {
	cell_area +=areas[iclust];

	//leave default slopes for a sinle subcell
	//else take computed least squared
	if (clustSize<=1) continue;
	const RealVect& normal = normals[iclust];
	Vector<Real> slopes;
	for (int idir=0; idir<2; idir++)
	{
	    if (fabs(normal[2])>1.e-12)
	    {
		slopes.push_back(-normal[idir]*a_vectDx[2]/(normal[2]*a_vectDx[idir]));
	    }
	    else
	    {
		slopes.push_back(0.0);
		slopes.push_back(0.0);
	    }
	}
	slopesCell.push_back(slopes);
    }
    
    //store data for underresolved cell
    a_bndryArea(a_vof,0) = cell_area;
    if (clustSize>1) 
    {
	for (int is=0;is<a_vofStencil2D.size();is++)
	{
	    a_vofStencil2D[is].slopes = slopesCell[indexCluster[is]];
	}    
    }
}

bool  
EBPetscCompGridRichards:: computeClustersData(Vector<RealVect>&         a_normals,
					      Vector<Real>&             a_areas,
					      Vector<int>&              a_indexCluster, 
					      const Vector<PairPoints>& a_pairs,
					      const RealVect&           a_vectDx)
{
    CH_TIME("EBPetscCompGridRichards::computeClusterData");
    Vector<Vector<RealVect> > clusters;
    const int numPairs = a_pairs.size();
    for (int is=0;is<numPairs;is++)
    {
	const PairPoints& pair = a_pairs[is];
	const RealVect& first  = pair.first;
	const RealVect& second = pair.second;
	bool in_cluster = false;
	const int numClust = clusters.size();
	for (int iclust=0;iclust<numClust; iclust++)
	{
	    Vector<RealVect>& clust = clusters[iclust];
	    const int clustsize = clusters[iclust].size();
	    bool addfirst  =true;
	    bool addsecond =true;
	    for(int ipoint=0;ipoint<clustsize;ipoint++)
	    {
		if(first == clust[ipoint])
		{
		    in_cluster = true;
		    addfirst   = false;
		}
		    
		if(second==clust[ipoint])
		{
		    in_cluster = true;
		    addsecond = false;
		}
		if (!addfirst && !addsecond) break;
	    }
	    if (in_cluster)
	    {
		a_indexCluster.push_back(iclust);
		if (addfirst)  clust.push_back(first);
		if (addsecond) clust.push_back(second);
		break;
	    }
		

	    std::stack<RealVect> Stack;
	    Vector<RealVect>     visited;
	    const RealVect& clustpoint = clust[0];

	    Stack.push(first);
	    visited.push_back(first);
	    in_cluster = isPointsConnected(Stack, visited, clustpoint, a_pairs, RealVect(-100.,-100., -100.));
	    // pout()<<"\t \t first="<<first<<"\t second"<<clustpoint<<"\t isConn="<<in_cluster<<endl;
		
	    if (in_cluster) 
	    {
		clust.push_back(first);
		clust.push_back(second);
		const int stacksize=Stack.size();
		for (int s=0;s<stacksize;s++)
		{
		    const RealVect point = Stack.top();
		    Stack.pop();
		    bool addPoint=true;
		    const int clustsize1 = clust.size();
		    for(int ip=0;ip<clustsize1;ip++)
		    {
			if(point==clust[ip])
			{
			    addPoint = false;
			    break;
			}
		    }
		    if (addPoint) clust.push_back(point);
		}
		a_indexCluster.push_back(iclust);
		in_cluster = true;
		break;
	    }
	}
	if (!in_cluster)
	{
	    Vector<RealVect> clust;
	    clust.push_back(first);
	    clust.push_back(second);
	    a_indexCluster.push_back(clusters.size());
	    clusters.push_back(clust);
	}
    }

    const int numClust = clusters.size();
    int totalSize=0;
    for (int iclust=0;iclust<numClust;iclust++)
    {
	const Vector<RealVect>& clust = clusters[iclust];

	const Real area =  computeShoeLaceArea(clust, m_vectDxs[0]);//a_vectDx);
	a_areas.push_back(area);

	RealVect normal;
	const int clustsize = clusters[iclust].size();
	totalSize +=clustsize;
	
	if(clustsize==3)
	{
	    const RealVect d1 = clust[1]-clust[0];
	    const RealVect d2 = clust[2]-clust[0];
	    normal[0] = d1[1]*d2[2]-d1[2]*d2[1];
	    normal[1] = d1[2]*d2[0]-d1[0]*d2[2];
	    normal[2] = d1[0]*d2[1]-d1[1]*d2[0];
	    normal /=normal.vectorLength();
	}
	else if (clustsize>3)
	{
	    RealVect centroid=RealVect::Zero;
	    for(int ipoint=0;ipoint<clustsize;ipoint++)
	    {
		centroid +=clust[ipoint];
	    }
	    centroid /=clustsize;
	    Real xy=0.;
	    Real xz=0.;
	    Real yz=0.;
	    Real xx=0.;
	    Real yy=0.;
	    Real zz=0.;
	    for(int ipoint=0;ipoint<clustsize;ipoint++)
	    {
		const RealVect d = clust[ipoint]-centroid;
		xy += d[0]*d[1];
		xz += d[0]*d[2];
		yz += d[1]*d[2];
		xx += d[0]*d[0];
		yy += d[1]*d[1];
		zz += d[2]*d[2];
	    }
	    normal[0] = yz*xy-xz*yy;
	    normal[1] = xz*xy-yz*xx;
	    normal[2] = xx*yy-xy*xy;
	    normal /=normal.vectorLength();
	}
	else 
	{
	    const int clustsize = clusters[iclust].size();
	    pout()<<"\n \t clusterSize="<<clustsize<<"\t  cluster="<<clusters[iclust]<<endl;
	    MayDay::Abort("Wrong number of points in the unstructured grid subcell!");
	}
	a_normals.push_back(normal);
    }
    if (totalSize!=numPairs) return false;

    // CH_assert(totalSize==numPairs);
    
//    pout()<<"\t num clusters="<<clusters.size()<<"  sum_sizes="<<totalSize<<endl;
    // for (int iclust=0;iclust<clusters.size();iclust++)
    // {
    // 	const int clustsize = clusters[iclust].size();
    // 	pout()<<"\t clusterSize="<<clustsize<<"\t normal="<<a_normals[iclust]<<"\t area="<<a_areas[iclust]<<endl;
	
    // 	pout()<<"\t "<<clusters[iclust]<<endl;
    // }
    // pout()<<endl;
    return true;
}


void  
EBPetscCompGridRichards:: resolveCell(Vector<PairPoints>&       a_resolvedPairs,
				      IFData2<3>&               a_ifdata,
				      const IntVect&            a_vofIV,
				      const int                 a_lev,
				      const Vector<int>&        a_ratios)
{
    CH_TIME("EBPetscCompGridRichards::resolveCell");
    int lev = a_lev;
    Vector<int> ratios = a_ratios;
    Box refCell(a_vofIV, a_vofIV);
    refCell.refine(2);
    ratios.push_back(2);
    lev++;
    // pout()<<endl;
    // pout()<<"\t  ilev="<<lev<<" iv="<<a_vofIV;
    const RealVect cellCorner = RealVect(a_vofIV);
	
    for (BoxIterator boxit(refCell); boxit.ok(); ++boxit)
    {
	Vector<RealVect> resolvedPoints;
	Vector<IntVect>  resolvedSides;
	Vector<PairPoints> resolvedPairs;
	
	const IntVect& ivFine = boxit();
	const VolIndex vofFine=VolIndex(ivFine,0);
	
	const int numPoints=findIntersectionPoints(a_ifdata, resolvedPoints, resolvedSides,vofFine, lev, ratios);
	if (!numPoints) continue;
	const int numPairs=findPairPoints(resolvedPairs,resolvedPoints,resolvedSides,ivFine);
	if (numPairs!=numPoints)
	{
	    resolvedPairs.clear();
	    resolveCell(resolvedPairs, a_ifdata, ivFine, lev, ratios);
	}
	//assumes that there is only one resolved point on each edge


	a_resolvedPairs.append(resolvedPairs);
//	pout()<<" numPairs="<<numPairs;
    }
}


void  
EBPetscCompGridRichards:: setUnderresolvedPairs(Vector<PairPoints>&       a_pairPoints,
						std::stack<RealVect>&     a_allStacks,
						Vector<int>&              a_stenIndex,
						IFData2<3>&               a_ifdata,
						const Vector<UGData>&     a_vofStencil2D,
						const Vector<int>&        a_savedIndex,
						const VolIndex&           a_vof,
						const int                 a_lev)
    
{
    CH_TIME("EBPetscCompGridRichards::setUnderresolvedPairs");
    const Real eps=1.e-11 / pow(2.,a_lev);
    // const IntVect domainSize = a_domain.size();
    Vector<PairPoints> pairPoints = a_pairPoints;
    const IntVect& vofIV = a_vof.gridIndex();
    const RealVect cellCorner = RealVect(vofIV);
    //find store all resolved pairs 
    Vector<PairPoints> returnPairs;
    Vector<PairPoints> allResolvedPairs;
    //pout()<<"check_1  iv="<<vofIV<<endl;

    {
	Vector<int> ratios;
	for (int ilev = 0;ilev <a_lev; ++ilev)
	{
	    ratios.push_back(m_refRatios[ilev]);
	}
	
	//recursive call to resolve a cell	
	resolveCell(allResolvedPairs, a_ifdata, vofIV, a_lev, ratios);

	//pout()<<"Underresolved IV="<<vofIV<<"\t all_resolved pairs="<<allResolvedPairs.size()<<endl;

	//multiple roots on edge can result in different edge intersection on coarse and fine grid
       	//need to update all edge intersections to avoid multiple roots problem
	const int sizep=allResolvedPairs.size();
	for (int ipair=0;ipair<sizep;ipair++)
	{
	    for (int ipoint=0;ipoint<2;ipoint++)
	    {
		RealVect& point  = (ipoint==0 ? allResolvedPairs[ipair].first : allResolvedPairs[ipair].second);

		//check if point is an edge point
		IntVect edgeInds;
		if(isPointOnCoarseEdge(edgeInds, point, cellCorner, eps)!=2) continue;
		
		int varDir=-1;
		for(int idir=0;idir<3;idir++)
		{
		    if (edgeInds[idir]==-1)
		    {
			varDir=idir;
			break;
		    }
		}
		
		//find point on edge from underresolved grid
		for (int is=0;is<a_vofStencil2D.size();is++)
		{
		    const int ind=a_savedIndex[is];
		    PairPoints& pair = pairPoints[ind];
		    RealVect& first  = pair.first;
		    RealVect& second = pair.second;
		    IntVect edgeFirst;
		    if(isPointOnCoarseEdge(edgeFirst, first, cellCorner, eps)==2)
		    {
			bool needsCorrect=true;
			for(int idir=0;idir<3;idir++)
			{
			    if (idir==varDir) continue;
			    if (edgeFirst[idir]!=edgeInds[idir]) 
			    {
				needsCorrect=false;
				break;
			    }
			}
			if (needsCorrect)
			{
			    first[varDir] = point[varDir];
			}
		    }
		    
		    IntVect edgeSecond;
		    if(isPointOnCoarseEdge(edgeSecond, second, cellCorner, eps)==2)
		    {
			bool needsCorrect=true;
			for(int idir=0;idir<3;idir++)
			{
			    if (idir==varDir) continue;
			    if (edgeSecond[idir]!=edgeInds[idir]) 
			    {
				needsCorrect=false;
				break;
			    }
			}
			if (needsCorrect)
			{
			    second[varDir] = point[varDir];
			}
		    }
		}//vofstencil_iterator
	    }//ipoint
	}//allresoved_iterator
    }//local scope

    // if (a_vof.gridIndex()==IntVect(40,43,70))
    // {
	 
    // 	const int sizep=allResolvedPairs.size();
    // 	a_dumpFile<<"=> vof="<<a_vof
    // 		  <<"numPoints="<<sizep
    // 		  <<"\t numPairs="<<sizep
    // 		  <<"\t size="<<sizep 
    // 		  <<endl;
    // 	for (int is=0;is<sizep;is++)
    // 	{
    // 	    const PairPoints& pair = allResolvedPairs[is];
    // 	    const RealVect& first  = pair.first;
    // 	    const RealVect& second = pair.second;
	    
    // 	    RealVect grad1 = RealVect::Zero;
    // 	    RealVect grad2 = RealVect::Zero;
	    
    // 	    a_dumpFile<<"\t first="<<pair.first
    // 		      <<"\t second="<<pair.second
    // 		      <<"\t conVoF="<< VolIndex()
    // 		      <<"\t normal="<< grad1
    // 		      <<"\t bndryDir="<<-1
    // 		      <<"\t gradFirst=" <<grad1
    // 		      <<"\t gradSecond="<<grad2
    // 		      <<"\t dotProduct="<<1
    // 		      <<endl;
    // 	}
    // 	a_dumpFile<<endl;
    // }
    
    //check connections
    for (int is=0;is<a_vofStencil2D.size();is++)
    {
	const int ind=a_savedIndex[is];
	const PairPoints& pair = pairPoints[ind];
	const RealVect& first  = pair.first;
	const RealVect& second = pair.second;
	RealVect tempv = RealVect::Zero;
	//first exclude pair on the same top or bottom edge 
	if (fabs(first[2]-second[2])<eps)
	{
	    if (fabs(first[0]-second[0])<eps)
	    {
		if((fabs(first[1]-cellCorner[1])>eps     || fabs(first[1]-cellCorner[1]-1)>eps)
		   && (fabs(second[1]-cellCorner[1])>eps || fabs(second[1]-cellCorner[1]-1)>eps)) continue;
	    }
	    if (fabs(first[1]-second[1])<eps)
	    {
		if((fabs(first[0]-cellCorner[0])>eps     || fabs(first[0]-cellCorner[0]-1)>eps)
		   &&(fabs(second[0]-cellCorner[0])>eps  || fabs(second[0]-cellCorner[0]-1)>eps)) continue;
	    }
	}
	
	//make internal face are always included
	if (fabs(first[2]-second[2])>eps)
	{
	    returnPairs.push_back(a_pairPoints[ind]);
	    a_stenIndex.push_back(is);
	}
	else
	{
	    std::stack<RealVect> Stack;
	    Vector<RealVect> visited;
	    Stack.push(first);
	    visited.push_back(first);

		
	    const bool is_connected = isPointsConnected(Stack,
							visited,
							second,
							allResolvedPairs,
							cellCorner);
	    // if (a_vof.gridIndex()==IntVect(40,43,70))
	    // {
	    // 	if (fabs(first[2]-70)<eps && fabs(second[2]-70)<eps)
	    // 	{
	    // 	    pout()<< "\t Check Pair:  first= "<<first<<"\t second= "<<second<<"\t isConnected="<<is_connected<<endl;
		
	    // 	    // const int sizePairs = allResolvedPairs.size();
	    // 	    // pout()<<"\t \t Check case size="<<sizePairs<<endl;
	    // 	    // for (int ipoint=0;ipoint<sizePairs;ipoint++)
	    // 	    // {
	    // 	    // 	pout()<< "Pair:  first= "<<allResolvedPairs[ipoint].first<<"\t second= "<<allResolvedPairs[ipoint].second;
	    // 	    // 	pout()<<endl;
	    // 	    // 	a_allStacks.push(allResolvedPairs[ipoint].first);
	    // 	    // 	a_allStacks.push(allResolvedPairs[ipoint].second);
	    // 	    // }
	    // 	    // pout()<<endl;
	    // 	}
	    // }

	    if (is_connected)
	    {
		returnPairs.push_back(a_pairPoints[ind]);
		a_stenIndex.push_back(is);
		const int stacksize=Stack.size();
		for (int is=0;is<stacksize;is++)
		{
		    a_allStacks.push(Stack.top());
		    Stack.pop();
		}
	    }
	}
    }//vofStencilIterator
    a_pairPoints.clear();
    a_pairPoints=returnPairs;
}


int
EBPetscCompGridRichards:: isPointOnCoarseEdge(IntVect& a_edgeInds,
					      const RealVect& a_point,
					      const RealVect& a_cellCorner,
					      const Real a_eps)
{
    CH_TIME("EBPetscCompGridRichards::isPointOnCoarseEdge");
    a_edgeInds = IntVect(-1,-1,-1);
    int count=0;
    for(int idir=0;idir<3;idir++)
    {
	if(fabs(a_point[idir]-a_cellCorner[idir])<a_eps)
	{
	    a_edgeInds[idir]=0;
	    count++;
	}
	if(fabs(a_point[idir]-a_cellCorner[idir]-1)<a_eps)
	{
	    a_edgeInds[idir]=1;
	    count++;
	}
    }
    return count;
}

bool  
EBPetscCompGridRichards:: isUGDataInStencil(const UGData&         a_ugdata,
					    const Vector<UGData>& a_vofStencil)
{
    CH_TIME("EBPetscCompGridRichards::isUGDatainStencil");
    const int size = a_vofStencil.size();
    for (int ivof=0;ivof<size;ivof++)
    {
	if (a_ugdata == a_vofStencil[ivof]) return true;
    }

    return false;
}

bool  
EBPetscCompGridRichards:: isPointsConnected(std::stack<RealVect>& a_stack,
					    Vector<RealVect>& a_visited,
					    const RealVect& a_secondPoint,
					    const Vector<PairPoints>& a_resolvedPairs,
					    const RealVect& a_cellCorner)
{
    CH_TIME("EBPetscCompGridRichards::isPointsConnected");
    if (a_stack.empty()) return false;

    const Real eps=1e-11;
    const RealVect& front = a_stack.top();
    
    //use dfs search (better to use trees for that)
    const int size=a_resolvedPairs.size();
    for (int ipair=0;ipair<size;ipair++)
    {
	const RealVect& first  = a_resolvedPairs[ipair].first;
	const RealVect& second = a_resolvedPairs[ipair].second;

	if ((second-front).vectorLength()<eps || (first-front).vectorLength()<eps)
	{
	    if ((second-a_secondPoint).vectorLength()<eps || (first-a_secondPoint).vectorLength()<eps)
	    {
		a_stack.push(a_secondPoint);
		return true;
	    }
	    bool includeFirst=true;
	    bool includeSecond=true;
	    const int sizeVis = a_visited.size();
	    for(int ivis=0;ivis<sizeVis;ivis++)
	    {
		// if ((first-a_visited[ivis]).vectorLength()<eps) includeFirst=false;
		// if ((second-a_visited[ivis]).vectorLength()<eps) includeSecond=false;
		if (first==a_visited[ivis]) includeFirst=false;
		if (second==a_visited[ivis]) includeSecond=false;
	    }
	    if (includeFirst)
	    {
		// we consider only internal points
		if((fabs(first[0]-a_cellCorner[0])>eps) && (fabs(first[1]-a_cellCorner[1])>eps) && (fabs(first[0]-a_cellCorner[0]-1)>eps) && (fabs(first[1]-a_cellCorner[1]-1)>eps))
		{
		    a_stack.push(first);
		    a_visited.push_back(first);
		    if ( isPointsConnected(a_stack,
					   a_visited,
					   a_secondPoint,
					   a_resolvedPairs,
					   a_cellCorner) ) return true;
		}
	    }
	    if (includeSecond)
	    {
		if((fabs(second[0]-a_cellCorner[0])>eps) && (fabs(second[1]-a_cellCorner[1])>eps) && (fabs(second[0]-a_cellCorner[0]-1)>eps) && (fabs(second[1]-a_cellCorner[1]-1)>eps))
		{
		    a_stack.push(second);
		    a_visited.push_back(second);
		    if (isPointsConnected(a_stack,
					  a_visited,
					  a_secondPoint,
					  a_resolvedPairs,
					  a_cellCorner) ) return true;
		}
	    }
	}
	else
	{
	    continue; 
	}
    }
	
    a_stack.pop();
    return false; 
}

Real  
EBPetscCompGridRichards:: computeShoeLaceArea(const Vector<RealVect>&   a_points,
					      const RealVect&           a_vectDx)
{
    CH_TIME("EBPetscCompGridRichards::computeSoeLaceArea");
    const int size = a_points.size();
    if (!size) return 0.;

    //first we need to sort them to get polygon order from center
    //perhaps there is a better algorithm to do that, but there are usually few points
    //I used here polar coord convertion
    Vector<RealVect>  points = a_points;
    Vector<Real> thetas;
    thetas.resize(size,0.);
    RealVect center = RealVect::Zero;
    for (int i=0;i<size;i++)
    {
	center +=points[i];
//	pout()<<points[i]<<"\t";
    }

    center /=size;
    
    for (int i=0;i<size;i++)
    {
	points[i][0] -= center[0];
	points[i][1] -= center[1];
	
	if (points[i][0]==0 && points[i][1]>0)
	{
	    thetas[i]=0.5*M_PI;
	}
	else if (points[i][0]==0 && points[i][1]<0)
	{
	    thetas[i]=-0.5*M_PI;
	}
	else
	{
	    thetas[i]=atan(points[i][1]/points[i][0]);
	    if (points[i][0]<0 && points[i][1]>=0)
	    {
		thetas[i] +=M_PI;
	    }
	    else if (points[i][0]<0 && points[i][1]<0)
	    {
		thetas[i] -=M_PI;
	    }
	}
	while(thetas[i]<0.) thetas[i] +=2.*M_PI;
	CH_assert(thetas[i]<=2.*M_PI);
    }
    for (int i=0;i<size;i++)
    {
	Real mintheta=2.*M_PI;
	int minIndex=0;
	for (int j=i;j<size;j++)
	{
	    if (thetas[j]<mintheta)
	    {
		mintheta=thetas[j];
		minIndex=j;
	    }
	}
	if (minIndex!=i)
	{
	    RealVect tempVect = points[i];
	    Real temptheta    = thetas[i];
	    
	    points[i] = points[minIndex];
	    thetas[i] = thetas[minIndex];
	    
	    points[minIndex]=tempVect;
	    thetas[minIndex]=temptheta;
	}
    }

    Real sum=0.;
    const Real dxdy = a_vectDx[0]*a_vectDx[1];
    
    for (int i=0;i<size-1;i++)
    {
	sum += points[i][0]*points[i+1][1] - points[i][1]*points[i+1][0];
    }
    sum += points[size-1][0]*points[0][1] - points[0][0]*points[size-1][1];
    
    sum *=0.5*dxdy;

//    pout()<<sum;
//    pout()<<endl;

    return fabs(sum);
}

    
bool  
EBPetscCompGridRichards:: setUGData(UGData&           a_ugdata,
				    const PairPoints& a_pair,
				    const VolIndex&   a_vofCenter,
				    const VolIndex&   a_vofConnect,
				    const EBISBox&    a_ebisBox,
				    const RealVect&   a_vectDx,
				    const int         a_bndryDir)
{
    CH_TIME("EBPetscCompGridRichards::setupUGData");

    RealVect faceData=RealVect::Zero;
    VolIndex vof=VolIndex();
	
    RealVect vofCenterLoc  = EBArith::getVoFLocation(a_vofCenter, a_vectDx, RealVect::Zero);
    vofCenterLoc          += a_ebisBox.bndryCentroid(a_vofCenter)*a_vectDx;
    
    RealVect point0 = a_pair.first;
    RealVect point1 = a_pair.second;
    point0 *=m_vectDxs[0];//a_vectDx;
    point1 *=m_vectDxs[0];//a_vectDx;

    //check if face belongs to this vof
    // const RealVect reld0 = point0-vofCenterLoc;
    // const RealVect reld1 = point1-vofCenterLoc;
    // if ((reld0[0]*reld1[0]>0.)&&(reld0[1]*reld1[1]>0.)) return false;

    const RealVect delta = point0-point1;
    const RealVect centroid = 0.5*((point0+point1)/m_vectDxs[0]);//a_vectDx);//we need centroid in reference frame
    const Real length2d = sqrt(delta[0]*delta[0] + delta[1]*delta[1]);
		    
    //save face length
    faceData[2] = length2d;
    CH_assert(fabs(length2d) > 1e-16);

    if (a_bndryDir<0)
    {
	RealVect vofCentConnectVect = RealVect(a_vofConnect.gridIndex()-a_vofCenter.gridIndex());
	const Real norm0 = -delta[1]/length2d;
	const Real norm1 =  delta[0]/length2d;
	if ((!vofCentConnectVect[0])&&(!vofCentConnectVect[1]))
	{
	    const Real sign = sgn(vofCentConnectVect[2]);
	    computeGradImplicitFunction(vofCentConnectVect, centroid);
	    vofCentConnectVect *=-sign;
	}
	Real dotnorm = vofCentConnectVect[0]*norm0 + vofCentConnectVect[1]*norm1;

#ifdef UGDEBUGOUT
 	pout()<<"\t vofCentConnectVectc="<<vofCentConnectVect<<"\t dotnorm2D="<<dotnorm<<endl;
#endif	
	    
	// if(!dotnorm)
	// {
	//     MayDay::Error(" setUGData: bad initialization for unstructured grid - dotproduct is zero");
	// }

	faceData[0] = sgn(dotnorm)*norm0;
	faceData[1] = sgn(dotnorm)*norm1;

	vof = a_vofConnect;
	// const Real d0 = sqrt(vofCenterVect[0]*vofCenterVect[0]+vofCenterVect[1]*vofCenterVect[1]);
	// const Real d1 = sqrt(vofConnectVect[0]*vofConnectVect[0]+vofConnectVect[1]*vofConnectVect[1]);
	// a_ugdata.interpWeight = 1.0/d0/(1.0/d0+1.0/d1);
    }
    else
    {
	switch (a_bndryDir)
	{
	case 0:
	    faceData[0] = -1.0;
	    faceData[1] = 0.0;
	    break;
	case 1:
	    faceData[0] = 1.0;
	    faceData[1] = 0.0;
	    break;
	case 2:
	    faceData[0] = 0.0;
	    faceData[1] = -1.0;
	    break;
	case 3:
	    faceData[0] = 0.0;
	    faceData[1] = 1.0;
	    break;
	case 4:
	    faceData[0] = 0.0;
	    faceData[1] = 0.0;
	    break;
	default:
	    faceData[0] = 0.0;
	    faceData[1] = 0.0;
	    break;
	}	    
    }
    a_ugdata.vof = vof;
    a_ugdata.faceData=faceData;
    a_ugdata.bndryDir=a_bndryDir;
    //a_ugdata.numFaces++;
    
    return true;
}

bool
EBPetscCompGridRichards:: checkVoFPairs(const Vector<PairPoints>& a_checkPairs,
					const PairPoints&         a_vofPair,
					const IntVect&            a_vofIV)
{
    CH_TIME("EBPetscCompGridRichards::checkVoFPairs");
    const Real eps = 1.e-15;

    const int size=a_checkPairs.size();
    for (int i=0;i<size;i++)
    {
	if(((a_checkPairs[i].first - a_vofPair.first).vectorLength()<eps  && (a_checkPairs[i].second -a_vofPair.second).vectorLength()<eps) ||
	   ((a_checkPairs[i].first - a_vofPair.second).vectorLength()<eps && (a_checkPairs[i].second - a_vofPair.first).vectorLength()<eps))
	{
	    return true;
	    const int size2 = a_checkPairs[i].alternIVs.size();
	    for (int ii=0;ii<size2;ii++)
	    {
		//pout()<<a_vofIV<<"-"<<a_checkPairs[i].alternIVs[ii]<<endl;
		if(a_vofIV==a_checkPairs[i].alternIVs[ii])
		{
		    return true;
		}
	    }
	}
    }
    return false;
}

int  
EBPetscCompGridRichards:: findPairPoints(Vector<PairPoints>& a_pairPoints,
					 const Vector<RealVect>& a_points,
					 const Vector<IntVect>&  a_sides,
					 const IntVect& a_iv)
				{
    CH_TIME("EBPetscCompGridRichards::findPairPoints");
    int size=a_points.size();
    for (int ipoint=0;ipoint<size;ipoint++)
    {
	//const RealVect& p0 = a_points[ipoint];
	
	for (int ii=ipoint+1;ii<size;ii++)
	{
	    //const RealVect& p1 = a_points[ii];
	    //check only faces
	    //if (p0[0]==p1[0] || p0[1]==p1[1] || p0[2]==p1[2])
	    {
		const IntVect& side0 = a_sides[ipoint];
		const IntVect& side1 = a_sides[ii];
		
		IntVect pairIV=a_iv;
		for (int idir=0;idir<3;idir++)
		{
		    pairIV[idir] +=(side0[idir]+side1[idir])/2;
		}
		PairPoints pair;
		pair.first  = a_points[ipoint];
		pair.second = a_points[ii];

		//currently no support for multicells (more than one vofs in one EBcell)
		if (pairIV!=a_iv) pair.alternIVs.push_back(pairIV);

		//we need alternative ivs, since sometimes it shifts
		for (int idir=0;idir<3;idir++)
		{
		    const int shift=(side0[idir]+side1[idir])/2;
		    //  pout()<< shift<<" ";
		    if (shift)
		    {
			IntVect alternIV = pairIV;
			alternIV[idir] -=shift;
			if (alternIV!=a_iv) pair.alternIVs.push_back(alternIV);
			// pout()<<"\t\t check altern"<<alternIV;
			// for (int idir2=0;idir2<3;idir2++)
			// {
			//     pout()<<"\t idir="<<idir2<<" s0="<<side0[idir2]<<" s1="<<side1[idir2];
			//     if (idir2==idir) continue;
			//     const int shift2=(side0[idir2]+side1[idir2])/2;
			//     if (shift2)
			//     {
			// 	IntVect alternIV2 = alternIV;
			// 	alternIV2[idir2] -=shift2;
			// 	pair.alternIVs.push_back(alternIV2);
			//     }
			// }
			//	pout()<<alternIV<<" ";
		    }
		}
		//pout() <<endl;
		if (pair.alternIVs.size()>0) a_pairPoints.push_back(pair);
	    }//if
	}
    }
    return a_pairPoints.size();
}

void
EBPetscCompGridRichards:: computeGradImplicitFunction(RealVect& a_grad,
						      const RealVect& a_point,
						      const int a_dim,
						      const Real a_delta)
{
    CH_TIME("EBPetscCompGridRichards::computeGradImplicitFunction");
    CH_assert(a_dim > 0 && a_dim <= 3);
    //Real valPoint = m_implicitBaseIF->value(a_point);
    for (int idir=0;idir<a_dim;idir++)
    {
	RealVect pointDeltaPlus  = a_point;
	RealVect pointDeltaMinus = a_point;
	pointDeltaPlus [idir] +=a_delta;
	pointDeltaMinus[idir] -=a_delta;
	a_grad[idir]=2.*(m_implicitBaseIF->value(pointDeltaPlus)-m_implicitBaseIF->value(pointDeltaMinus))/a_delta;
    }
}

int  
EBPetscCompGridRichards:: findIntersectionPoints(IFData2<3>& a_ifdata,
						 Vector<RealVect>& a_points,
						 Vector<IntVect>&  a_sides,
						 const VolIndex& a_vof,
						 const int a_lev,
						 const Vector<int>& a_refRatios,
						 const bool a_skipDuplicates)
{
    CH_TIME("EBPetscCompGridRichards::findIntersectionPoints");
    IndexTM<Real, 3> vectdx;
    Real dx=1.0;
    //add these lines for underresolved cells call
    // Vector<int> ratios=m_refRatios;
    //ratios[m_finestLevel]=2;
    //for (int ilev = m_finestLevel;ilev < a_lev; ++ilev)	ratios.push_back(2);
    if (a_lev > 0)
    {
	for (int ilev = 0;ilev < a_lev; ++ilev)
	{
	    dx /= a_refRatios[ilev];
	}
    }
    
    vectdx.setAll(dx);

    //RealVect vofLocation = RealVect(a_vof.gridIndex())+0.5*RealVect::Unit;
    RealVect vofLocation = EBArith::getVoFLocation(a_vof, dx*RealVect::Unit, RealVect::Zero);
    //RealVect physLocation = RealVect(a_vof.gridIndex())*dx;
    //const RealVect centroid = (vofLocation+a_ebisBox.centroid(a_vof))*m_vectDxs[a_lev];

    RvgDim cellCenter;
    int degreeP = 2;
    for (int idir = 0;idir < 3; ++idir)
    {
	cellCenter[idir] = vofLocation[idir];
	//cellCenter[idir] = physLocation[idir];
    }
    //pout()<<a_vof<<"\t"<<vofLocation<<"\t"<< vectdx<<endl;

    a_ifdata.defineCoordinateSystem(vectdx, cellCenter, degreeP);
    const EdgeIntersections& edgeIn= a_ifdata.m_intersections;
    //ifdata.print(pout());
    for (typename EdgeIntersections::const_iterator it = edgeIn.begin();
	 it != edgeIn.end(); ++it)
    {
#ifdef UGDEBUGOUT
	pout() << a_vof<<"\t"<< "\t EdgeIndex "
	       << it->first
	       << " = "
	       << setprecision(5)
	    // << setiosflags(ios::showpoint)
	    // << setiosflags(ios::scientific)
	       << it->second
	       << "\n";
#endif
	RealVect point;
	IntVect  side = IntVect::Zero;
	//point[2]=0.;
	const Real half = 0.5*dx;
	const Real halfDx = 0.5*dx;
	const Real itSecond = it->second; 
	if (it->first[0] == 0)
	{
	    point[0]=vofLocation[0]+itSecond;
	    if (itSecond == -halfDx)
	    {
		side[0]=-1;
	    }
	    else if (itSecond == halfDx)
	    {
		side[0]=1;
	    }
	    
	    if (it->first[1] == 0)
	    {
		point[1] = vofLocation[1]-half;
		side[1]  = -1; 
	    }else
	    {
		point[1]=vofLocation[1]+half;
		side[1]  = 1; 
	    }
	    
	    if (it->first[2] == 0)
	    {
		side[2] = -1;
		point[2]=vofLocation[2]-half;

	    }else
	    {
		side[2] = 1;
		point[2]=vofLocation[2]+half;
	    }
	    
	}
	else if (it->first[0] == 1)
	{
	    point[1]=vofLocation[1]+itSecond;
	    if (itSecond == -halfDx)
	    {
		side[1]=-1;
	    }
	    else if (itSecond == halfDx)
	    {
		side[1]=1;
	    }
		    
	    if (it->first[1] == 0)
	    {
		point[0]=vofLocation[0]-half;
		side[0]  = -1; 
	    }else
	    {
	        point[0]=vofLocation[0]+half;
		side[0]  = 1; 
	    }
	    
	    if (it->first[2] == 0)
	    {
		side[2] = -1;
		point[2]=vofLocation[2]-half;
	    }else
	    {
		side[2] = 1;
		point[2]=vofLocation[2]+half;
	    }
	    
	}else
	{
	    point[2]=vofLocation[2]+itSecond;
	    if (it->first[1] == 0)
	    {
		point[0]=vofLocation[0]-half;
		side[0]  = -1; 
	    }else
	    {
		point[0]=vofLocation[0]+half;
		side[0]  = 1; 
	    }
	    if (it->first[2] == 0)
	    {
		point[1]=vofLocation[1]-half;
		side[1]  = -1; 
	    }else
	    {
		point[1]=vofLocation[1]+half;
		side[1]  = 1; 
	    }
	    
	    side[2] = 0;

	    if (itSecond == -halfDx)
	    {
		side[2]=-1;
	    }
	    else if (itSecond == halfDx)
	    {
		side[2]=1;
	    }
	}

	//check if point belongs to this vof
	//const RealVect phys_point =point*m_vectDxs[a_lev];
	// const bool checkIn = checkPointBelongToVoF(phys_point, centroid, m_vectDxs[a_lev]);
	// if (!checkIn) continue;

	//pout() << "\t" << point<<"\n";
	//point *=vectdx;
	bool notThere=true;
	int size=a_points.size();
	if (a_skipDuplicates)
	{
	    for (int ipoint=0;ipoint<size;ipoint++)
	    {
		if(a_points[ipoint]==point)
		{
		    // //we need correct sides
		    // for (int idir = 0; idir < 3; idir++)
		    // {
		    // 	if(!a_sides[ipoint][idir] && side[idir]) a_sides[ipoint][idir]=side[idir];
		    // }
		    notThere = false;
		    break;
		}
	    }
	}
	if (notThere)
	{
	    a_points.push_back(point);
	    a_sides.push_back(side);
	}
    }
    // pout() << "\t" << "\n";

    return a_points.size();
}

void
EBPetscCompGridRichards::define2DGradStencil(const int a_startLevel, 
					     const int a_finestLevel)
{
    CH_TIME("EBPetscCompGridRichards::define2DGradStencil");
    //const int nghostTemp = 2;

    int startLevel = 0;
    // if (a_startLevel > startLevel)
    // {
    // 	startLevel=a_startLevel;
    // }
    m_finestLevel = a_finestLevel;
     //gets VofStencil for Irr surface
    for (int ilev = startLevel; ilev <= a_finestLevel; ilev++)
    {
	//const ProblemDomain& domain  = m_eblgsPtr[ilev]->getDomain();
	//const IntVect domainSize = domain.size();
	const RealVect& vectDx = m_vectDxs[ilev];
	const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
//// #pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& datInd = dit[mybox];
	    const EBISBox&   ebisBox = ebisl[datInd];

	    // Box grid;
	    // if (ilev)
	    // {
	    // 	grid = grow(dbl.get(datInd), nghostTemp);
	    // 	grid &= domain;
	    // }
	    // else
	    // {

	    // 	grid = dbl.get(datInd);
	    // }
	    const Box& 	grid = dbl.get(datInd);

	    IntVectSet irregSet =  ebisBox.getIrregIVS(grid);
	    
	    for (VoFIterator vofit(irregSet,ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
		const VolIndex& vof = vofit();
		
		if (ebisBox.bndryArea(vof)<Tol) continue;
		const RealVect vofLoc     = EBArith::getVoFLocation(vof, vectDx, RealVect::Zero) + ebisBox.bndryCentroid(vof)*vectDx;

		// VoFStencil& xStencil = (*m_2DGradStencil[ilev])[datInd](vof,0);
		// VoFStencil& yStencil = (*m_2DGradStencil[ilev])[datInd](vof,1);
		Vector<VolIndex> vecVoFs[2];
		Vector<VolIndex> closestVoFs;
		Real xcoords[3];
		Real ycoords[3];
		Real weights[4];
		    
		//we need 3 additional vofs for bilinear stencil derivative
		const Vector<UGData>& ugDataVec = (*m_surfaceGridData[ilev])[datInd](vof,0); //this is a copy
		const int ugSize =  ugDataVec.size();

		//fill containers
		const int minSize = std::min(3,ugSize);
		if (ugSize<=3)
		{
		    for (int idata=0; idata<minSize; idata++)
		    {
			const VolIndex& conVoF = ugDataVec[idata].vof;
			if (!conVoF.isDefined()) continue;
			for (int idir=0;idir<2;idir++)
			{
			    vecVoFs[idir].push_back(conVoF);
			}
		    }
		}
		else
		{
		    //fill with closest vof in coordinate direction in case there are more than 3 vofs in UGData
		    int skipIndx[6] = {-1,-1,-1,-1,-1,-1};
		    for (int ivof=0;ivof<minSize;ivof++)
		    {
			Real mindist[2] = {1e+20, 1e+20};
			int minIndx[2]  = {-1, -1};
			for (int idata=0; idata<ugSize; idata++)
			{
			    const VolIndex& conVoF = ugDataVec[idata].vof;
			    if (!conVoF.isDefined()) continue;

			    const RealVect conVoFLoc = EBArith::getVoFLocation(conVoF, vectDx, RealVect::Zero) + ebisBox.bndryCentroid(conVoF)*vectDx;
			    const RealVect dist = conVoFLoc-vofLoc;
			    for (int idir=0;idir<2;idir++)
			    {
				bool skip = false;
				for (int iskip=0;iskip<3;iskip++)
				{
				    if (skipIndx[iskip*2+idir] == idata) skip = true;
				}
				if (skip) continue;
				    
				if (fabs(dist[idir])>Tol && fabs(dist[idir])<mindist[idir])
				{
				    mindist[idir]=fabs(dist[idir]);
				    minIndx[idir]=idata;
				}
			    }
			}
			for (int idir=0;idir<2;idir++)
			{
			    if (minIndx[idir]>=0)
			    {
				vecVoFs[idir].push_back(ugDataVec[minIndx[idir]].vof);
				skipIndx[ivof*2+idir]=minIndx[idir];
			    }
			}
		    }
		    //stenSize=vecVoFs[0].size();
		    //CH_assert(vecVoFs[0].size()==vecVoFs[1].size());
		}
		
		// pout()<<vof<<"\t size="<<ugSize<<"\t stenSizeX="<<vecVoFs[0].size()<<"\t stenSizeY="<<vecVoFs[1].size()<<endl;
		// for (int idir=0;idir<2;idir++)
		// {
		//     pout()<<"  idir="<<idir;
		//     for (int ivof=0;ivof<vecVoFs[idir].size();ivof++)
		//     {
		// 	pout()<<"\t"<<vecVoFs[idir][ivof];
		//     }
		//     pout()<<endl;
		// }
		//we must select from closest vofs to get 4 vofs
		if (vecVoFs[0].size()<3 || vecVoFs[1].size()<3)
		{
		    findClosestIrrVoFs(closestVoFs,vof,vectDx,ebisBox, 3, 2);
		    
		    //fill the containers with closest vof not in the stencil
		    for (int idir=0;idir<2;idir++)
		    {
			if (vecVoFs[idir].size()>=3) continue;
			int vecSize = vecVoFs[idir].size();
			const int closeSize=closestVoFs.size();
			for (int ivof=0;ivof<closeSize;ivof++)
			{
			    bool isIn = false;
			    for (int idata=0; idata<vecSize; idata++)
			    {
				if (vecVoFs[idir][idata] == closestVoFs[ivof])
				{
				    isIn = true;
				    break;
				}
			    }
			    if (!isIn)
			    {
				vecVoFs[idir].push_back(closestVoFs[ivof]);
				vecSize++;
				if (vecSize>=3) break;
			    }
			}
		    }
		}
		
		CH_assert(vecVoFs[0].size()==3);
		CH_assert(vecVoFs[1].size()==3);

                //create xstencil
		for (int idir=0;idir<2;idir++)
		{
		    VoFStencil& curStencil = (*m_2DGradStencil[ilev])[datInd](vof,idir);
		    for (int ivof=0;ivof<3;ivof++)
		    {
			const VolIndex& curVoF=vecVoFs[idir][ivof];
			const RealVect conVoFLoc = EBArith::getVoFLocation(curVoF, vectDx, RealVect::Zero) + ebisBox.bndryCentroid(curVoF)*vectDx;
			const RealVect coord = conVoFLoc - vofLoc;
			xcoords[ivof] = coord[0];
			ycoords[ivof] = coord[1];
		    }
		    scatt_bilinear_derivative_stencil(weights, xcoords, ycoords, idir);
		    if (fabs(weights[0])>Tol) curStencil.add(vof, weights[0]);
		
		    for (int ivof=0;ivof<3;ivof++)
		    {
			const VolIndex& curVoF=vecVoFs[idir][ivof];
			if (fabs(weights[ivof+1])>Tol) curStencil.add(curVoF, weights[ivof+1]);
		    }
		    
		    // for (int ivof=0;ivof<curStencil.size();ivof++)
		    // {
		    // 	pout()<<"\t idir="<<idir<<"\t "<<curStencil.vof(ivof)<<"\t w="<<curStencil.weight(ivof)<<endl;
		    // }
		}
	    }//vof
	}//mybox
    }
}

void  
EBPetscCompGridRichards::printOutputRate(const Vector<LevelData<EBCellFAB>* >& a_psi, 
					 const Real a_time,
					 const Real a_dt)
{
    CH_TIME("EBPetscCompGridRichards::printOutputRate");
    // const int rank = procID();
    // pout<<"rank="<<rank<<endl;
    //const Real Tol = 1.e-15;

    Real output=0.0, output_real=0.0;
    Real upslopeX;
    //int ilev=0;
    const int ilev = m_finestLevel;
    //for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	const Real dy = m_vectDxs[ilev][1];

	const int ivEnd = m_eblgsPtr[ilev]->getDomain().size(0)-1;
	upslopeX = (1.0+ivEnd)*m_vectDxs[ilev][0];
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();

	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
// #pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& datInd = dit[mybox];
	    const Box& grid = dbl.get(datInd);
	    IntVectSet ivsIrreg = ebisl[datInd].getIrregIVS(grid);
	    if(ivsIrreg.isEmpty()) continue;

	    IntVect small = IntVect::Zero;
	    small[1] =  grid.smallEnd(1);
	    IntVect big = IntVect::Zero;
	    big[1]=grid.bigEnd(1);
	    Box sliceBox(small, big);
	    FArrayBox bndryValue(sliceBox, 4);
	    bndryValue.setVal(0.0);
	    //bndryValue.setVal(-1.0,3);

	    IntVectSet ivsGrid;

	    //  Real max_psi=0.0;//test it
	    
	    const EBISBox&   ebisBox = ebisl[datInd];
	    const EBGraph& ebgraph = ebisl[datInd].getEBGraph();
	    for (VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
	    {
		const VolIndex& vof = vofit();
		const IntVect& iv = vof.gridIndex();
	 	const Real psi = (*a_psi[ilev])[datInd](vof,0);
		
		RealVect vofCenterLoc  = EBArith::getVoFLocation(vof, m_vectDxs[ilev], RealVect::Zero);
		vofCenterLoc          += ebisBox.bndryCentroid(vof)*m_vectDxs[ilev];		//get min upslope X
		if (psi>=0. && ebisBox.bndryArea(vof))
		{
		  
		  if (vofCenterLoc[0]<upslopeX)
		    {
		      upslopeX=vofCenterLoc[0];
		    }
		  
		}
		//pout()<<vof<<"\t psi="<<psi<<"\t x="<<vofCenterLoc<<"\t min="<<upslopeX<<"\t bndrycentroid="<<ebisBox.bndryCentroid(vof)<<"\t area="<<ebisBox.bndryArea(vof)<<"\t n="<<ebisBox.normal(vof)<<"\t k="<<ebisBox.volFrac(vof)<<endl;

		if (iv[0] != ivEnd) continue;
		const Vector<UGData>& ugDataVec = (*m_surfaceGridData[ilev])[datInd](vof,0);
		const int ugVecSize =  ugDataVec.size();
		
		Real xSlope = (*m_domainSlopes[ilev])[datInd](vof,0);
		const Real nmann =  (*m_manningCoeff[ilev])[datInd](vof,0);
		const RealVect normal0 = -ebisBox.normal(vof);
		const Real area = ebisBox.bndryArea(vof)*fabs(normal0[2]);
		IntVect ivSlice(0,iv[1],0);

  		
		//get an average
		{   
		    bndryValue(ivSlice,0) += psi*area;
		    bndryValue(ivSlice,1) += xSlope*area;
		    bndryValue(ivSlice,2) += nmann*area;
		    bndryValue(ivSlice,3) += area;//ebisBox.bndryCentroid(vof)[0];
		    ivsGrid |= ivSlice;
		    
		    for (int iug=0; iug<ugVecSize; ++iug)
		    {
			if (ugDataVec[iug].slopes.size()) xSlope=ugDataVec[iug].slopes[0];
			if (m_surface_solver_type==DiffusionWaveUnstructured) xSlope += (*m_2DGrad[ilev])[datInd](vof,0);

			//for (int iface=0; iface<ugDataVec[iug].numFaces;iface++)
			{
			    const int& bndryDir = ugDataVec[iug].bndryDir;//s[iface];
		    
			    if (psi>0.0 && bndryDir==1)
			    {
				const Real& length  = ugDataVec[iug].faceData[2];//[iface][2];
				output_real += pow(psi,5./3.)*sqrt(fabs(xSlope))/nmann*length;
			    }
			}
		    }
		}
	    }
	    for (IVSIterator ivit(ivsGrid); ivit.ok(); ++ivit)
	    {
		const IntVect& iv = ivit();
		for (int jj=0; jj<4; ++jj) bndryValue(iv,jj) /= bndryValue(iv,3);
		if (bndryValue(iv,0) >0.0) output += pow(bndryValue(iv,0),5./3.)*sqrt(fabs(bndryValue(iv,1)))/bndryValue(iv,2)*dy;
	    }
	    
	}
    }
    Real TotalOutput=output; 
    Real TotalOutput2=output_real; 
    //pout()<<output<<"\t Realoutput=";
    Real minupslopeX=upslopeX;

#ifdef CH_MPI
    const int rank = procID();
    //  pout()<<"rank="<<rank<<endl;

    MPI_Reduce(&output, &TotalOutput, 1, MPI_CH_REAL, MPI_SUM, 0, Chombo_MPI::comm);
    MPI_Reduce(&output_real, &TotalOutput2, 1, MPI_CH_REAL, MPI_SUM, 0, Chombo_MPI::comm);
    MPI_Reduce(&upslopeX, &minupslopeX, 1, MPI_CH_REAL, MPI_MIN, 0, Chombo_MPI::comm);
    
    if (!rank)
#endif
    {
	pout()<<TotalOutput<<"\t RealOutput="<<TotalOutput2<<" \t X_Upslope="<<minupslopeX<<endl;
	static bool isOpened = false;
	std::ofstream ofs;
	if (!isOpened)
	{
	    ofs.open ("OutputRate.dat", std::ofstream::out);
	    isOpened = true;
	}
	else
	{
	    ofs.open ("OutputRate.dat", std::ofstream::out | std::ofstream::app);
	}
	ofs <<scientific;
	ofs.precision(14);
	ofs<<a_time+a_dt<<"\t"<<TotalOutput<<"\t"<<TotalOutput2<<"\t"<<minupslopeX<<endl;
	ofs.close();
    }
}

void  
EBPetscCompGridRichards::projectFluxOnXYPlane(Vector<RefCountedPtr<LevelData<BaseIVFAB<Real> > > >& a_val)
{
    CH_TIME("EBPetscCompGridRichards::projectFluxOnXYPlane");
    //const int nghostTemp = 1;
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
// #pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& datInd = dit[mybox];
	    const EBISBox&   ebisBox = ebisl[datInd];
	    //has to correspond to number of ghost cells
	    // Box grid;
	    // if (ilev)
	    // {
	    // 	grid = grow(dbl.get(datInd), nghostTemp);
	    // 	grid &= m_eblgsPtr[ilev]->getDomain();
	    // }
	    // else
	    // {
	    // 	grid = dbl.get(datInd);
	    // }
	    const Box& 	grid = dbl.get(datInd);
	    IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
	    for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
	    {
		const VolIndex& vof = vofit();
		//const RealVect normal0 = ebisBox.normal(vof);
		const RealVect normal = getAnisotropicNormal(ebisBox.normal(vof), m_vectDxs[ilev]);
		// Real anisotr_corr=0.0;
		// for (int idir = 0; idir < SpaceDim; idir++)
		// {
		//     anisotr_corr += normal0[idir]*normal0[idir] 
		// 	/(m_vectDxs[ilev][idir]*m_vectDxs[ilev][idir]);
		// }

		if (ebisBox.bndryArea(vof))
		    (*a_val[ilev])[datInd](vof, 0) *=fabs(normal[2]);///sqrt(anisotr_corr);
		else
		    (*a_val[ilev])[datInd](vof, 0) = 0.0;
		
		//pout()<<vof<<"\t"<<fabs(normal[2])<<"\t"<<(*a_val[ilev])[datInd](vof, 0)<<"\t"<<anisotr_corr<<endl;
	    }
	}
	a_val[ilev]->exchange();
    }
}

//computes exchange flux from surface wave equation
void
EBPetscCompGridRichards::addExchangeFlux(BaseIVFAB<Real>&       a_exflux,
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
					 const bool a_isDiffusionWave,
					 const bool a_useRegularization,
					 const Real a_regConstant)
{
    CH_TIME("EBPetscCompGridRichards::addExchangeFlux");
    const IntVect domainSize = a_domain.size();

    ParserFunc surfaceBCParser_Lo[2];
    ParserFunc surfaceBCParser_Hi[2];
    surfaceBCParser_Lo[0].define(m_paramExpr.surfaceBC[0]);
    surfaceBCParser_Lo[1].define(m_paramExpr.surfaceBC[1]);
    surfaceBCParser_Hi[0].define(m_paramExpr.surfaceBC[2]);
    surfaceBCParser_Hi[1].define(m_paramExpr.surfaceBC[3]);
    for (int idir=0; idir<2; idir++)
    {
	surfaceBCParser_Lo[idir].setGravitationalConstant(m_gconst);
	surfaceBCParser_Hi[idir].setGravitationalConstant(m_gconst);
    }

    IntVectSet irregSet = a_ebisBox.getIrregIVS(a_box);
    for (VoFIterator vofit(irregSet, a_ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
 	const VolIndex& vof = vofit();

	// if (vof.gridIndex()==IntVect(99,0,0))
	// {
	//     pout()<<endl;
	// }
	
	Real& exflux = a_exflux(vof, 0);

	Real psinew, psiold;
	// if (a_useRegularization)
	// {
	//     //adding regularization for surface max function
	//     if (a_psinew(vof, 0) > 0.0)
	//     {
	// 	psinew  = exp(-a_regConstant/a_psinew(vof, 0))*a_psinew(vof, 0);
	//     }
	//     else
	//     {
	// 	psinew = 0.0;
	//     }
	    
	//     if (a_psiold(vof, 0) > 0.0)
	//     {
	// 	psiold  = exp(-a_regConstant/a_psiold(vof, 0))*a_psiold(vof, 0);
	//     }
	//     else
	//     {
	// 	psiold = 0.0;
	//     }
	// }
	// else
	{
	    psinew  = std::max(a_psinew(vof, 0), 0.0);
	    psiold  = std::max(a_psiold(vof, 0), 0.0);
	}
	
	
	const Vector<UGData>& UGDataVec = a_UGData(vof,0);
	const int sizeUG = UGDataVec.size();
	if (sizeUG)
	{
	    Real flux[2] {0.,0.};
	    for (int idata=0; idata<sizeUG; idata++)
	    {
		const UGData& ug_point = UGDataVec[idata];
		const VolIndex& conVoF = ug_point.vof;
		//for (int iface=0; iface<ug_point.numFaces; iface++)
		{
		    const RealVect& normal = ug_point.faceData;//[iface];
		    const int bndryDir     = ug_point.bndryDir;//s[iface];    
		    if (bndryDir>=0)//!iv[idir] || iv[idir] == domainSize[idir]-1)
		    {
			Real fluxBC[2];
			ParserInput parserinput;
			RealVect centLoc  = EBArith::getVoFLocation(vof, a_vectDx, RealVect::Zero) + a_ebisBox.bndryCentroid(vof)*a_vectDx;
			// centLoc += EBArith::getVoFLocation(conVoF, a_vectDx, RealVect::Zero) + a_ebisBox.bndryCentroid(conVoF)*a_vepctDx;
			// centLoc /=2.0;

			parserinput.point = centLoc;
			parserinput.time  = m_time;
			parserinput.depth = 0.0; //no depth for surface here 
			parserinput.psi   = psinew;
			if (ug_point.slopes.size())
			{
			    parserinput.slope_x = ug_point.slopes[0];
			    parserinput.slope_y = ug_point.slopes[1];

			}
			else
			{
			    parserinput.slope_x = a_domSl(vof, 0);
			    parserinput.slope_y = a_domSl(vof, 1);
			}
			parserinput.manncoef = a_mannC(vof,0);
			//Here we assume that the origin is zero
			switch (bndryDir)
			{
			case 0:
			    parserinput.point[0] = 0.0;
			    parserinput.slope = parserinput.slope_x;
			    if (m_surface_solver_type==DiffusionWaveUnstructured) parserinput.slope += a_2DGrad(vof,0);
			    fluxBC[0] = surfaceBCParser_Lo[0].Eval(parserinput);
			    fluxBC[1] = 0.;
			    break;
			case 1:
			    parserinput.point[0] = domainSize[0]*a_vectDx[0];
			    parserinput.slope = parserinput.slope_x;
			    if (m_surface_solver_type==DiffusionWaveUnstructured) parserinput.slope += a_2DGrad(vof,0);
			    fluxBC[0] = surfaceBCParser_Hi[0].Eval(parserinput);
			    fluxBC[1] = 0.;
			    break;
			case 2:
			    parserinput.point[1] = 0.0;
			    parserinput.slope = parserinput.slope_y;
			    if (m_surface_solver_type==DiffusionWaveUnstructured) parserinput.slope += a_2DGrad(vof,1);
			    fluxBC[0] = 0.;
			    fluxBC[1] = surfaceBCParser_Lo[1].Eval(parserinput);
			    break;
			case 3:
			    parserinput.point[1] = domainSize[1]*a_vectDx[1];
			    parserinput.slope = parserinput.slope_y;
			    if (m_surface_solver_type==DiffusionWaveUnstructured) parserinput.slope += a_2DGrad(vof,1);
			    fluxBC[0] = 0.;
			    fluxBC[1] = surfaceBCParser_Hi[1].Eval(parserinput);
			    break;
			case 4:
			    fluxBC[0] = 0.;
			    fluxBC[1] = 0.;
			    break;
			    
			default:
			    fluxBC[0] = 0.;
			    fluxBC[1] = 0.;
			    break;
			}
			for (int idir=0; idir<2; idir++)
			{
			    flux[idir] -=fluxBC[idir]*normal[idir]*normal[2];
			}
		    }else
		    {
			//Real velc[2];
			//Real psi;
			for (int idir=0; idir<2; idir++)
			{
			    if (fabs(normal[idir])<Tol) continue;
			    
			    if (conVoF.isDefined())
			    {
				
				//velc = -(weight*a_domSl(vof, idir)+(1.0-weight)*a_domSl(conVoF, idir));
				// psi     = weight*psinew+(1.-weight)*std::max(a_psinew(conVoF,0), 0.0);
				// flux[idir] = copysign(sqrt(fabs(velc)), velc)*pow(psi,(5.0/3.0))/a_mannC(vof,0);
				Real velc=0;
				if (ug_point.slopes.size())
				    velc = ug_point.slopes[idir];
				else
				    velc = a_domSl(vof, idir);

				if (m_surface_solver_type==DiffusionWaveUnstructured) velc += a_2DGrad(vof,idir);
				velc *=normal[idir];
				// const IntVect& vofIV=conVoF.gridIndex();
				// pout()<< "isIrregular="<<a_ebisBox.isIrregular(vofIV)
				//       <<"\t isCovered="<<a_ebisBox.isCovered(vofIV)
				//       <<"\t isMultiValued="<<a_ebisBox.isMultiValued (vofIV)
				//       <<"\t bndryAreaEbis="<<a_ebisBox.bndryArea(conVoF)
				//       <<"\t normalEbis="<<a_ebisBox.normal(conVoF)
				//       <<"\t volFrac="<<a_ebisBox.volFrac(conVoF)
				//       <<"\t psiCon="<<a_psinew(conVoF,0)
				//       <<endl;

							
				Real velCon = a_domSl(conVoF, idir);

				if (m_surface_solver_type==DiffusionWaveUnstructured) velCon += a_2DGrad(conVoF,idir);
				velCon *= normal[idir];
				
				const Real psiCon   = std::max(a_psinew(conVoF,0), 0.0);
				if (velc >=0. && velCon > 0.) 
				{
				    if (psiCon>0.0)
					flux[idir] += sgn(velCon)*sqrt(fabs(velCon))*pow(psiCon,(5.0/3.0))/a_mannC(conVoF,0)*normal[2];
				    
				}
				
				if (velc<=0. && velCon <= 0.) 
				{
				    if (psinew>0.0)
					flux[idir] += sgn(velc)*sqrt(fabs(velc))*pow(psinew,(5.0/3.0))/a_mannC(vof,0)*normal[2];
				}
			    }
			    else
			    {
				MayDay::Error("addExchangeFlux: Undefined conVoF");

				// if (psinew>0.0)
				// 	flux[idir] = copysign(sqrt(fabs(a_domSl(vof, idir))), -a_domSl(vof, idir))*pow(psinew,(5.0/3.0))/a_mannC(vof,0);
				// else
				// 	flux[idir] = 0.0;
			    }
			}
		    }
		}//iface
	    }//idata
	    
	    for (int idir=0; idir<2; idir++)
	    {
		//const Real flux = copysign(sqrt(fabs(velc[idir])), velc[idir])*pow(psi,(5.0/3.0))/a_mannC(vof,0);
		exflux += flux[idir]/a_bndryArea(vof,0);
	    }
	}//sizeUG
	exflux -= (psinew -psiold)/a_dt;
    }//vofIterator
}



void  
EBPetscCompGridRichards::floorIrregular(Vector<LevelData<BaseIVFAB<Real> >* >& a_psi)
{
    CH_TIME("EBPetscCompGridRichards::floorIrregular");
    //const int nghostTemp = 2;
    for (int ilev = 0; ilev <= m_finestLevel; ilev++)
    {
	const DisjointBoxLayout& dbl = m_eblgsPtr[ilev]->getDBL();
	const EBISLayout& ebisl = m_eblgsPtr[ilev]->getEBISL();
	DataIterator dit = dbl.dataIterator();
	int nbox=dit.size();
#pragma omp parallel for
	for (int mybox=0;mybox<nbox; mybox++)
	{
	    const DataIndex& datInd = dit[mybox];
	    //has to correspond to number of ghost cells
	    // Box grownBox = grow(dbl.get(datInd), nghostTemp);
	    // grownBox &= m_eblgsPtr[ilev]->getDomain();
	    const Box& grownBox = dbl.get(datInd);
	    IntVectSet ivsIrreg = ebisl[datInd].getIrregIVS(grownBox);
	    for (VoFIterator vofit(ivsIrreg, ebisl[datInd].getEBGraph()); vofit.ok(); ++vofit)
	    {
		const VolIndex& vof = vofit();
		if ((*a_psi[ilev])[datInd](vof, 0) < 0.0) 
		    (*a_psi[ilev])[datInd](vof, 0) = 0.0;
	    }
	}
	a_psi[ilev]->exchange();
      }
}

#include "NamespaceFooter.H"

