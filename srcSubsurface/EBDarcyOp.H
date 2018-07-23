#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBDARCYOP_H_
#define _EBDARCYOP_H_

#include "EBAMRPoissonOp.H"
#include "EBLevelDataOps.H"
#include "RichardsEBBC.H"
#include "RichardsDomainBC.H"

#include "NamespaceHeader.H"

#define EBAMRPOVC_NUMSTEN 8

//! \class EBDarcyOp
//! This class implements an operator that solves the equation
//! the class is a copy paste copy of ConductivityOp class with changes in flux terms and BC functions 
//! (beta div (b grad) (phi+z) = rhs
//! using the AMRLevelOp interface.
class EBDarcyOp
{
  public:

  //! Constructs a darcy operator similar to conductivity operator
  EBDarcyOp(const EBLevelGrid &                                  a_eblgFine,
	    const EBLevelGrid &                                  a_eblg,
	    const EBLevelGrid &                                  a_eblgCoar,
	    const RefCountedPtr<RichardsDomainBC>&               a_domainBC,
	    const RefCountedPtr<RichardsEBBC>&                   a_ebBC,
	    const RefCountedPtr<EBQuadCFInterp>&                 a_quadCFI,
	    const int&                                           a_refToFine,
	    const int&                                           a_refToCoar,
	    const RealVect&                                      a_vectDx,
	    const Real&                                          a_beta,
	    const RefCountedPtr<LevelData<EBFluxFAB> >&          a_bcoef,
	    const RefCountedPtr<LevelData<BaseIVFAB<Real> > >&   a_bcoIrreg,
	    const RefCountedPtr<LevelData<EBFluxFAB> >&          a_darcyFlux,
	    const IntVect&                                       a_ghostCellsPhi,
	    const IntVect&                                       a_ghostCellsRHS);
  

  //! Destructor.
  virtual ~EBDarcyOp();


  //only weights by kappa.  time and tga have their demands.
  virtual void kappaScale(LevelData<EBCellFAB> & a_rhs);

 
  void setBeta(const Real& a_beta);

  void diagonalScale(LevelData<EBCellFAB> & a_rhs,
                     bool a_kappaWeighted);

  
  void fillGrad(const LevelData<EBCellFAB>& a_phi);

  void getFlux(EBFluxFAB&                    a_flux,
               const LevelData<EBCellFAB>&   a_data,
               const Box&                    a_grid,
               const DataIndex&              a_dit,
               Real                          a_scale);

  void getFlux(EBFaceFAB&                    a_fluxCentroid,
               const EBCellFAB&              a_phi,
               const Box&                    a_ghostedBox,
               const Box&                    a_fabBox,
               const ProblemDomain&          a_domain,
               const EBISBox&                a_ebisBox,
               const RealVect&               a_vectDx,
               const DataIndex&              a_datInd,
               const int&                    a_idir);

  void getFluxIrregular(EBCellFAB&                  a_lhs,
			const LevelData<EBCellFAB>& a_phi,
			const DataIndex&            a_datInd);
  // This operator may be time-dependent for BC.
  void setTime(Real a_Time);

  ///
  /**
     This function assumes that coarse-fine boundary condtions have
     been dealt with.
  */
  virtual void applyOp(LevelData<EBCellFAB>&             a_opPhi,
                       const LevelData<EBCellFAB>&       a_phi,
                       const LevelData<EBCellFAB>* const a_phiCoarse,
                       const bool&                       a_homogeneousPhysBC,
                       const bool&                       a_homogeneousCFBC);

  /// virtual function called by LevelTGA
  virtual void applyOpNoBoundary(LevelData<EBCellFAB>&        a_opPhi,
                                 const LevelData<EBCellFAB>&  a_phi);

  ///
  /**
     this is the linearop function.  CFBC is set to homogeneous.  phic is null
  */
  virtual void applyOp(LevelData<EBCellFAB>&             a_opPhi,
                       const LevelData<EBCellFAB>&       a_phi,
                       bool                              a_homogeneousPhysBC);

  ///


  virtual void applyAMROp(LevelData<EBCellFAB>&       a_LofPhi,
			  const LevelData<EBCellFAB>& a_phiFine,
			  const LevelData<EBCellFAB>& a_phi,
			  const LevelData<EBCellFAB>& a_phiCoar,
			  const bool                  a_homogeneousPhysBC,
			  EBDarcyOp*                  a_finerOp);

  /**
   */
  virtual void create(LevelData<EBCellFAB>&       a_lhs,
                      const LevelData<EBCellFAB>& a_rhs);


  ///
  /**
   */
  virtual void assign(LevelData<EBCellFAB>&       a_lhs,
                      const LevelData<EBCellFAB>& a_rhs);

  ///
  /**
   */
  virtual Real dotProduct(const LevelData<EBCellFAB>& a_1,
                          const LevelData<EBCellFAB>& a_2);

  ///
  /**
   */
  virtual void incr(LevelData<EBCellFAB>&       a_lhs,
                    const LevelData<EBCellFAB>& a_x,
                    Real                        a_scale);

  ///
  /**
   */
  virtual void axby(LevelData<EBCellFAB>&       a_lhs,
                    const LevelData<EBCellFAB>& a_x,
                    const LevelData<EBCellFAB>& a_y,
                    Real                        a_a,
                    Real                        a_b);

  ///
  /**
   */
  virtual void scale(LevelData<EBCellFAB>& a_lhs,
                     const Real&           a_scale);

  ///
  /**
   */
  virtual Real norm(const LevelData<EBCellFAB>& a_rhs,
                    int                         a_ord);

  virtual Real localMaxNorm(const LevelData<EBCellFAB>& a_rhs);
  ///
  /**
   */
  virtual void setToZero(LevelData<EBCellFAB>& a_lhs);

  ///
  /**
   */
  virtual void setVal(LevelData<EBCellFAB>& a_lhs, const Real& a_value);


  void getDivFStencil(VoFStencil&      a_vofStencil,
                      const VolIndex&  a_vof,
                      const DataIndex& a_dit);

  void getJacobianStencil(VoFStencil&      a_vofStencil, //stencil for pressurehead derivative
			  const VolIndex&  a_vof,
			  const DataIndex& a_dit);

  void getFluxStencil(VoFStencil&      a_fluxStencil,
                      const FaceIndex& a_face,
                      const DataIndex& a_dit,
		      const int a_dir);

  void getFaceCenteredFluxStencil(VoFStencil&      a_fluxStencil,
                                  const FaceIndex& a_face,
                                  const DataIndex& a_dit,
				  const int a_dir);

  void incrOpRegularDir(EBCellFAB&             a_lhs,
                        const EBCellFAB&       a_phi,
                        const bool&            a_homogeneous,
                        const int&             a_dir,
                        const DataIndex&       a_datInd);
  void applyOpIrregular(EBCellFAB&             a_lhs,
                        const EBCellFAB&       a_phi,
                        const bool&            a_homogeneous,
                        const DataIndex&       a_datInd);


  ///do not call this one unless you really know what you are doing
  void defineStencils();
protected:
  void incrOpRegularAllDirs(Box * a_loBox,
                            Box * a_hiBox,
                            int * a_hasLo,
                            int * a_hasHi,
                            Box & a_curDblBox,
                            Box & a_curPhiBox,
                            int a_nComps,
                            BaseFab<Real> & a_curOpPhiFAB,
                            const BaseFab<Real> & a_curPhiFAB,
                            bool a_homogeneousPhysBC,
                            const DataIndex& a_dit);

  void applyDomainFlux(Box * a_loBox,
                       Box * a_hiBox,
                       int * a_hasLo,
                       int * a_hasHi,
                       Box & a_dblBox,
                       int a_nComps,
                       BaseFab<Real> & a_phiFAB,
                       bool a_homogeneousPhysBC,
                       const DataIndex& a_dit);


  static bool                     s_turnOffBCs;
    //static bool                     s_forceNoEBCF;

  void calculateAlphaWeight();

  void defineEBCFStencils();

  void getFluxEBCF(EBFaceFAB&                    a_flux,
                   const EBCellFAB&              a_phi,
                   const Box&                    a_ghostedBox,
                   Vector<FaceIndex>&            a_faceitEBCF,
                   Vector<VoFStencil>&           a_stenEBCF);

  void getFluxRegOnly(EBFaceFAB&                    a_fluxCentroid,
                      const EBCellFAB&              a_phi,
                      const Box&                    a_ghostedBox,
                      const RealVect&               a_vectDx,
                      const DataIndex&              a_datInd,
                      const int&                    a_idir);


  const IntVect                   m_ghostCellsPhi;
  const IntVect                   m_ghostCellsRHS;
  RefCountedPtr<EBQuadCFInterp>   m_quadCFIWithCoar;


  EBLevelGrid                     m_eblg;
  EBLevelGrid                     m_eblgFine;
  EBLevelGrid                     m_eblgCoar;

  RefCountedPtr<RichardsDomainBC>     m_domainBC;
  RefCountedPtr<RichardsEBBC>         m_ebBC;

  RealVect                        m_vectDx;
  RealVect                        m_vectDxFine;
  RealVect                        m_vectDxCoar;

  Real                            m_time;


  RefCountedPtr<LevelData<EBFluxFAB> >          m_bcoef;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >   m_bcoIrreg;
  RefCountedPtr<LevelData<EBFluxFAB> >          m_DarcyFlux;

  Real                            m_beta;
  //weights that get multiplied by alpha
  LayoutData<BaseIVFAB<Real> >       m_alphaDiagWeight;
  // //weights that get multiplied by beta
  // LayoutData<BaseIVFAB<Real> >       m_betaDiagWeight;
  bool                            m_hasEBCF;

  //stencils for operator evaluation
  LayoutData<RefCountedPtr<EBStencil> >  m_opEBStencil;


  //cache the vofiterators
  //for irregular cell iteration (includes buffer around multivalued cells)
  LayoutData<VoFIterator >                     m_vofIterIrreg;
  LayoutData<VoFIterator >                     m_vofIterMulti;
  //for domain boundary conditions at ir regular cells
  LayoutData<VoFIterator >                     m_vofIterDomLo[CH_SPACEDIM];
  LayoutData<VoFIterator >                     m_vofIterDomHi[CH_SPACEDIM];


  //stuff to make EBCF go faster
  LayoutData< Vector<FaceIndex>  >   m_faceitCoar[2*SpaceDim];
  LayoutData< Vector<VoFStencil> >  m_stencilCoar[2*SpaceDim];

  //flux register with finer level
  EBFastFR       m_fastFR;

  int                             m_refToFine;
  int                             m_refToCoar;
  bool                            m_hasFine;
  bool                            m_hasCoar;


  // Vector<IntVect> m_colors;

private:

  void reflux(LevelData<EBCellFAB>& a_residual,
              const LevelData<EBCellFAB>& a_phiFine,
              const LevelData<EBCellFAB>& a_phi,
              EBDarcyOp*                  a_finerOp);

  void getFlux(FArrayBox&                    a_flux,
               const FArrayBox&              a_phi,
               const Box&                    a_faceBox,
               const int&                    a_idir,
               const RealVect&               a_vectDx,
               const DataIndex&              a_datInd);

  void incrementFRCoar(EBFastFR&                   a_fluxReg,
                       const LevelData<EBCellFAB>& a_phiFine,
                       const LevelData<EBCellFAB>& a_phi);
  
  void incrementFRFine(EBFastFR&                   a_fluxReg,
                       const LevelData<EBCellFAB>& a_phiFine,
                       const LevelData<EBCellFAB>& a_phi,
                       EBDarcyOp*                  a_finerOp);

private:

  //! Default constructor. Creates an undefined conductivity operator.
  EBDarcyOp();

  //copy constructor and operator= disallowed for all the usual reasons
  EBDarcyOp(const EBDarcyOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }

  void operator=(const EBDarcyOp& a_opin)
  {
    MayDay::Error("invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
