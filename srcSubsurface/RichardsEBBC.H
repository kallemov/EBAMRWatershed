#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _RICHARDSEBBC_H_
#define _RICHARDSEBBC_H_

#include "RefCountedPtr.H"

#include "NeumannPoissonEBBC.H"
#include "BaseBCValue.H"
#include "NamespaceHeader.H"

///
/**
 */
class RichardsEBBC: public ConductivityBaseEBBC
{
public:
  ///
  virtual void define(const LayoutData<IntVectSet>& a_cfivs,
                      const Real&                   a_factor)
  {
    ;
    //no flux stencil for Neumann
  }

  ///
  RichardsEBBC(const ProblemDomain& a_domain,
	       const EBISLayout&    a_layout,
	       const RealVect&      a_vectDx);
  
  virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
  {
    return NULL;
  }

  ///
  virtual ~RichardsEBBC();

  ///
  virtual void setValue(Real a_value);

  ///
  virtual int whichBC()
  {
    return 0;
  };

  ///
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

  //include overland solver
  virtual void setOverlandFlowSolver(bool a_OverlandFlowSolver)
  {
    m_OverlandFlowSolver = a_OverlandFlowSolver;
  }

  ///
  virtual void getEBFlux(Real&                         a_flux,
                         const VolIndex&               a_vof,
                         const LevelData<EBCellFAB>&   a_phi,
                         const LayoutData<IntVectSet>& a_cfivs,
                         const DataIndex&              a_dit,
                         const RealVect&               a_probLo,
                         const RealVect&               a_vectDx,
                         const bool&                   a_useHomogeneous,
                         const Real&                   a_time,
                         const pair<int,Real>*         a_cacheHint=0 );

  virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                           const EBCellFAB&              a_phi,
                           VoFIterator&                  a_vofit,
                           const LayoutData<IntVectSet>& a_cfivs,
                           const DataIndex&              a_dit,
                           const RealVect&               a_probLo,
                           const RealVect&               a_vectDx,
                           const Real&                   a_factor,
                           const bool&                   a_useHomogeneous,
                           const Real&                   a_time);


private:
  NeumannPoissonEBBC m_bc;
  bool m_OverlandFlowSolver;
};

///
class RichardsEBBCFactory: public BaseEBBCFactory
{
public:
  ///
  RichardsEBBCFactory();

  ///
  virtual ~RichardsEBBCFactory();

  ///
  virtual void setValue(Real a_value);

  ///
  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

  ///
  virtual void setData( RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_data)
  {
    m_data = a_data;
    m_dataBased = true;
  }

  //include overland solver
  void setOverlandFlowSolver(bool a_OverlandFlowSolver)
  {
    m_OverlandFlowSolver = a_OverlandFlowSolver;
  }

  ///
  virtual RichardsEBBC* create(const ProblemDomain& a_domain,
			       const EBISLayout&    a_layout,
			       const RealVect&      a_dx,
			       const IntVect*       a_ghostCellsPhi=0,
			       const IntVect*       a_ghostCellsRhs=0);

private:
  bool m_onlyHomogeneous;
  bool m_isFunction;
  RefCountedPtr<LevelData<BaseIVFAB<Real> > >m_data;
  bool m_dataBased;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_flux;
  bool m_OverlandFlowSolver;
};

#include "NamespaceFooter.H"
#endif
