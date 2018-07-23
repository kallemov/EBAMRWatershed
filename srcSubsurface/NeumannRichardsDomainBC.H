#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEUMANNRICHARDSDOMAINBC_H_
#define _NEUMANNRICHARDSDOMAINBC_H_

#include "RefCountedPtr.H"

#include "NeumannPoissonDomainBC.H"
#include "ConductivityBaseDomainBC.H"
#include "NamespaceHeader.H"

class NeumannRichardsDomainBC: public ConductivityBaseDomainBC
{
public:
  NeumannRichardsDomainBC();

  virtual ~NeumannRichardsDomainBC();

  virtual void 
  fillPhiGhost(FArrayBox&     a_state,
               const Box&     a_valid,
               const Box&     a_domain,
               Real           a_dx,
               bool           a_homogeneous){;}


  virtual void setValue(Real a_value);

  virtual int whichBC(int                  a_idir,
                      Side::LoHiSide       a_side);

  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_dx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous);

  virtual void getFaceVel(Real&                 a_faceFlux,
                          const FaceIndex&      a_face,
                          const EBFluxFAB&      a_vel,
                          const RealVect&       a_probLo,
                          const RealVect&       a_dx,
                          const int&            a_idir,
                          const int&            a_icomp,
                          const Real&           a_time,
                          const Side::LoHiSide& a_side);

private:
  bool m_onlyHomogeneous;
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_flux;
};

class NeumannRichardsDomainBCFactory: public BaseDomainBCFactory
{
public:
  NeumannRichardsDomainBCFactory();

  virtual ~NeumannRichardsDomainBCFactory();

  virtual void setValue(Real a_value);

  virtual void setFunction(RefCountedPtr<BaseBCValue> a_flux);

  virtual NeumannRichardsDomainBC* create(const ProblemDomain& a_domain,
                                              const EBISLayout&    a_layout,
                                              const RealVect&      a_dx);

private:
  bool m_onlyHomogeneous;
  bool m_isFunction;

  Real m_value;
  RefCountedPtr<BaseBCValue> m_flux;
};

#include "NamespaceFooter.H"
#endif
