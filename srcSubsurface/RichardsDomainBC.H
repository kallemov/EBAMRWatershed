#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif


#ifndef _RICHARDSDOMAINBC_H_
#define _RICHARDSDOMAINBC_H_

#include "BaseDomainBC.H"
#include "ParserBCValue.H"
#include "ConductivityBaseDomainBC.H"

//this copy-paste part is taken from ConductivityIBC.H file
#include "NamespaceHeader.H"

class RichardsDomainBC: public ConductivityBaseDomainBC
{
public:

  ///
  /**
   */
  virtual ~RichardsDomainBC()
  {;}

  virtual void fillPhiGhost(FArrayBox&     a_state,
			    const Box&     a_valid,
			    const Box&     a_domain,
			    Real           a_dx,
			    bool           a_homogeneous) {;}
    

  virtual void getFluxStencil(      VoFStencil&      a_stencil,
                              const VolIndex&        a_vof,
                              const int&             a_comp,
                              const RealVect&        a_vectDx,
                              const int&             a_idir,
                              const Side::LoHiSide&  a_side,
                              const EBISBox&         a_ebisBox);

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_vectDx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  ///
  /**
     Elliptic solver flux.
  */
  virtual void getFaceFlux(Real&                 a_faceFlux,
                           const VolIndex&       a_vof,
                           const int&            a_comp,
                           const EBCellFAB&      a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_vectDx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous);

  virtual void getInhomFaceFlux(Real&                 a_faceFlux,
                                const VolIndex&       a_vof,
                                const int&            a_comp,
                                const EBCellFAB&      a_phi,
                                const RealVect&       a_probLo,
                                const RealVect&       a_vectDx,
                                const int&            a_idir,
                                const Side::LoHiSide& a_side,
                                const DataIndex&      a_dit,
                                const Real&           a_time);

  virtual void getFaceGradPhi(Real&                 a_faceFlux,
                              const FaceIndex&      a_face,
                              const int&            a_comp,
                              const EBCellFAB&      a_phi,
                              const RealVect&       a_probLo,
                              const RealVect&       a_vectDx,
                              const int&            a_idir,
                              const Side::LoHiSide& a_side,
                              const DataIndex&      a_dit,
                              const Real&           a_time,
                              const bool&           a_useAreaFrac,
                              const RealVect&       a_centroid,
                              const bool&           a_useHomogeneous);


  RichardsDomainBC(const IntVect&        a_domainLoBCType,
		   const IntVect&        a_domainHiBCType,
		   const Vector<RefCountedPtr<ParserBCValue> >& a_domainLoBCValueFunc,
		   const Vector<RefCountedPtr<ParserBCValue> >& a_domainHiBCValueFunc,
		   const ProblemDomain&  a_domain,
		   const RealVect&       a_vectDx)
  {
    m_domainLoBCType  = a_domainLoBCType;
    m_domainHiBCType  = a_domainHiBCType;
    m_domainLoBCValueFunc = a_domainLoBCValueFunc;
    m_domainHiBCValueFunc = a_domainHiBCValueFunc;
    m_domain          = a_domain;
    m_vectDx              = a_vectDx;
  }

  ///
  /**
     A query of whether a_jvof is Dirichlet Domain boundary to a_ivof. a_jvof is the ghost vof; a_ivof is the valid computational vof.
     For now this is mainly used in form_matrix, which needs to know the type of BC at the boundary.
     Default implementation is MayDay::Error();
   */
  virtual bool isDirichletDom(const VolIndex&   a_ivof,
                              const VolIndex&   a_jvof,
                              const EBCellFAB&  a_phi) const ;


  //get face value on domain BC
  virtual Real getBCValue(const FaceIndex&      a_faceidx,
			  const EBCellFAB&      a_phi,
			  const EBFaceFAB&      a_bcoeff,
			  const Side::LoHiSide& a_side,
			  const int&            a_dir,
			  const RealVect&       a_vectDx,
			  const Real            a_time,
			  const int             a_comp);

private:
  //weak construction forbidden
    // RichardsDomainBC(){;};

  IntVect        m_domainLoBCType;
  IntVect        m_domainHiBCType;
  Vector<RefCountedPtr<ParserBCValue> >  m_domainLoBCValueFunc;
  Vector<RefCountedPtr<ParserBCValue> >  m_domainHiBCValueFunc;
  ProblemDomain  m_domain;
  RealVect       m_vectDx;
};

///
/**
 **/
class RichardsDomainBCFactory: public BaseDomainBCFactory
{
public:

  ///
  RichardsDomainBCFactory(const IntVect&        a_domainLoBCType,
			  const IntVect&        a_domainHiBCType,
			  const Vector<RefCountedPtr<ParserBCValue> >& a_domainLoBCValueFunc,
			  const Vector<RefCountedPtr<ParserBCValue> >& a_domainHiBCValueFunc)
  {
    m_domainLoBCType  = a_domainLoBCType;
    m_domainHiBCType  = a_domainHiBCType;
    m_domainLoBCValueFunc = a_domainLoBCValueFunc;
    m_domainHiBCValueFunc = a_domainHiBCValueFunc;
  }

  ///
  /**
   */
  virtual ~RichardsDomainBCFactory()
  {
  }

  ///
  /**
   */
  virtual RichardsDomainBC* create(const ProblemDomain& a_domain,
				   const EBISLayout&    a_layout,
				   const RealVect&      a_vectDx)
  {
    RichardsDomainBC* newBC
      = new RichardsDomainBC(m_domainLoBCType,
			     m_domainHiBCType,
			     m_domainLoBCValueFunc,
			     m_domainHiBCValueFunc,
			     a_domain,
			     a_vectDx);

    return newBC;
  }

private:
  //weak construction
  RichardsDomainBCFactory()
  {
    MayDay::Error("RichardsDomainBC:: weak constructor is not supported");
  }

  IntVect   m_domainLoBCType;
  IntVect   m_domainHiBCType;
  Vector<RefCountedPtr<ParserBCValue> >  m_domainLoBCValueFunc;
  Vector<RefCountedPtr<ParserBCValue> >  m_domainHiBCValueFunc;
};

#include "NamespaceFooter.H"

#endif
