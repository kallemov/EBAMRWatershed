#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "RichardsDomainBC.H"
#include "VoFIterator.H"
#include "NeumannRichardsDomainBC.H"
#include "DirichletRichardsDomainBC.H"
#include "EBLevelDataOps.H"
#include "ParmParse.H"

#include "NamespaceHeader.H"

void RichardsDomainBC::getFluxStencil(VoFStencil&      a_stencil,
				      const VolIndex&        a_vof,
				      const int&             a_comp,
				      const RealVect&        a_vectDx,
				      const int&             a_idir,
				      const Side::LoHiSide&  a_side,
				      const EBISBox&         a_ebisBox)
{
  if(a_side == Side::Lo)
    {
      if(m_domainLoBCType[a_idir] == 1)
        {
          DirichletRichardsDomainBC diriBC;
          diriBC.getFluxStencil(a_stencil, a_vof, a_comp, a_vectDx, a_idir, a_side, a_ebisBox);
        }
      else
        {
          NeumannRichardsDomainBC neumannBC;
          neumannBC.getFluxStencil(a_stencil, a_vof, a_comp, a_vectDx, a_idir, a_side, a_ebisBox);
        }
    }
  else if (a_side == Side::Hi)
    {
      if(m_domainHiBCType[a_idir] == 1)
        {
          DirichletRichardsDomainBC diriBC;
          diriBC.getFluxStencil(a_stencil, a_vof, a_comp, a_vectDx, a_idir, a_side, a_ebisBox);
        }
      else
        {
          NeumannRichardsDomainBC neumannBC;
          neumannBC.getFluxStencil(a_stencil, a_vof, a_comp, a_vectDx, a_idir, a_side, a_ebisBox);
        }
    }
  else
    {
      MayDay::Error("RichardsDomainBC::getFluxStencil -- bad side");
    }
}

void RichardsDomainBC::
getFaceFlux(BaseFab<Real>&        a_faceFlux,
            const BaseFab<Real>&  a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_vectDx,
            const int&            a_dir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  if(a_side == Side::Lo)
    {
      if (m_domainLoBCType[a_dir] == 0)
        {
          NeumannRichardsDomainBC neumannBC;
          neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
          neumannBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainLoBCValueFunc[a_dir]));
          neumannBC.getFaceFlux(a_faceFlux,
                                a_phi,
                                a_probLo,
                                a_vectDx,
                                a_dir,
                                a_side,
                                a_dit,
                                a_time,
                                a_useHomogeneous);
        }
      else if (m_domainLoBCType[a_dir] == 1)
        {
          DirichletRichardsDomainBC diriBC;
          diriBC.setCoef(m_eblg, m_beta, m_bcoef);
	  diriBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainLoBCValueFunc[a_dir]));
          diriBC.getFaceFlux(a_faceFlux,
			     a_phi,
			     a_probLo,
			     a_vectDx,
			     a_dir,
			     a_side,
			     a_dit,
			     a_time,
			     a_useHomogeneous);
        }
      else
        {
          MayDay::Error("ConductivityIBC: bad BC type");
        }
    }
  else if(a_side == Side::Hi)
    {
      if (m_domainHiBCType[a_dir] == 0)
        {
          NeumannRichardsDomainBC neumannBC;
          neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
          neumannBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainHiBCValueFunc[a_dir]));
          neumannBC.getFaceFlux(a_faceFlux,
                                a_phi,
                                a_probLo,
                                a_vectDx,
                                a_dir,
                                a_side,
                                a_dit,
                                a_time,
                                a_useHomogeneous);
        }
      else if (m_domainHiBCType[a_dir] == 1)
        {
          DirichletRichardsDomainBC diriBC;
          diriBC.setCoef(m_eblg, m_beta, m_bcoef);
          diriBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainHiBCValueFunc[a_dir]));
          diriBC.getFaceFlux(a_faceFlux,
			     a_phi,
			     a_probLo,
			     a_vectDx,
			     a_dir,
			     a_side,
			     a_dit,
			     a_time,
			     a_useHomogeneous);
        }
      else
        {
          MayDay::Error("ConductivityIBC: bad BC type");
        }
    }
}
/***/
void RichardsDomainBC::
getFaceFlux(Real&                 a_faceFlux,
            const VolIndex&       a_vof,
            const int&            a_comp,
            const EBCellFAB&      a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_vectDx,
            const int&            a_dir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
  if(a_side == Side::Lo)
    {
      if (m_domainLoBCType[a_dir] == 0)
        {
          NeumannRichardsDomainBC neumannBC;
          neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
          neumannBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainLoBCValueFunc[a_dir]));
          neumannBC.getFaceFlux(a_faceFlux,
                                a_vof,
                                a_comp,
                                a_phi,
                                a_probLo,
                                a_vectDx,
                                a_dir,
                                a_side,
                                a_dit,
                                a_time,
                                a_useHomogeneous);
        }
      else if (m_domainLoBCType[a_dir] == 1)
        {
          DirichletRichardsDomainBC diriBC;
          diriBC.setCoef(m_eblg, m_beta, m_bcoef);
	  diriBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainLoBCValueFunc[a_dir]));
          diriBC.getFaceFlux(a_faceFlux,
			     a_vof,
			     a_comp,
			     a_phi,
			     a_probLo,
			     a_vectDx,
			     a_dir,
			     a_side,
			     a_dit,
			     a_time,
			     a_useHomogeneous);
        }
      else
        {
          MayDay::Error("ConductivityIBC: bad BC type");
        }
    }
  else if(a_side == Side::Hi)
    {
      if (m_domainHiBCType[a_dir] == 0)
        {
          NeumannRichardsDomainBC neumannBC;
          neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
          neumannBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainHiBCValueFunc[a_dir]));
          neumannBC.getFaceFlux(a_faceFlux,
                                a_vof,
                                a_comp,
                                a_phi,
                                a_probLo,
                                a_vectDx,
                                a_dir,
                                a_side,
                                a_dit,
                                a_time,
                                a_useHomogeneous);
        }
      else if (m_domainHiBCType[a_dir] == 1)
        {
          DirichletRichardsDomainBC diriBC;
          diriBC.setCoef(m_eblg, m_beta, m_bcoef);
          diriBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainHiBCValueFunc[a_dir]));
          diriBC.getFaceFlux(a_faceFlux,
			     a_vof,
			     a_comp,
			     a_phi,
			     a_probLo,
			     a_vectDx,
			     a_dir,
			     a_side,
			     a_dit,
			     a_time,
			     a_useHomogeneous);
        }
      else
        {
          MayDay::Error("ConductivityIBC: bad BC type");
        }
    }
}

/***/
void RichardsDomainBC::
getInhomFaceFlux(Real&                 a_faceFlux,
                 const VolIndex&       a_vof,
                 const int&            a_comp,
                 const EBCellFAB&      a_phi,
                 const RealVect&       a_probLo,
                 const RealVect&       a_vectDx,
                 const int&            a_dir,
                 const Side::LoHiSide& a_side,
                 const DataIndex&      a_dit,
                 const Real&           a_time)
{
  if(a_side == Side::Lo)
    {
      if (m_domainLoBCType[a_dir] == 0)
        {
          NeumannRichardsDomainBC neumannBC;
          neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
          neumannBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainLoBCValueFunc[a_dir]));
          neumannBC.getInhomFaceFlux(a_faceFlux,
                                     a_vof,
                                     a_comp,
                                     a_phi,
                                     a_probLo,
                                     a_vectDx,
                                     a_dir,
                                     a_side,
                                     a_dit,
                                     a_time);
        }
      else if (m_domainLoBCType[a_dir] == 1)
        {
          DirichletRichardsDomainBC diriBC;
          diriBC.setCoef(m_eblg, m_beta, m_bcoef);
          diriBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainLoBCValueFunc[a_dir]));

          diriBC.getInhomFaceFlux(a_faceFlux,
				  a_vof,
				  a_comp,
				  a_phi,
				  a_probLo,
				  a_vectDx,
				  a_dir,
				  a_side,
				  a_dit,
				  a_time);
        }
      else
        {
          MayDay::Error("ConductivityIBC: bad BC type");
        }
    }
  else if(a_side == Side::Hi)
    {
      if (m_domainHiBCType[a_dir] == 0)
        {
          NeumannRichardsDomainBC neumannBC;
          neumannBC.setCoef(m_eblg, m_beta, m_bcoef);
          neumannBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainHiBCValueFunc[a_dir]));
          neumannBC.getInhomFaceFlux(a_faceFlux,
                                     a_vof,
                                     a_comp,
                                     a_phi,
                                     a_probLo,
                                     a_vectDx,
                                     a_dir,
                                     a_side,
                                     a_dit,
                                     a_time);
        }
      else if (m_domainHiBCType[a_dir] == 1)
        {
          DirichletRichardsDomainBC diriBC;
          diriBC.setCoef(m_eblg, m_beta, m_bcoef);
          diriBC.setFunction(RefCountedPtr<BaseBCValue>(m_domainHiBCValueFunc[a_dir]));
          diriBC.getInhomFaceFlux(a_faceFlux,
				  a_vof,
				  a_comp,
				  a_phi,
				  a_probLo,
				  a_vectDx,
				  a_dir,
				  a_side,
				  a_dit,
				  a_time);
        }
      else
        {
          MayDay::Error("ConductivityIBC: bad BC type");
        }
    }
}

/***/
void
RichardsDomainBC::
getFaceGradPhi(Real&                 a_faceFlux,
               const FaceIndex&      a_face,
               const int&            a_comp,
               const EBCellFAB&      a_phi,
               const RealVect&       a_probLo,
               const RealVect&       a_vectDx,
               const int&            a_dir,
               const Side::LoHiSide& a_side,
               const DataIndex&      a_dit,
               const Real&           a_time,
               const bool&           a_useAreaFrac,
               const RealVect&       a_centroid,
               const bool&           a_useHomogeneous)
{
  MayDay::Error("RichardsDomainBC::getFaceGradPhi: not needed");
}

/***/
bool
RichardsDomainBC::
isDirichletDom(const VolIndex&   a_ivof,
               const VolIndex&   a_jvof,
               const EBCellFAB&  a_phi) const
{
  IntVect diff = a_jvof.gridIndex() - a_ivof.gridIndex();
  for(int idir = 0; idir < SpaceDim; idir++)
    {
      if(diff[idir] == 1)
        {
          if(m_domainHiBCType[idir] == 1)
            {
              return true;
            }
          else
            {
              return false;
            }
        }
      else if(diff[idir] == -1)
        {
          if(m_domainLoBCType[idir] == 1)
            {
              return true;
            }
          else
            {
              return false;
            }
        }
    }
  MayDay::Error("ConductivityIBC::isDirichletDom -- not a domain boundary");
  return true;
}


Real RichardsDomainBC::
getBCValue(const FaceIndex&      a_faceidx,
	   const EBCellFAB&      a_phi,
	   const EBFaceFAB&      a_bcoeff,
	   const Side::LoHiSide& a_side,
	   const int&            a_dir,
	   const RealVect&       a_vectDx,
	   const Real            a_time,
	   const int             a_comp)
{
    const Real Tol = 1e-15;

    RealVect point = EBArith::getFaceLocation(a_faceidx, a_vectDx, RealVect::Zero);
  const RealVect normal = EBArith::getDomainNormal(a_dir,a_side);
  VolIndex vof0;
  Real value=0;
  if(a_side == Side::Lo)
  {
      if (m_domainLoBCType[a_dir] == 0)
      {
	  const VolIndex& vof = a_faceidx.getVoF(Side::flip(a_side));
	  Real flux = -m_domainLoBCValueFunc[a_dir]->parserValue(point,normal, a_time,a_comp);
	  if (fabs(a_bcoeff(a_faceidx,0))>Tol)
	  {
	      flux /=a_bcoeff(a_faceidx,0);
	  }
	  else
	  {
	      flux = 0.0;
	  }

	  value = a_phi(vof, a_comp) - 0.5*flux*a_vectDx[a_dir]; 
	  if (a_dir==2) value += 0.5*a_vectDx[a_dir];

	  vof0=vof;
      }
      else if (m_domainLoBCType[a_dir] == 1)
      {
	  value = m_domainLoBCValueFunc[a_dir]->parserValue(point,normal, a_time,a_comp);
      }
      else
      {
          MayDay::Error("ConductivityIBC: bad BC type");
      }
  }
  else if(a_side == Side::Hi)
  {
      if (m_domainHiBCType[a_dir] == 0)
      {
	  const VolIndex& vof = a_faceidx.getVoF(Side::flip(a_side));
	  Real flux = -m_domainHiBCValueFunc[a_dir]->parserValue(point,normal, a_time,a_comp);

	  if (fabs(a_bcoeff(a_faceidx,0))>Tol)
	  {
	      flux /=a_bcoeff(a_faceidx,0);
	  }
	  else
	  {
	      flux = 0.0;
	  }

	  value = a_phi(vof, a_comp) + 0.5*flux*a_vectDx[a_dir]; 
	  if (a_dir==2) value -= 0.5*a_vectDx[a_dir];

	  vof0=vof;
      }
      else if (m_domainHiBCType[a_dir] == 1)
      {
	  value = m_domainHiBCValueFunc[a_dir]->parserValue(point,normal, a_time,a_comp);
      }
      else
      {
          MayDay::Error("ConductivityIBC: bad BC type");
      }
  }
  // pout()<<a_faceidx<<"\t point="<<point<<"\t vof="<<vof0<<"\t value="<<value<<"\t"<< a_phi(vof0, a_comp)<<endl;
  return value;
}
#include "NamespaceFooter.H"
