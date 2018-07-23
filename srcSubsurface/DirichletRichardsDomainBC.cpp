#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "DirichletPoissonDomainBCF_F.H"
#include "DirichletRichardsDomainBC.H"
#include "RichardsDomainBCF_F.H"
#include "EBArithF_F.H"
#include "NamespaceHeader.H"
/*****/
DirichletRichardsDomainBC::
DirichletRichardsDomainBC()
{
}
/*****/
DirichletRichardsDomainBC::
~DirichletRichardsDomainBC()
{
}

/*****/
void
DirichletRichardsDomainBC::
setValue(Real a_value)
{
   m_onlyHomogeneous = false;
   m_isFunctional = false;
   m_value = a_value;
   m_func = RefCountedPtr<BaseBCValue>();
}

/*****/
int
DirichletRichardsDomainBC::
whichBC(int                  a_idir,
        Side::LoHiSide       a_side)
{
  return 0;
}

/*****/
void
DirichletRichardsDomainBC::
setFunction(RefCountedPtr<BaseBCValue> a_func)
{
  m_value = 12345.6789;
  m_func = a_func;

  m_onlyHomogeneous = false;
  m_isFunctional = true;
}

/*****/
void
DirichletRichardsDomainBC::
getFaceFlux(BaseFab<Real>&        a_faceFlux,
            const BaseFab<Real>&  a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{

  for (int comp=0; comp<a_phi.nComp(); comp++)
    {
      const Box& box = a_faceFlux.box();

      int iside;

      if (a_side == Side::Lo)
        {
          iside = 1;
        }
      else
        {
          iside = -1;
        }

      if (a_useHomogeneous)
        {
          Real value = 0.0;

          FORT_SETDIRICHLETRICHARDSFACEFLUX(CHF_FRA1(a_faceFlux,comp),
					    CHF_CONST_FRA1(a_phi,comp),
					    CHF_CONST_REAL(value),
					    CHF_CONST_REALVECT(a_dx),
					    CHF_CONST_INT(a_idir),
					    CHF_CONST_INT(iside),
					    CHF_BOX(box));
        }
      else
        {
          if (m_isFunctional)
            {
              Real ihdx;

              ihdx = 2.0 / a_dx[a_idir];

              BoxIterator bit(box);

              for (bit.begin(); bit.ok(); ++bit)
                {
                  IntVect iv = bit();
                  IntVect ivNeigh = iv;
                  ivNeigh[a_idir] += sign(a_side);
                  const VolIndex vof      = VolIndex(iv,     0);
                  const VolIndex vofNeigh = VolIndex(ivNeigh,0);
                  const FaceIndex face = FaceIndex(vof,vofNeigh,a_idir);
                  const RealVect  point = EBArith::getFaceLocation(face,a_dx,a_probLo);
                  const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
                  Real value = m_func->value(face,a_side,a_dit,point,normal,a_time,comp);

                  Real phiVal = a_phi(iv,comp);
                  a_faceFlux(iv,comp) = iside * ihdx * (phiVal - value);
                }
            }
          else
            {
              if (m_onlyHomogeneous)
                {
                  MayDay::Error("DirichletPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
                }

              Real value = m_value;

              FORT_SETDIRICHLETRICHARDSFACEFLUX(CHF_FRA1(a_faceFlux,comp),
					       CHF_CONST_FRA1(a_phi,comp),
					       CHF_CONST_REAL(value),
					       CHF_CONST_REALVECT(a_dx),
					       CHF_CONST_INT(a_idir),
					       CHF_CONST_INT(iside),
					       CHF_BOX(box));
            }
        }
    }

  //again, following the odd convention of EBAMRPoissonOp
  //(because I am reusing its BC classes),
  //the input flux here is CELL centered and the input box
  //is the box adjacent to the domain boundary on the valid side.
  //because I am not insane (yet) I will just shift the flux's box
  //over and multiply by the appropriate coefficient
  a_faceFlux.shiftHalf(a_idir, -sign(a_side));
  const Box& faceBox = a_faceFlux.box();
  const BaseFab<Real>&   regCoef = (*m_bcoef)[a_dit][a_idir].getSingleValuedFAB();
  int  isrc = 0;
  int  idst = 0;
  int  inum = 1;
  FORT_MULTIPLYTWOFAB(CHF_FRA(a_faceFlux),
                      CHF_CONST_FRA(regCoef),
                      CHF_BOX(faceBox),
                      CHF_INT(isrc),CHF_INT(idst),CHF_INT(inum));

  //shift flux back to cell centered land
  a_faceFlux.shiftHalf(a_idir,  sign(a_side));
}

/*****/
void

DirichletRichardsDomainBC::
getFaceFlux(Real&                 a_faceFlux,
            const VolIndex&       a_vof,
            const int&            a_comp,
            const EBCellFAB&      a_phi,
            const RealVect&       a_probLo,
            const RealVect&       a_dx,
            const int&            a_idir,
            const Side::LoHiSide& a_side,
            const DataIndex&      a_dit,
            const Real&           a_time,
            const bool&           a_useHomogeneous)
{
    BaseFab<Real>&   regCoef = (*m_bcoef)[a_dit][a_idir].getSingleValuedFAB();
    regCoef.shiftHalf(a_idir, sign(a_side));

    a_faceFlux = 0.0;
    IntVect iv = a_vof.gridIndex();
    IntVect ivNeigh = iv;
    ivNeigh[a_idir] += sign(a_side);
    const VolIndex vofNeigh = VolIndex(ivNeigh,0);
    const FaceIndex face = FaceIndex(a_vof,vofNeigh,a_idir);

    getFaceGradPhi(a_faceFlux,face,a_comp,a_phi,a_probLo,a_dx,a_idir,
		   a_side,a_dit,a_time,false,RealVect::Zero,a_useHomogeneous);
    
    a_faceFlux *=regCoef(iv, 0);
    regCoef.shiftHalf(a_idir, -sign(a_side));
}

/*****/
void
DirichletRichardsDomainBC::
getFaceGradPhi(Real&                 a_faceFlux,
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
               const bool&           a_useHomogeneous)
{
  int iside = -sign(a_side);
  const Real ihdx = 2.0 / a_dx[a_idir];
  const EBISBox& ebisBox = a_phi.getEBISBox();

  Real value = -1.e99;
  RealVect point = EBArith::getFaceLocation(a_face,a_dx,a_probLo);
  if (a_useHomogeneous)
    {
      value = 0.0;
    }
  else if (m_isFunctional)
    {
      const RealVect normal = EBArith::getDomainNormal(a_idir,a_side);
      value = m_func->value(a_face,a_side,a_dit,point,normal,a_time,a_comp);
    }
  else
    {
      if (m_onlyHomogeneous)
        {
          MayDay::Error("DirichletPoissonDomainBC::getFaceFlux called with undefined inhomogeneous BC");
        }
      value = m_value;
    }

  const VolIndex& vof = a_face.getVoF(flip(a_side));
  

  a_faceFlux = iside * ihdx * (a_phi(vof,a_comp) - value);

  if (a_useAreaFrac)
    {
      MayDay::Error("DirichletPoissonDomainBC::getFaceFlux -- useAreaFrac=TRUE");
      a_faceFlux *= ebisBox.areaFrac(a_face);
    }
}

/*****/
void
DirichletRichardsDomainBC::
getFaceVel(Real&                 a_faceFlux,
           const FaceIndex&      a_face,
           const EBFluxFAB&      a_vel,
           const RealVect&       a_probLo,
           const RealVect&       a_dx,
           const int&            a_idir,
           const int&            a_icomp,
           const Real&           a_time,
           const Side::LoHiSide& a_side)
{
  CH_assert(a_idir == a_face.direction());
  Real value;
  if (m_isFunctional)
    {
      RealVect pt;
      IntVect iv = a_face.gridIndex(Side::Hi);
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (idir != a_face.direction())
            {
              Real ptval = a_dx[a_idir]*(Real(iv[idir]) + 0.5);
              pt[idir] = ptval;
            }
          else
            {
              pt[idir] = a_dx[a_idir]*(Real(iv[idir]));
            }
        }
      RealVect normal = EBArith::getDomainNormal(a_idir, a_side);

      value = m_func->value(pt, normal, a_time,a_icomp);

    }
  else
    {
      value = m_value;
    }
  a_faceFlux = value;
}
/******/
DirichletRichardsDomainBCFactory::
DirichletRichardsDomainBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = true;
  m_isFunction = false;
}

/******/
DirichletRichardsDomainBCFactory::
~DirichletRichardsDomainBCFactory()
{
}
/******/
void
DirichletRichardsDomainBCFactory::
setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}
/******/
void
DirichletRichardsDomainBCFactory::
setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
/******/
DirichletRichardsDomainBC*
DirichletRichardsDomainBCFactory::
create(const ProblemDomain& a_domain,
       const EBISLayout&    a_layout,
       const RealVect&      a_dx)
{
 DirichletRichardsDomainBC* newBC = new DirichletRichardsDomainBC();
  if (!m_onlyHomogeneous)
    {
      if (m_isFunction)
        {
          newBC->setFunction(m_flux);
        }
      else
        {
          newBC->setValue(m_value);
        }
    }
  return newBC;
}
#include "NamespaceFooter.H"
