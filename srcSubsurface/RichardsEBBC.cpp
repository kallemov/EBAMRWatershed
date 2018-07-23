#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "BoxIterator.H"

#include "RichardsEBBC.H"
#include "EBStencil.H"
#include "EBSUBSurfaceOps.H"
#include "NamespaceHeader.H"

/*****************/
RichardsEBBC::RichardsEBBC(const ProblemDomain& a_domain,
			   const EBISLayout&    a_layout,
			   const RealVect&      a_vectDx)
  :m_bc(a_domain, a_layout, a_vectDx)
{
  m_dataBased = false;
  m_OverlandFlowSolver = true;
}
/*****************/
RichardsEBBC::~RichardsEBBC()
{
}
/*****************/
void RichardsEBBC::setValue(Real a_value)
{
  m_bc.setValue(a_value);
}
/*****************/
void RichardsEBBC::setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_bc.setFunction(a_flux);
}
/*****************/
void RichardsEBBC::getEBFlux(Real&                         a_flux,
			     const VolIndex&               a_vof,
			     const LevelData<EBCellFAB>&   a_phi,
			     const LayoutData<IntVectSet>& a_cfivs,
			     const DataIndex&              a_dit,
			     const RealVect&               a_probLo,
			     const RealVect&               a_vectDx,
			     const bool&                   a_useHomogeneous,
			     const Real&                   a_time,
			     const pair<int,Real>*         a_cacheHint )
{
    if (m_dataBased)
    {
	a_flux = (*m_data)[a_dit](a_vof, 0);
    }
    else
    {
	m_bc.getEBFlux(a_flux,
		       a_vof,
		       a_phi,
		       a_cfivs,
		       a_dit,
		       a_probLo,
		       a_vectDx,
		       a_useHomogeneous,
		       a_time,
		       a_cacheHint );
    }

  // Real bcoef = (*m_bcoe)[a_dit](a_vof,0);
  // a_flux *= bcoef;
}
/*****************/
void RichardsEBBC::applyEBFlux(EBCellFAB&                    a_lphi,
			       const EBCellFAB&              a_phi,
			       VoFIterator&                  a_vofit,
			       const LayoutData<IntVectSet>& a_cfivs,
			       const DataIndex&              a_dit,
			       const RealVect&               a_probLo,
			       const RealVect&               a_vectDx,
			       const Real&                   a_factor,
			       const bool&                   a_useHomogeneous,
			       const Real&                   a_time)
{
  CH_TIME("RichardsEBBC::applyEBFlux");
  CH_assert(a_lphi.nComp() == 1 );
  CH_assert(a_phi.nComp()  == 1);

  Real flux = 0.0;
  const Real Tol = 1.e-15;

  const EBISBox&   ebisBox = a_phi.getEBISBox();
  for (a_vofit.reset(); a_vofit.ok(); ++a_vofit)
    {
      const VolIndex& vof = a_vofit();
      const RealVect&   normal = ebisBox.normal(vof);

      if (m_dataBased)
        {
//          if ((*m_data)[a_dit].getIVS().contains(vof.gridIndex()))
//             {
//               flux = (*m_data)[a_dit](vof, 0);
//             }
//          else
//            {
//              flux = 0.;
//            }
               flux = (*m_data)[a_dit](vof, 0);
        }
      else if (m_bc.m_isFunction)
        {
          const RealVect& centroid = ebisBox.bndryCentroid(vof);
	  // const RealVect&   normal = ebisBox.normal(vof);

          Real value = m_bc.m_flux->value(vof,centroid,normal,a_vectDx,a_probLo,a_dit,a_time,0);
          flux = -value;
        }
      else
        {
          if (m_bc.m_onlyHomogeneous)
            {
              MayDay::Error("RichardsEBBC::getFaceFlux called with undefined inhomogeneous BC");
            }

          flux = m_bc.m_value;
        }

      const Real areaFrac = ebisBox.bndryArea(vof);
      flux *= areaFrac;
      //      pout()<< areaFrac<<endl;
      // if(m_OverlandFlowSolver)
      // 	{
      // 	  flux *= m_beta;//*bcoef;
      // 	}
      // else
      // 	{
      // 	  //Real bcoef = (*m_bcoe)[a_dit](vof,0);
      // 	  flux *= m_beta;//*bcoef;
      // 	}

      flux *= m_beta;
      //const RealVect anormal = getAnisotropicNormal(normal, a_vectDx);
      Real anisotr_corr=0.0;
      for (int idir = 0; idir < SpaceDim; idir++)
      {
	  anisotr_corr += normal[idir]*normal[idir] 
	                 /(a_vectDx[idir]*a_vectDx[idir]);
      }
      if (areaFrac>Tol)
      {
	  a_lphi(vof,0) += flux * a_factor*sqrt(anisotr_corr);//*fabs(anormal[2]);
      }
    }
}


/*****************/
RichardsEBBCFactory::RichardsEBBCFactory()
{
  m_value = 12345.6789;
  m_flux = RefCountedPtr<BaseBCValue>();
  m_onlyHomogeneous = false; //don't need onlyhomogeneous type here
  m_isFunction = false;
  m_dataBased = false; //will be set by setData function
  m_OverlandFlowSolver = true;
}
/*****************/
RichardsEBBCFactory::~RichardsEBBCFactory()
{
}
/*****************/
void RichardsEBBCFactory::setValue(Real a_value)
{
  m_value = a_value;
  m_flux = RefCountedPtr<BaseBCValue>();

  m_onlyHomogeneous = false;
  m_isFunction = false;
}
/*****************/
void RichardsEBBCFactory::setFunction(RefCountedPtr<BaseBCValue> a_flux)
{
  m_value = 12345.6789;
  m_flux = a_flux;

  m_onlyHomogeneous = false;
  m_isFunction = true;
}
/*****************/
RichardsEBBC* RichardsEBBCFactory::create(const ProblemDomain& a_domain,
					  const EBISLayout&    a_layout,
					  const RealVect&      a_vectDx,
					  const IntVect*       a_ghostCellsPhi /*=0*/,
					  const IntVect*       a_ghostCellsRhs /*=0*/)
{
  RichardsEBBC* fresh = new RichardsEBBC(a_domain,a_layout,a_vectDx);

  if (!m_onlyHomogeneous)
    {
      if (m_OverlandFlowSolver)
	{
	  fresh->setOverlandFlowSolver(m_OverlandFlowSolver); 
          fresh->setData(m_data);
	  //overland ebbc based only on data from surface solver
	} 
      else if (m_dataBased)
        {
          fresh->setData(m_data);
        }
      else if (!m_isFunction)
        {
          fresh->setValue(m_value);
        }
      else
        {
          fresh->setFunction(m_flux);
        }
    }

  return fresh;
}
#include "NamespaceFooter.H"
