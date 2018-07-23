#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "VoFIterator.H"
#include "EBLevelDataOps.H"
#include "WatershedIBC.H"
#include "SpreadingCopier.H"
#include "ReductionCopier.H"
#include "ReductionOps.H"
#include "LevelDataOps.H"
#include "EBAMRDataOps.H"
#include "ParserFunc.H"
#include "EBSUBSurfaceOps.H"

#include "UsingNamespace.H"

WatershedIBC::WatershedIBC(const bool  a_includeSurfaceSolver,
			   const IntVect&  a_scalarDomainLoBCType,
			   const IntVect&  a_scalarDomainHiBCType,
			   const Vector<RefCountedPtr<ParserBCValue> >& a_scalarDomainLoBCValueFunc,
			   const Vector<RefCountedPtr<ParserBCValue> >& a_scalarDomainHiBCValueFunc,
			   const RefCountedPtr<ParserBCValue>&   a_scalarEBBCValueFunc)
{
  m_finestLevel = 0;

  m_scalarDomainLoBCType  = a_scalarDomainLoBCType;
  m_scalarDomainHiBCType  = a_scalarDomainHiBCType;
  m_scalarDomainLoBCValueFunc = a_scalarDomainLoBCValueFunc;
  m_scalarDomainHiBCValueFunc = a_scalarDomainHiBCValueFunc;
  m_scalarEBBCValueFunc = a_scalarEBBCValueFunc;

  //include surface solver
  m_includeSurfaceSolver = a_includeSurfaceSolver; 
  // if (m_includeSurfaceSolver)
  // {
  //     //allocate data on EB exchange between surface and subsurface solver
  //     const int nlevels = m_watershedPtr->m_params.m_maxLevel + 1;
  //      m_TempFlux.resize(nlevels);
  //     for (int ilev = 0; ilev<nlevels; ilev++)
  //     {
  // 	  m_TempFlux[ilev]   = RefCountedPtr<LevelData<EBFluxFAB> >(new LevelData<EBFluxFAB>());
  //     }
  // }
}

WatershedIBC::~WatershedIBC()
{

}


RefCountedPtr<RichardsDomainBCFactory>
WatershedIBC::getScalarBC() const
{
  //allocate BC holder
  RefCountedPtr<RichardsDomainBCFactory> richardsDomainBCFac = 
    RefCountedPtr<RichardsDomainBCFactory>(new RichardsDomainBCFactory(m_scalarDomainLoBCType,
								       m_scalarDomainHiBCType,
								       m_scalarDomainLoBCValueFunc,
								       m_scalarDomainHiBCValueFunc));

  return richardsDomainBCFac;
}


RefCountedPtr<RichardsEBBCFactory>
WatershedIBC::getScalarEBBC() const
{
  //allocate EBBC holder
  //this ebbc is neumann type always based on data
  //and must be set if surface solver is not used 
  RichardsEBBCFactory*    ebbcRichards =  new RichardsEBBCFactory();
  if (!m_includeSurfaceSolver) 
  {
      ebbcRichards->setOverlandFlowSolver(false);
      ebbcRichards->setFunction(RefCountedPtr<BaseBCValue>(m_scalarEBBCValueFunc));
      // ebbcRichards->setValue(0.0);
  }
  else 
  {
      ebbcRichards->setOverlandFlowSolver(true);
  }
  
  RefCountedPtr<RichardsEBBCFactory>  ebbc = RefCountedPtr<RichardsEBBCFactory>(ebbcRichards);

  return ebbc;
}


