#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "EBQuadCFStencilColor.H"
#include "EBDebugOut.H"
#include <algorithm>

#include "NamespaceHeader.H"

/************************************/
EBQuadCFStencilColor::EBQuadCFStencilColor()
{
    //QuadCFStencil::QuadCFStencil();
    //EBQuadCFInterp::EBQuadCFInterp();
}
/************************************/
EBQuadCFStencilColor::~EBQuadCFStencilColor()
{
}
/************************************/
EBQuadCFStencilColor::
EBQuadCFStencilColor(const DisjointBoxLayout& a_gridsFine,
               const DisjointBoxLayout& a_gridsCoar,
               const EBISLayout&        a_ebislFine,
               const EBISLayout&        a_ebislCoar,
               const ProblemDomain&     a_domainCoar,
               const int&               a_nref,
               const int&               a_nvar,
               const LayoutData<IntVectSet>&  a_cfivs,
               const EBIndexSpace* const a_ebisPtr,
               bool a_doEBCFCrossing)
    :EBQuadCFInterp()

{
  define(a_gridsFine,   a_gridsCoar,
         a_ebislFine,   a_ebislCoar, a_domainCoar,
         a_nref, a_nvar, a_cfivs, a_ebisPtr, a_doEBCFCrossing);
}
/************************************/
void
EBQuadCFStencilColor::
define(const DisjointBoxLayout& a_gridsFine,
       const DisjointBoxLayout& a_gridsCoar,
       const EBISLayout&        a_ebislFine,
       const EBISLayout&        a_ebislCoar,
       const ProblemDomain&     a_domainCoar,
       const int&               a_nref,
       const int&               a_nvar,
       const LayoutData<IntVectSet>&  a_cfivs,
       const EBIndexSpace* const a_ebisPtr,
       bool a_doEBCFCrossing)
{
    EBQuadCFInterp::define(a_gridsFine,   a_gridsCoar,
			   a_ebislFine,   a_ebislCoar, a_domainCoar,
			   a_nref, a_nvar, a_cfivs, a_ebisPtr, a_doEBCFCrossing);
//  QuadCFStencil::define(a_gridsFine, &a_gridsCoar, dxfine, a_nref, a_nvar, m_domainFine);
}

void
EBQuadCFStencilColor::getCFVoFStencils(const VolIndex& a_vof,
				       const DataIndex& a_dit,
				       VoFStencil& a_fineStencil,
				       VoFStencil& a_coarStencil)
{
    a_fineStencil.clear();
    a_coarStencil.clear();
    const EBISBox& ebisBoxFine =  m_ebcfdata->m_ebislFine[a_dit];
    const EBISBox& ebisBoxCoar =  m_ebcfdata->m_ebislCoarsenedFine[a_dit];
    const  VolIndex vofCoar = m_ebcfdata->m_ebislFine.coarsen(a_vof, m_refRatio, a_dit);
    for (int idir = 0; idir < 3; idir++)
    {
	const IntVectSet& ebcfivsLo = m_ebcfdata->m_ebcfivsLo[idir][a_dit];
	const IntVectSet& ebcfivsHi = m_ebcfdata->m_ebcfivsHi[idir][a_dit];
	const IntVect& iv =  a_vof.gridIndex();
	if (m_doEBCFCrossing && ebisBoxFine.isIrregular(iv))
	{

	    if (ebcfivsLo.contains(iv))
	    {
		a_fineStencil +=m_fineStencilLo[idir][a_dit](a_vof, 0);
		a_coarStencil +=m_coarStencilLo[idir][a_dit](a_vof, 0);
	    }
	    
	    if (ebcfivsHi.contains(iv))
	    {
		a_fineStencil +=m_fineStencilHi[idir][a_dit](a_vof, 0);
		a_coarStencil +=m_coarStencilHi[idir][a_dit](a_vof, 0);
	    }
	}
	
	//CF for regular cells
	const IntVect& ivf = a_vof.gridIndex();
	const IntVect ivc = coarsen(ivf, m_refRatio);
	Box coarfine(ivc, ivc);
	coarfine.grow(idir, 2);
	coarfine &= m_ebcfdata->m_domainCoar;
	IntVectSet ivsCoar(coarfine);
	for (VoFIterator vofit(ivsCoar, ebisBoxCoar.getEBGraph()); vofit.ok(); ++vofit)
	{
	    a_coarStencil.add(vofit(), 1.0);
	}
	
	Box finefine(ivf, ivf);
	finefine.grow(idir, 2);
	finefine &= m_ebcfdata->m_domainFine;
	
	IntVectSet ivsFine(finefine);
	for (VoFIterator vofit(ivsFine, ebisBoxFine.getEBGraph()); vofit.ok(); ++vofit)
	{
	    a_fineStencil.add(vofit(), 1.0);
	}
	
	if (m_ebcfdata->m_edgeIVS[a_dit].contains(a_vof.gridIndex()))
	{
	    a_fineStencil += m_stencilEdges[a_dit](a_vof, 0);
	}
	if (m_ebcfdata->m_cornerIVS[a_dit].contains(a_vof.gridIndex()))
	{
	    a_fineStencil += m_stencilCorners[a_dit](a_vof, 0);
	}
    }
}
#include "NamespaceFooter.H"
