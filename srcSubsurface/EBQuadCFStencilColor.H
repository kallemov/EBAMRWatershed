#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// dtgraves fri, sept 14, 2001
// mfbarad made it eb aware in June, 2005
// and on and on...

#ifndef _EBQUADCFSTENCIL_H_
#define _EBQUADCFSTENCIL_H_

#include "EBQuadCFInterp.H"
#include "QuadCFInterp.H"
#include "NamespaceHeader.H"

///
/**
   Fills one layer of fine level ghost cells
   by piecewise linear interpolation
   from the coarse level.
 */
class EBQuadCFStencilColor: public EBQuadCFInterp
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  EBQuadCFStencilColor();

  ///
  virtual ~EBQuadCFStencilColor();

  EBQuadCFStencilColor(const DisjointBoxLayout& a_dblFine,
                 const DisjointBoxLayout& a_dblCoar,
                 const EBISLayout&        a_ebislFine,
                 const EBISLayout&        a_ebislCoar,
                 const ProblemDomain&     a_domainCoar,
                 const int&               a_nref,
                 const int&               a_nvar,
                 const LayoutData<IntVectSet>&  a_cfivs,
                 const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance(),
                 bool a_doEBCFCrossing = true);

  void define(const DisjointBoxLayout& a_dblFine,
              const DisjointBoxLayout& a_dblCoar,
              const EBISLayout&        a_ebislFine,
              const EBISLayout&        a_ebislCoar,
              const ProblemDomain&     a_domainCoar,
              const int&               a_nref,
              const int&               a_nvar,
              const LayoutData<IntVectSet>&  a_cfivs,
              const EBIndexSpace* const a_ebisPtr = Chombo_EBIS::instance(),
              bool a_doEBCFCrossing = true);
    
    void getCFVoFStencils(const VolIndex& a_vof,
			  const DataIndex& a_dit,
			  VoFStencil& a_fineStencil,
			  VoFStencil& a_coarStencil);



private:

  //disallowed for all the usual reasons
  EBQuadCFStencilColor(const EBQuadCFStencilColor& ebcin)
  {
    MayDay::Error("ebqcfi 2 invalid operator");
  }
  void operator=(const EBQuadCFStencilColor& fabin)
  {
    MayDay::Error("ebqcfi 3 invalid operator");
  }
};

#include "NamespaceFooter.H"
#endif
