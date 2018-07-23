#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBSUBSURFACEOps_H_
#define _EBSUBSURFACEOps_H_
#include "EBLevelGrid.H"

#include "NamespaceHeader.H"

extern void getSaturation(EBCellFAB&        a_F,
			  const EBCellFAB& a_psi,
			  const Real a_Ssat,
			  const Real a_Sres,
			  const Real a_alpha,
			  const Real a_nsoil,
			  const ProblemDomain& a_domain,
			  const Box&       a_box);

//computes and adding saturation terms to Newton-Krylov vector F(x)
extern void addSaturationTerms(EBCellFAB&        a_F,
			       const EBCellFAB& a_psi,
			       const EBCellFAB& a_psi_old,
			       const EBCellFAB& a_Ss,
			       const EBCellFAB& a_Phi,
			       //const EBCellFAB& a_Rho,
			       const Real a_Rho,
			       const Real a_Ssat,
			       const Real a_Sres,
			       const Real a_alpha,
			       const Real a_nsoil,
			       const ProblemDomain& a_domain,
			       const Box&       a_box);

extern void irregSaturationTerms(BaseIVFAB<Real>&   a_F,
				 const BaseIVFAB<Real>& a_psi,
				 const EBCellFAB& a_psi_old,
				 const EBCellFAB& a_Ss,
				 const EBCellFAB& a_Phi,
				 //const EBCellFAB& a_Rho,
				 const Real rho,
				 const Real a_Ssat,
				 const Real a_Sres,
				 const Real a_alpha,
				 const Real a_nsoil,
				 const ProblemDomain& a_domain,
				 const Box&       a_box);

extern void ccpHarmonicAverageToFaces(EBFaceFAB &      a_face,
			       const EBCellFAB &       a_cell,
			       const EBGraph &         a_ebGraph,
			       const Box &             a_grid,
			       const int &             a_idir,
			       const ProblemDomain &   a_domain,
			       const int &             a_comp);

extern void ccpHarmonicAverageToFaces(LevelData<EBFluxFAB> &  a_Flux,
				      const LevelData<EBCellFAB> &   a_CCcell,
				      const DisjointBoxLayout &      a_grids,
				      const EBISLayout &             a_ebisl,
				      const ProblemDomain &          a_domain,
				      const int                      a_comp,
      			              const int                      a_lev);


extern void ccpHarmonicAverageToFaces(LevelData<EBFluxFAB> &         a_Flux,
				      const LevelData<EBCellFAB> &   a_CCcell,
				      const DisjointBoxLayout &      a_grids,
				      const EBISLayout &             a_ebisl,
				      const ProblemDomain &          a_domain,
    			              const int                      a_lev);


extern void scaleLevelFluxByKr(LevelData<EBFluxFAB> &         a_Flux,
			       const LevelData<EBCellFAB> &   a_psi,
			       const LevelData<EBFluxFAB> &   a_darcy,
			       const Real a_alpha,
			       const Real a_nsoil,
			       const Real a_lambda,
			       const DisjointBoxLayout &      a_grids,
			       const EBISLayout &             a_ebisl,
			       const ProblemDomain &          a_domain,
			       const int& a_method,
    			       const int                      a_lev);


extern void scaleFaceByKr(EBFaceFAB &             a_face,
			  const EBCellFAB &       a_psi,
			  const EBFaceFAB &       a_darcy,
			  const Real a_alpha,
			  const Real a_nsoil,
			  const Real a_lambda,
			  const EBGraph &         a_ebGraph,
			  const Box &             a_grid,
			  const int &             a_idir,
			  const ProblemDomain &   a_domain,
			  const int &             a_method);

extern Real getRelativePermeability(const Real a_psi,
				    const Real a_alpha,
				    const Real a_nsoil,
				    const Real a_lambda);

extern Real getDerivRelPermeability(const Real a_psi,
				    const Real a_alpha,
				    const Real a_nsoil,
				    const Real a_lambda);

extern Real getDerivSaturation(const Real a_psi,
			       const Real a_Ssat,
			       const Real a_Sres,
			       const Real a_alpha,
			       const Real a_nsoil);


extern void
ccpAverageFaceToCells(LevelData<EBCellFAB>&         a_cellData,
                      const LevelData<EBFluxFAB>&   a_macData,
                      const DisjointBoxLayout&      a_grids,
                      const EBISLayout&             a_ebisl,
                      const ProblemDomain&          a_domain,
                      const RealVect&               a_dx);

extern void
ccpLevelccpAverageFaceToCells(LevelData<EBCellFAB>&       a_cellData,
			      const LevelData<EBFluxFAB>& a_fluxData,
			      const DisjointBoxLayout &      a_grids,
			      const EBISLayout &             a_ebisl,
			      const ProblemDomain &          a_domain,
			      const RealVect &        a_dx,
    			       const int                      a_lev);


extern RealVect 
getAnisotropicNormal(const RealVect& a_normal,
		     const RealVect& a_vectDx);


extern void 
copyBIVF2EBCellFAB(const Vector<LevelData<BaseIVFAB<Real> >*>& a_bivf, 
		   Vector<LevelData<EBCellFAB>*>& a_ebcellf,
		   const Vector<EBLevelGrid>&        a_ebisl);

extern void 
copyBIVF2EBCellFAB(const Vector<LevelData<BaseIVFAB<Real> >*>& a_bivf, 
		   Vector<LevelData<EBCellFAB>*>& a_ebcellf,
		   const Vector<const EBLevelGrid*>&        a_ebisl);


extern void 
copyEBCellFAB2BIVF(const Vector<LevelData<EBCellFAB>*>& a_ebcellf,
		   Vector<LevelData<BaseIVFAB<Real> >*>& a_bivf,
		   const Vector<EBLevelGrid>&        a_ebisl);

extern Real 
getAnisotropicCorrection(const VolIndex& a_vof ,
			 const RealVect& a_vectDx,
			 const EBISBox& a_ebisBox);


extern bool
findClosestIrrVoF(VolIndex& a_vof,
		  const Vector<VolIndex>& a_VoFs,
		  const RealVect& a_center,
		  const RealVect& a_dx,
		  const EBISBox& a_ebisBox);

extern bool
findClosestIrrVoF(VolIndex& a_vof,
		  const Vector<VolIndex>& a_VoFs,
		  const RealVect& a_center,
		  const int a_dir,
		  const Side::LoHiSide& a_side,
		  const RealVect& a_dx,
		  const EBISBox& a_ebisBox);

extern bool
findClosestIrrVoFs(Vector<VolIndex>& a_closestVoFs,
		   const VolIndex& a_vof,
		   const int a_dir,
		   const Side::LoHiSide& a_side,
		   const RealVect& a_dx,
		   const EBISBox& a_ebisBox,
		   const int a_maxNum,
		   const int a_maxDepth);
extern bool
findClosestIrrVoFs(Vector<VolIndex>& a_closestVoFs,
		   const VolIndex& a_vof,
		   const RealVect& a_dx,
		   const EBISBox& a_ebisBox,
		   const int a_maxNum,
		   const int a_maxDepth);

#include "NamespaceFooter.H"
#endif
