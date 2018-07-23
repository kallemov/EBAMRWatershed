#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
//#include <algorithm>
#include "EBArith.H"
#include "FaceIterator.H"
#include "VoFIterator.H"
#include "EBSUBSurfaceOps.H"
#include "EBSUBSurfaceOpsF_F.H"
#include "PolyGeom.H"

//const static Real regConstant = 1.e-2;
const static int s_nghostTemp=0;
#include "NamespaceHeader.H"


void getSaturation(EBCellFAB&        a_F,
		   const EBCellFAB& a_psi,
		   const Real a_Ssat,
		   const Real a_Sres,
		   const Real a_alpha,
		   const Real a_nsoil,
		   const ProblemDomain& a_domain,
		   const Box&       a_box)
{
    const EBISBox& ebisBox = a_psi.getEBISBox();
    IntVectSet ivsTot(a_box);
    for (VoFIterator vofit(ivsTot, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
	const VolIndex& vof = vofit();
	Real& Freal    = a_F(vof, 0);
	const Real psireal    = a_psi(vof, 0);
	
	FORT_GETSATURATIONPOINT(CHF_REAL(Freal),
				CHF_CONST_REAL(psireal),
				CHF_CONST_REAL(a_Ssat),
				CHF_CONST_REAL(a_Sres),
				CHF_CONST_REAL(a_alpha),
				CHF_CONST_REAL(a_nsoil));
    }
}

void addSaturationTerms(EBCellFAB&        a_F,
			const EBCellFAB& a_psi,
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
			const Box&       a_box)
{
  CH_TIME("EBSUBSurfaceOps::addSaturationTerms");

 
  BaseFab<Real>& F_vec              = a_F.getSingleValuedFAB();
  const BaseFab<Real>& psi            = a_psi.getSingleValuedFAB();
  const BaseFab<Real>& psiold         = a_psi_old.getSingleValuedFAB();
  const BaseFab<Real>& Ss           = a_Ss.getSingleValuedFAB();
  const BaseFab<Real>& Phi          = a_Phi.getSingleValuedFAB();
  //  const BaseFab<Real>& rho          = a_Rho.getSingleValuedFAB();

  //these cells have unit volume fraction
  FORT_ADDSATURATIONTERMS(CHF_FRA1(F_vec,0),
			  CHF_CONST_FRA1(psi,0),
			  CHF_CONST_FRA1(psiold,0),
			  CHF_CONST_FRA1(Ss,0),
			  CHF_CONST_FRA1(Phi,0),
			  // CHF_CONST_FRA1(rho,0),
			  CHF_CONST_REAL(rho),
			  CHF_CONST_REAL(a_Ssat),
			  CHF_CONST_REAL(a_Sres),
			  CHF_CONST_REAL(a_alpha),
			  CHF_CONST_REAL(a_nsoil),
			  CHF_BOX(a_box));
  //pointwise operation so just have to do the multi valued cells
  //as an irregular iteration
  const EBISBox& ebisBox = a_psi.getEBISBox();
  IntVectSet multiIVS = ebisBox.getMultiCells(a_box);
  for (VoFIterator vofit(multiIVS, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
      const VolIndex& vof = vofit();
      Real& Freal    = a_F(vof, 0);
      Real psireal    = a_psi(vof, 0);
      Real psioldreal = a_psi_old(vof, 0);
      Real Ssreal   = a_Ss(vof, 0);
      Real Phireal   = a_Phi(vof, 0);
      //Real Rhoreal   = a_Rho(vof, 0);
      // const Real volFrac = ebisBox.volFrac(vof);

      Real addterm;
      FORT_GETSATURATIONTERMSPOINT(CHF_REAL(addterm),
				   CHF_CONST_REAL(psireal),
				   CHF_CONST_REAL(psioldreal),
				   CHF_CONST_REAL(Ssreal),
				   CHF_CONST_REAL(Phireal),
				   //CHF_CONST_REAL(Rhoreal),
				   CHF_CONST_REAL(rho),
				   CHF_CONST_REAL(a_Ssat),
				   CHF_CONST_REAL(a_Sres),
				   CHF_CONST_REAL(a_alpha),
				   CHF_CONST_REAL(a_nsoil));

      
      Freal +=addterm;
    }

  //the DarcyOp evaluate it as kappa*L(h), then we have to multiply saturation term by volume fraction
  //now multiply all irregular cells by its volume fraction 

  // IntVectSet irrIVS = ebisBox.getIrregIVS(a_box);
  // for (VoFIterator vofit(irrIVS, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
  // {
  //     const VolIndex& vof = vofit();
  //     a_F(vof, 0) *=ebisBox.volFrac(vof);
  // }
}


void irregSaturationTerms(BaseIVFAB<Real>&   a_F,
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
			  const Box&       a_box)
{
    CH_TIME("EBSUBSurfaceOps::addSaturationTerms");
    
    
    const EBISBox& ebisBox = a_psi_old.getEBISBox();
    IntVectSet ivsIrreg = ebisBox.getIrregIVS(a_box);
    for (VoFIterator vofit(ivsIrreg, ebisBox.getEBGraph()); vofit.ok(); ++vofit)
    {
	const VolIndex& vof = vofit();
	Real& Freal    = a_F(vof, 0);
	Real psireal    = a_psi(vof, 0);
	Real psioldreal = a_psi_old(vof, 0);
	Real Ssreal   = a_Ss(vof, 0);
	Real Phireal   = a_Phi(vof, 0);
	//Real Rhoreal   = a_Rho(vof, 0);
	// const Real volFrac = ebisBox.volFrac(vof);
	
	Real addterm;
	FORT_GETSATURATIONTERMSPOINT(CHF_REAL(addterm),
				     CHF_CONST_REAL(psireal),
				     CHF_CONST_REAL(psioldreal),
				     CHF_CONST_REAL(Ssreal),
				     CHF_CONST_REAL(Phireal),
				     //CHF_CONST_REAL(Rhoreal),
				     CHF_CONST_REAL(rho),
				     CHF_CONST_REAL(a_Ssat),
				     CHF_CONST_REAL(a_Sres),
				     CHF_CONST_REAL(a_alpha),
				     CHF_CONST_REAL(a_nsoil));
	
	
	Freal +=addterm;
	// Freal *=ebisBox.volFrac(vof);
    }
}


void ccpHarmonicAverageToFaces(EBFaceFAB &             a_face,
			       const EBCellFAB &       a_cell,
			       const EBGraph &         a_ebGraph,
			       const Box &             a_grid,
			       const int &             a_idir,
			       const ProblemDomain &   a_domain,
			       const int &             a_comp)
{
  CH_TIME("EBSUBSurfaceOps::ccpHarmonicAverageToFaces");
  FaceStop::WhichFaces stopCrit;
  if (a_domain.isPeriodic(a_idir))
    {
      stopCrit = FaceStop::SurroundingWithBoundary;
    }
  else
    {
      stopCrit = FaceStop::SurroundingNoBoundary;
    }
  //initially set to zero so the boundary faces will all be zero
  a_face.setVal(0.);

  BaseFab<Real> &       regFace = a_face.getSingleValuedFAB();
  const BaseFab<Real> & regCell = a_cell.getSingleValuedFAB();

  //only want non-domain boundary faces
  Box faceBox = a_grid;
  faceBox.grow(a_idir, 1);
  faceBox  &= a_domain;
  faceBox.grow(a_idir, -1);
  faceBox.surroundingNodes(a_idir);
  //do regular cells in fortran
  FORT_CCPHARMMEANCELLTOFACE(CHF_FRA1(regFace, 0),
			     CHF_CONST_FRA1(regCell, a_comp),
			     CHF_CONST_INT(a_idir),
			     CHF_BOX(faceBox));
  
  //fix up irregular cells.  Only multiValued cells should and their neighbors
  //need fixing
  IntVectSet ivsIrreg = a_ebGraph.getMultiCells(a_grid);
  ivsIrreg.grow(1);
  ivsIrreg  &= a_grid;
  for (FaceIterator faceit(ivsIrreg, a_ebGraph, a_idir, stopCrit); faceit.ok(); ++faceit)
    {
      const FaceIndex & face = faceit();

      a_face(face, 0) = 2.0/(1.0/a_cell(face.getVoF(Side::Hi), a_comp) +
			     1.0/a_cell(face.getVoF(Side::Lo), a_comp));
      // pout()<<a_face(face, 0)<<"\n";
    }
}


void ccpHarmonicAverageToFaces(LevelData<EBFluxFAB> &         a_Flux,
			       const LevelData<EBCellFAB> &   a_CCcell,
			       const DisjointBoxLayout &      a_grids,
			       const EBISLayout &             a_ebisl,
			       const ProblemDomain &          a_domain,
			       const int                      a_comp,
			       const int                      a_lev)
{
  CH_TIME("EBSUBSurfaceOps::ccpHarmonicAverageToFaces(LevelData)");
  Interval interv(a_comp,a_comp);
  LevelData<EBCellFAB>& CCval = (LevelData<EBCellFAB>&) a_CCcell;
  CCval.exchange(interv);
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      Box grid;
      if (a_lev>0)
      {
	  grid = grow(a_grids.get(dit()), s_nghostTemp);
	  grid &= a_domain;
      }
      else
      {
	  grid = a_grids.get(dit());
      }
     
      //const Box grid = a_grids.get(dit());
      // Box grid = grow(a_grids.get(dit()), 1);
      // grid &= a_domain;

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB &       faceVal = a_Flux[dit()][idir];
          const EBCellFAB & cellVal = a_CCcell[dit()];

          ccpHarmonicAverageToFaces(faceVal, cellVal, ebgraph, grid, idir, a_domain, a_comp);
        }
    }
}


void ccpHarmonicAverageToFaces(LevelData<EBFluxFAB> &         a_Flux,
			       const LevelData<EBCellFAB> &   a_CCcell,
			       const DisjointBoxLayout &      a_grids,
			       const EBISLayout &             a_ebisl,
			       const ProblemDomain &          a_domain,
    			       const int                      a_lev)
{
  CH_TIME("EBSUBSurfaceOps::ccpHarmonicAverageToFaces(LevelData)");
  // Interval interv(a_comp,a_comp);
  // LevelData<EBCellFAB>& CCval = (LevelData<EBCellFAB>&) a_CCcell;
  // a_CCcell.exchange();
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
       Box grid;
      if (a_lev>0)
      {
	  grid = grow(a_grids.get(dit()), s_nghostTemp);
	  grid &= a_domain;
      }
      else
      {
	  grid = a_grids.get(dit());
      }

      //const Box grid = a_grids.get(dit());
      // Box grid = grow(a_grids.get(dit()), 1);
      // grid &= a_domain;

      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB &       faceVal = a_Flux[dit()][idir];
          const EBCellFAB & cellVal = a_CCcell[dit()];

          ccpHarmonicAverageToFaces(faceVal, cellVal, ebgraph, grid, idir, a_domain, idir);
        }
    }
}


void scaleLevelFluxByKr(LevelData<EBFluxFAB> &         a_Flux,
			const LevelData<EBCellFAB> &   a_psi,
			const LevelData<EBFluxFAB> &   a_darcy,
			const Real a_alpha,
			const Real a_nsoil,
			const Real a_lambda,
			const DisjointBoxLayout &      a_grids,
			const EBISLayout &             a_ebisl,
			const ProblemDomain &          a_domain,
			const int& a_method,
			const int                      a_lev)
{
  CH_TIME("EBSUBSurfaceOps::scaleLevelFluxByKr");
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      // Box grid;
      // if (a_lev>0)
      // {
      // 	  grid = grow(a_grids.get(dit()), s_nghostTemp);
      // 	  grid &= a_domain;
      // }
      // else
      // {
      // 	  grid = a_grids.get(dit());
      // }
      const Box& grid = a_grids.get(dit());
      // Box grid = grow(a_grids.get(dit()), 2);
      // grid &= a_domain;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          EBFaceFAB &       faceVal = a_Flux[dit()][idir];
          const EBFaceFAB & darcyVal = a_darcy[dit()][idir];
          const EBCellFAB & psifab    = a_psi[dit()];

          scaleFaceByKr(faceVal, psifab, darcyVal,
			a_alpha, a_nsoil, a_lambda,
			ebgraph, grid, idir, a_domain, 
			a_method);
	  
        }
    }
}

void scaleFaceByKr(EBFaceFAB &             a_face,
		   const EBCellFAB &       a_psi,
		   const EBFaceFAB &       a_darcy,
		   const Real a_alpha,
		   const Real a_nsoil,
		   const Real a_lambda,
		   const EBGraph &         a_ebGraph,
		   const Box &             a_grid,
		   const int &             a_idir,
		   const ProblemDomain &   a_domain,
		   const int &             a_method)
{
  CH_TIME("EBSUBSurfaceOps::scaleFaceByKr");
  const Real fluxTol=1e-16;

  FaceStop::WhichFaces stopCrit;
  if (a_domain.isPeriodic(a_idir))
    {
      stopCrit = FaceStop::SurroundingWithBoundary;
    }
  else
    {
      stopCrit = FaceStop::SurroundingNoBoundary;
    }

  BaseFab<Real> &       regFace   = a_face.getSingleValuedFAB();
  const BaseFab<Real> & psireg    = a_psi.getSingleValuedFAB();

  //only want non-domain boundary faces
  Box faceBox = a_grid;
  faceBox.grow(a_idir, 1);
  faceBox  &= a_domain;
  faceBox.grow(a_idir, -1);
  faceBox.surroundingNodes(a_idir);
  //do regular cells in fortran
  //harmonic averages
  if (a_method==0)
    {
      //do harmonic averaging of Kr(h)
      FORT_MULTBYHARMMEANKR(CHF_FRA1(regFace, 0),
			    CHF_CONST_FRA1(psireg, 0),
			    CHF_CONST_REAL(a_alpha),
			    CHF_CONST_REAL(a_nsoil),
			    CHF_CONST_REAL(a_lambda),
			    CHF_CONST_INT(a_idir),
			    CHF_BOX(faceBox));


      //fix up irregular cells.  Only multiValued cells should and their neighbors
      //need fixing
      IntVectSet ivsIrreg = a_ebGraph.getMultiCells(a_grid);
      ivsIrreg.grow(1);
      ivsIrreg  &= a_grid;
      for (FaceIterator faceit(ivsIrreg, a_ebGraph, a_idir, stopCrit); faceit.ok(); ++faceit)
	{
	  const FaceIndex & face = faceit();
	  const Real& psiHi = a_psi(face.getVoF(Side::Hi), 0);
	  const Real& psiLo = a_psi(face.getVoF(Side::Lo), 0);
	  Real Krhi, Krlo;  

	  FORT_RELATIVEPERMEABILITY(CHF_REAL(Krhi),
				    CHF_CONST_REAL(psiHi),
				    CHF_CONST_REAL(a_alpha),
				    CHF_CONST_REAL(a_nsoil),
				    CHF_CONST_REAL(a_lambda));

	  FORT_RELATIVEPERMEABILITY(CHF_REAL(Krlo),
				    CHF_CONST_REAL(psiLo),
				    CHF_CONST_REAL(a_alpha),
				    CHF_CONST_REAL(a_nsoil),
				    CHF_CONST_REAL(a_lambda));


	  a_face(face, 0) *= 2.0/(1.0/Krhi + 1.0/Krlo);
	}

    }
  //upwind 
  else  if (a_method==1)
    {
      stopCrit = FaceStop::SurroundingNoBoundary;
      //need fixing
      IntVectSet ivsTot(a_grid);
      //ivsTot.grow(1);
      //ivsTot &= a_grid;
      for (FaceIterator faceit(ivsTot, a_ebGraph, a_idir, stopCrit); faceit.ok(); ++faceit)
	{
	  const FaceIndex & face = faceit();
	  Real psi, psi2;
	  bool isZero=false;


	  if (fabs(a_darcy(face,0))<fluxTol) 
	  {
	    psi  = a_psi(face.getVoF(Side::Hi), 0);
	    psi2 = a_psi(face.getVoF(Side::Lo), 0);
	    isZero=true;
	  }
    	  else if (a_darcy(face,0)>0.0) 
	  {
	    psi= a_psi(face.getVoF(Side::Lo), 0);
	  }
	  else if (a_darcy(face,0)<0.0) 
	  {
	    psi = a_psi(face.getVoF(Side::Hi), 0);
	  }
 	     
	  Real Kr, Kr2;  
	      
	  FORT_RELATIVEPERMEABILITY(CHF_REAL(Kr),
				    CHF_CONST_REAL(psi),
				    CHF_CONST_REAL(a_alpha),
				    CHF_CONST_REAL(a_nsoil),
				    CHF_CONST_REAL(a_lambda));
	  if (isZero)
	    {
	      FORT_RELATIVEPERMEABILITY(CHF_REAL(Kr2),
					CHF_CONST_REAL(psi2),
					CHF_CONST_REAL(a_alpha),
					CHF_CONST_REAL(a_nsoil),
					CHF_CONST_REAL(a_lambda));

	      a_face(face, 0) *= 2.0/(1.0/Kr + 1.0/Kr2);
	    }
	  else
	    {
	      a_face(face, 0) *=Kr;
	    }
	  //  pout()<<a_idir<<"-"<<face<<"\t"<<side<<"\t"<<a_darcy(face,0)<<"\t"<<Kr<<endl;
	}
      // pout()<<"****************"<<endl;

    }
  else
    {
      MayDay::Error("EBSUBSurfaceOps::scaleFaceByKr:: unsupported method of averaging Kr to Faces");
    }
}

Real getRelativePermeability(const Real a_psi,
			     const Real a_alpha,
			     const Real a_nsoil,
			     const Real a_lambda)
{
  Real Kr;
  FORT_RELATIVEPERMEABILITY(CHF_REAL(Kr),
			    CHF_CONST_REAL(a_psi),
			    CHF_CONST_REAL(a_alpha),
			    CHF_CONST_REAL(a_nsoil),
			    CHF_CONST_REAL(a_lambda));
  return Kr;
}

Real getDerivRelPermeability(const Real a_psi,
			     const Real a_alpha,
			     const Real a_nsoil,
			     const Real a_lambda)
{
    Real ah = fabs(a_psi*a_nsoil);
    Real ahn= pow(ah,a_nsoil);
    Real inv_ahn_one = pow((1. + ahn), 1./a_nsoil );

    Real derKr= -(0.5/(ah*ah*fabs(a_psi)))
	*ahn
	*sqrt( inv_ahn_one /  pow( (1.+ahn), 7. ) ) 
	*(-ah * (1. + ahn) + (-4. + ahn)*inv_ahn_one)
        *(-ah + ahn * (-ah + inv_ahn_one))
	*(-1. + a_nsoil);

  return derKr;
}

Real getDerivSaturation(const Real a_psi,
			const Real a_Ssat,
			const Real a_Sres,
			const Real a_alpha,
			const Real a_nsoil)
{
    Real ah = fabs(a_psi*a_nsoil);
    Real ahn= pow(ah,a_nsoil);
    Real ahn_one = (1. + ahn);

    Real derSw = 
      ahn*pow(ahn_one, 1./a_nsoil)/pow(ahn_one, 2)*(-1. + a_nsoil)*(a_Sres - a_Ssat)/fabs(a_psi);
    
    return derSw;
}

void
ccpAverageFaceToCells(EBCellFAB &             a_cellData,
                      const EBFluxFAB &       a_fluxData,
                      const EBGraph &         a_ebGraph,
                      const Box &             a_grid,
                      const ProblemDomain &   a_domain,
                      const RealVect &        a_dx)
{
  IntVectSet ivsIrreg = a_ebGraph.getIrregCells(a_grid);
  int icomp = 0;
  
  const EBISBox&  ebisBox = a_cellData.getEBISBox();
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      const EBFaceFAB & faceData = a_fluxData[idir];
      const BaseFab<Real> & regFaceData =   faceData.getSingleValuedFAB();
      BaseFab<Real> &       regCellData = a_cellData.getSingleValuedFAB();
      FORT_CCPAVEFACETOCELL( CHF_FRA1(regCellData, idir),
                             CHF_CONST_FRA1(regFaceData, icomp),
                             CHF_CONST_INT(idir),
                             CHF_BOX(a_grid));

      for (VoFIterator vofit(ivsIrreg, a_ebGraph); vofit.ok(); ++vofit)
        {
          const VolIndex & vof = vofit();
          int numFaces = 0;
          Real cellVal = 0.;
          Real sumAreaFrac = 0.;
          for (SideIterator sit; sit.ok(); ++sit)
            {
              Vector<FaceIndex> faces = a_ebGraph.getFaces(vof, idir, sit());
              //if we have faces, then use them.  otherwise, need to extrapolate to covered
              if (faces.size() > 0)
                {
                  for (int iface = 0; iface < faces.size(); iface++)
                    {
                      Real areaFrac = ebisBox.areaFrac(faces[iface]);
                      // cellVal += areaFrac*faceData(faces[iface], icomp);
                      cellVal += faceData(faces[iface], icomp);
                      sumAreaFrac += areaFrac;
                      numFaces++;
                    }
                }
              // else
              //   {
              //     Real extrapValCov = ccpGetCoveredExtrapValue(vof, idir, sit(),
              //                                                  faceData,
              //                                                  ebisBox,
              //                                                  a_grid,
              //                                                  a_domain,
              //                                                  a_dx,
              //                                                  icomp);


              //     //add in covered face value as if it were a normal value
              //     cellVal += extrapValCov;
              //     numFaces++;
              //   }
            }
          if (numFaces > 0)
            {
              cellVal /= Real(numFaces);
              // cellVal /= sumAreaFrac;
              a_cellData(vof, idir) = cellVal;
            }
          else
            {
              a_cellData(vof, idir) = 0;
            }
        }
    } //end loop over directions
}

void ccpLevelccpAverageFaceToCells(LevelData<EBCellFAB>&       a_cellData,
				   const LevelData<EBFluxFAB>& a_fluxData,
				   const DisjointBoxLayout &      a_grids,
				   const EBISLayout &             a_ebisl,
				   const ProblemDomain &          a_domain,
				   const RealVect &        a_dx,
				   const int                      a_lev)
{
  for (DataIterator dit = a_grids.dataIterator(); dit.ok(); ++dit)
    {
      const EBGraph & ebgraph = a_ebisl[dit()].getEBGraph();
      // Box grid;
      // if (a_lev>0)
      // {
      // 	  grid = grow(a_grids.get(dit()), s_nghostTemp);
      // 	  grid &= a_domain;
      // }
      // else
      // {
      // 	  grid = a_grids.get(dit());
      // }
      const Box& grid = a_grids.get(dit());
      //for (int idir = 0; idir < SpaceDim; idir++)
        {
	  EBCellFAB& cellFab       = a_cellData[dit()];
          const EBFluxFAB& fluxFab = a_fluxData[dit()];
         
	  ccpAverageFaceToCells(cellFab,
				fluxFab,
				ebgraph,
				grid,
				a_domain,
				a_dx);
	  
	}
    }
}


void copyBIVF2EBCellFAB(const Vector<LevelData<BaseIVFAB<Real> >*>& a_bivf, 
			Vector<LevelData<EBCellFAB>*>& a_ebcellf,
			const Vector<EBLevelGrid>&        a_eblg)
{  
  CH_TIME("copyBIVF2EBCelFAB");
  int ncomp = a_bivf[0]->nComp();

  for (int ilev = 0; ilev < a_ebcellf.size(); ilev++)
  {
      const DisjointBoxLayout& dbl = a_eblg[ilev].getDBL();
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
	  EBCellFAB &FabEB = (*a_ebcellf[ilev])[dit()];
	  const EBISBox& ebisBox = a_eblg[ilev].getEBISL()[dit()];
	  const EBGraph& ebgraph = ebisBox.getEBGraph();
	  const Box& grid = dbl.get(dit());
	  IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
	  for(VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      for (int comp = 0; comp < ncomp; comp++)
	      {
		  FabEB(vof,comp) = (*a_bivf[ilev])[dit()](vof,comp);
	      }
	  }
      }
      a_ebcellf[ilev]->exchange();
  }
}

void copyBIVF2EBCellFAB(const Vector<LevelData<BaseIVFAB<Real> >*>& a_bivf, 
			Vector<LevelData<EBCellFAB>*>& a_ebcellf,
			const Vector<const EBLevelGrid*>&        a_eblg)
{
  CH_TIME("copyBIVF2EBCelFAB");
  int ncomp = a_bivf[0]->nComp();

  for (int ilev = 0; ilev < a_ebcellf.size(); ilev++)
  {
      const DisjointBoxLayout& dbl = a_eblg[ilev]->getDBL();
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
	  EBCellFAB &FabEB = (*a_ebcellf[ilev])[dit()];
	  const EBISBox& ebisBox = a_eblg[ilev]->getEBISL()[dit()];
	  const EBGraph& ebgraph = ebisBox.getEBGraph();
	  const Box& grid = dbl.get(dit());
	  IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
	  for(VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      for (int comp = 0; comp < ncomp; comp++)
	      {
		  FabEB(vof,comp) = (*a_bivf[ilev])[dit()](vof,comp);
	      }
	  }
      }
      a_ebcellf[ilev]->exchange();
  }
}

void copyEBCellFAB2BIVF(const Vector<LevelData<EBCellFAB>*>& a_ebcellf,
			Vector<LevelData<BaseIVFAB<Real> >*>& a_bivf,
			const Vector<EBLevelGrid>&        a_eblg)
{  
  CH_TIME("copyEBCelFAB2BIVF");

  int ncomp = a_ebcellf[0]->nComp();
  for (int ilev = 0; ilev < a_ebcellf.size(); ilev++)
  {
      const DisjointBoxLayout& dbl = a_eblg[ilev].getDBL();
      for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit)
      {
	  const EBCellFAB &FabEB = (*a_ebcellf[ilev])[dit()];
	  const EBISBox& ebisBox = a_eblg[ilev].getEBISL()[dit()];
	  const EBGraph& ebgraph = ebisBox.getEBGraph();
	  const Box& grid = dbl.get(dit());
	  IntVectSet ivsIrreg = ebisBox.getIrregIVS(grid);
	  for(VoFIterator vofit(ivsIrreg, ebgraph); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      for (int comp = 0; comp < ncomp; comp++)
	      {
		  (*a_bivf[ilev])[dit()](vof,comp) = FabEB(vof,comp);
	      }
	  }
      }
      a_bivf[ilev]->exchange();
  }
}

RealVect getAnisotropicNormal(const RealVect& a_normal,
			      const RealVect& a_vectDx)
{
    RealVect anisotNormal;
    
    for (int idir = 0; idir < SpaceDim; idir++)
    {
	anisotNormal[idir] = a_normal[idir]/a_vectDx[idir];
    }

    anisotNormal /=anisotNormal.vectorLength();
    return anisotNormal;
}

Real getAnisotropicCorrection(const VolIndex& a_vof ,
			      const RealVect& a_vectDx,
			      const EBISBox& a_ebisBox)
{
    const Real areaFrac = a_ebisBox.bndryArea(a_vof);
    const RealVect normal = a_ebisBox.normal(a_vof);
    Real anisotr_corr=0.0;
    for (int idir = 0; idir < SpaceDim; idir++)
    {
	anisotr_corr += normal[idir]*normal[idir] 
	    /(a_vectDx[idir]*a_vectDx[idir]);
    }
    return areaFrac*sqrt(anisotr_corr);
}


bool  findClosestIrrVoF(VolIndex& a_vof,
			const Vector<VolIndex>& a_VoFs,
			const RealVect& a_center,
			const RealVect& a_dx,
			const EBISBox& a_ebisBox)
{
    if (!a_VoFs.size()) return false;

    const Real Tol = 1e-15;
    bool foundOne = false;
    Real minDist=1e+20;
    
    for (int ivof = 0; ivof < a_VoFs.size(); ivof++)
    {
	const VolIndex& curVoF = a_VoFs[ivof];
	if (curVoF==a_vof) continue;
	if (a_ebisBox.bndryArea(curVoF) >= Tol)
	{
	    const RealVect coordLo = RealVect::Unit*0.5+curVoF.gridIndex()+a_ebisBox.bndryCentroid(curVoF);
	    RealVect dist = (coordLo-a_center);
	    dist *= a_dx;
	    
	    if (dist.vectorLength()<minDist)
	    {
		foundOne = true;
		a_vof = curVoF;
		minDist = fabs(dist.vectorLength());
	    }
	}
	
    }
    return foundOne;
}

bool  findClosestIrrVoF(VolIndex& a_vof,
			const Vector<VolIndex>& a_VoFs,
			const RealVect& a_center,
			const int a_dir,
			const Side::LoHiSide& a_side,
			const RealVect& a_dx,
			const EBISBox& a_ebisBox)
{
    if (!a_VoFs.size()) return false;

    const Real Tol = 1e-15;
    bool foundOne = false;
    Real minDist=1e+20;
    
    for (int ivof = 0; ivof < a_VoFs.size(); ivof++)
    {
	const VolIndex& curVoF = a_VoFs[ivof];
	
	if (a_ebisBox.bndryArea(curVoF) >= Tol)
	{
	    const RealVect coordLo = RealVect::Unit*0.5+curVoF.gridIndex()+a_ebisBox.bndryCentroid(curVoF);
	    //const RealVect coordLo = RealVect::Unit*0.5+curVoF.gridIndex()+a_ebisBox.centroid(curVoF);
	    Real dist = coordLo[a_dir]-a_center[a_dir];
	    //dist *= a_dx;
	    
	    if ((a_side == Side::Lo && dist <-Tol) || (a_side == Side::Hi && dist > Tol))
	    {
		if (fabs(dist)<minDist)
		{
		    foundOne = true;
		    a_vof = curVoF;
		    minDist = fabs(dist);
		}
	    }
	}
	
    }
    return foundOne;
}


bool  findClosestIrrVoFs(Vector<VolIndex>& a_closestVoFs,
			 const VolIndex& a_vof,
			 const int a_dir,
			 const Side::LoHiSide& a_side,
			 const RealVect& a_dx,
			 const EBISBox& a_ebisBox,
			 const int a_maxNum,
			 const int a_maxDepth)
{
    a_closestVoFs.clear();
    const Real Tol = 1e-15;
    if (a_ebisBox.bndryArea(a_vof) < Tol) return false;
    bool foundOne = false;
    Vector<VolIndex> vofList;
    EBArith::getAllVoFsWithinRadius (vofList, a_vof, a_ebisBox, a_maxDepth);
    RealVect centCoord = RealVect(a_vof.gridIndex());
    centCoord += a_ebisBox.bndryCentroid(a_vof);

    VolIndex* closestVoFs = new VolIndex[a_maxNum];
    Real* distances = new Real[a_maxNum];
    std::fill_n(distances, a_maxNum, 1e+123);
    for (int ivof = 0; ivof < vofList.size(); ivof++)
    {
	const VolIndex& curVoF = vofList[ivof];
	if (curVoF==a_vof) continue;
	
	if (a_ebisBox.bndryArea(curVoF) >= Tol)
	{
	    RealVect coord = (RealVect(curVoF.gridIndex())+a_ebisBox.bndryCentroid(curVoF))- centCoord;
	    coord *= a_dx;
	    //we need distance in 2D
	    coord[2]=0.0;
	    const Real dist = coord.vectorLength();
	    //if ((a_side == Side::Lo && coord[a_dir] <= 0.25) || (a_side == Side::Hi && coord[a_dir] >=-0.25))
	    {
		for (int c=0;c<a_maxNum;c++)
		{
		    if (dist<distances[c])
		    {
			foundOne = true;
			for (int cc=a_maxNum-1;cc>c;cc--)
			{
			    closestVoFs[cc]=closestVoFs[cc-1];
			    distances[cc]=distances[cc-1];
			}
			closestVoFs[c]=curVoF;
			distances[c]=dist;
			break;
		    }
		}
	    }
	}
    }
    for (int c=0;c<a_maxNum;c++)
    {
	if (!closestVoFs[c].isDefined()) break;
	a_closestVoFs.push_back(closestVoFs[c]);
    }
    delete[] closestVoFs;
    delete[] distances;
    return foundOne;
}


bool findClosestIrrVoFs(Vector<VolIndex>& a_closestVoFs,
			const VolIndex& a_vof,
			const RealVect& a_dx,
			const EBISBox& a_ebisBox,
			const int a_maxNum,
			const int a_maxDepth)
{
    a_closestVoFs.clear();
    const Real Tol = 1e-15;
    if (a_ebisBox.bndryArea(a_vof) < Tol) return false;
    bool foundOne = false;
    Vector<VolIndex> vofList;
    EBArith::getAllVoFsWithinRadius (vofList, a_vof, a_ebisBox, a_maxDepth);
    RealVect centCoord = RealVect(a_vof.gridIndex());
    centCoord += a_ebisBox.bndryCentroid(a_vof);

    VolIndex* closestVoFs = new VolIndex[a_maxNum];
    Real* distances = new Real[a_maxNum];
    std::fill_n(distances, a_maxNum, 1e+123);
    for (int ivof = 0; ivof < vofList.size(); ivof++)
    {
	const VolIndex& curVoF = vofList[ivof];
	if (curVoF==a_vof) continue;
	
	if (a_ebisBox.bndryArea(curVoF) >= Tol)
	{
	    RealVect coord = (RealVect(curVoF.gridIndex())+a_ebisBox.bndryCentroid(curVoF))- centCoord;
	    coord *= a_dx;
	    //we need distance in 2D
	    coord[2] = 0.0;
	    const Real dist = coord.vectorLength();
	    
	    for (int c=0;c<a_maxNum;c++)
	    {
		if (dist<distances[c])
		{
		    foundOne = true;
		    for (int cc=a_maxNum-1;cc>c;cc--)
		    {
			closestVoFs[cc]=closestVoFs[cc-1];
			distances[cc]=distances[cc-1];
		    }
		    closestVoFs[c]=curVoF;
		    distances[c]=dist;
		    break;
		}
	    }
	}
    }
    
    for (int c=0;c<a_maxNum;c++)
    {
	if (!closestVoFs[c].isDefined()) break;
	a_closestVoFs.push_back(closestVoFs[c]);
    }
    delete[] closestVoFs;
    delete[] distances;
    return foundOne;
}

#include "NamespaceFooter.H"

