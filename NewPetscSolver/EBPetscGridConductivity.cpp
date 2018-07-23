#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif
#ifdef CH_USE_PETSC
//#include "EBDebugOut.H"
#include "EBPetscGridConductivity.H"
#include "EBConductivityOp.H"
#include "BaseEBCellFactory.H"
#include "EBStencil.H"

#include "NamespaceHeader.H"

void EBPetscGridConductivity::setOp(LinearOp<LevelData<EBCellFAB> > *a_op)
{
    m_op = static_cast<EBConductivityOp*>(a_op);
}

PetscErrorCode
EBPetscGridConductivity::formPetscMat(Mat& a_petsc_mat,
				 const LevelData<EBCellFAB>& a_phi)
{
    CH_TIME("EBPetscGridConductivity::formPetscMat");
    const int lev = 0;
    PetscErrorCode ierr;

    //set to zero
    // ierr = MatZeroEntries(a_petsc_mat);

    //get number of components of ebcellfab 
    const PetscInt num_comp = a_phi.nComp();
    Real alpha,beta;
    m_op->getAlphaBeta(alpha,beta);
    LayoutData<BaseIVFAB<Real> > const*       alphaDiagWeight;
    m_op->getAlphaDiagWeight(alphaDiagWeight);
    LayoutData<BaseIVFAB<VoFStencil> >const* ebbcFluxStencil;
    m_op->getEBBCFluxStencil(ebbcFluxStencil);
    const RealVect vectDx = m_op->dx()*RealVect::Unit;
    const RefCountedPtr<BaseDomainBC>& domainBC = m_op->getDomainBC();
    const RefCountedPtr<LevelData<EBFluxFAB> >& bcoef = m_op->getBScalingCoefficients();
    //const RefCountedPtr<LevelData<EBCellFAB> >& acoef = m_op->getAScalingCoefficients();
    const EBLevelGrid& eblg = m_op->getEBLG();
    const ProblemDomain& domain = eblg.getDomain();

    DataIterator dit = a_phi.dataIterator();
    int nbox=dit.size();
    for (int mybox=0;mybox<nbox; mybox++)
    {
	const DataIndex& datInd = dit[mybox];
	const EBCellFAB& phiFab = a_phi[datInd];
	const EBISBox& ebisBox = phiFab.getEBISBox();
	const BaseIVFAB<Real>& alphaDiag = (*alphaDiagWeight)[datInd];

	VoFIterator& vofit = (*m_vofIterator[lev])[datInd];
	for(vofit.reset(); vofit.ok(); ++vofit)
	{
	    const VolIndex& vof = vofit();
	    const IntVect& iv = vof.gridIndex();
	    const PetscInt row_idx = (*m_VofGlobalMap[lev])[datInd](vof,0);
	    if(row_idx<0) continue;

	    //off diagonal
	    std::vector<PetscScalar> vals;
	    std::vector<PetscInt> cols;
	    VoFStencil opStencil;
	    
	    m_op->getDivFStencil(opStencil, vof, datInd);

	    if (ebisBox.isIrregular(iv) && ebbcFluxStencil !=NULL)
	    {    		
		opStencil += (*ebbcFluxStencil)[datInd](vof,0);
	    }
		
	    const int size = opStencil.size();
	    for (PetscInt ivof = 0; ivof < size; ivof++)
	    {
		const VolIndex& curVoF = opStencil.vof(ivof);
		PetscInt col_idx =(*m_VofGlobalMap[lev])[datInd](curVoF,0);
		if(col_idx>=0)
		{
		    const Real offDiag = opStencil.weight(ivof) * beta;
		    for (PetscInt icomp=0; icomp<num_comp;icomp++)
		    {
			{
			    vals.push_back(offDiag);
			    cols.push_back(col_idx);
			}
			col_idx++;
		    }
		}
	    }

	    //handle domainBC
	    const RealVect bdx2 = 2./(vectDx*vectDx) * beta;
	    for( int idir = 0 ; idir < CH_SPACEDIM ; ++idir)
	    {
		for (SideIterator sit; sit.ok(); ++sit)
		{
		    int sign = 1;
		    if  (sit() == Side::Lo) sign =-1;
		    const IntVect conIV = iv+sign*BASISV(idir);
		    //skip if vof is not boundary
		    if (domain.contains(conIV)) continue;
		    
		    Vector<FaceIndex> faces = ebisBox.getFaces(vof, idir, sit());
		    for (int iface=0; iface<faces.size(); iface++)
		    {
			const FaceIndex& face = faces[iface];
			if (face.isBoundary())
			{
			    const bool isDirichlet = domainBC->isDirichletDom(vof, VolIndex(conIV,0),phiFab);
			    if (isDirichlet)
			    {
				//adjust diagonal value
				{
				    Real bcValue = -bdx2[idir];
				    bcValue *= (*bcoef)[datInd][idir](face,0)*ebisBox.areaFrac(face);
				    
				    vals.push_back(bcValue);
				    cols.push_back(row_idx);
				}
			    }else
			    {
				//no need for Neumann domainBC adjustment
			    }
			}
		    }//iface
		}//sideIterator for domainBC
	    }//idir
	    
	    
	    //	diagonal part
	    if (alpha)
	    {
		const Real diagWeight = alpha* alphaDiag(vof,0);
		vals.push_back(diagWeight);
		cols.push_back(row_idx);
	    }
	    
	    for (PetscInt icomp=0; icomp<num_comp;icomp++)
	    {
		const PetscInt row = row_idx+icomp;
		ierr = MatSetValues(a_petsc_mat,1, &row, cols.size(), &(cols[0]), &(vals[0]), ADD_VALUES); CHKERRQ(ierr);
	    }
	}//vofIterator
    }//ditIterator
    
    
    ierr = MatAssemblyBegin(a_petsc_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
    ierr = MatAssemblyEnd(a_petsc_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);

    PetscFunctionReturn(0);
}

PetscErrorCode
EBPetscGridConductivity::addBCToRHS(Vec a_petsc_vec,
				    const LevelData<EBCellFAB>& a_phi,
				    const bool a_homogeneous)
{
    CH_TIME("EBPetscGridConductivity::addBCToRHS");

    if (a_homogeneous) PetscFunctionReturn(0);
    
    PetscErrorCode ierr;
    PetscScalar *arr;
    ierr = VecGetArray(a_petsc_vec, &arr);CHKERRQ(ierr);
    //we need to monkey with fluxes from domainBC to get a adjustment for rhs
    const int lev = 0;
    //get number of components of ebcellfab 
    const PetscInt num_comp = a_phi.nComp();
    const EBLevelGrid& eblg = m_op->getEBLG();
    
    const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
    
    Real alpha,beta;
    m_op->getAlphaBeta(alpha,beta);
    const RefCountedPtr<BaseDomainBC>& domainBC = m_op->getDomainBC();
    const RealVect vectDx = m_op->dx()*RealVect::Unit;
    const RefCountedPtr<LevelData<EBFluxFAB> >& bcoef = m_op->getBScalingCoefficients();
    //const RefCountedPtr<LevelData<EBCellFAB> >& acoef = m_op->getAScalingCoefficients();
    
    const RealVect idx1      = beta/vectDx;
    const RealVect idx2      = beta/(vectDx*vectDx);

    DataIterator dit = a_phi.dataIterator();
    int nbox=dit.size();
    for (int mybox=0;mybox<nbox; mybox++)
    {
	const DataIndex& datInd = dit[mybox];
	const EBCellFAB& phi = a_phi[datInd];
	const BaseFab<Real>& phiFab = phi.getSingleValuedFAB();
	const EBISBox& ebisBox = phi.getEBISBox();
	const Box& grid = dbl.get(datInd);
	//IntVectSet ivsTot(grid);
	RealVect origin = RealVect::Zero;
	
	for (int idir=0; idir<SpaceDim; idir++)
	{
	    Box loBox,hiBox;
	    int hasLo,hasHi;
	    EBArith::loHi(loBox, hasLo,
			  hiBox, hasHi,
			  eblg.getDomain(),grid, idir);

	    if (hasLo == 1)
            {
              Box lbox=loBox;
              lbox.shift(idir,-1);
              BaseFab<Real> loFaceFlux(loBox,num_comp);

              domainBC->getFaceFlux(loFaceFlux,phiFab, origin, vectDx,idir,Side::Lo,datInd, 0.,a_homogeneous);
	       
	       for(BoxIterator bit(lbox); bit.ok(); bit.next())
	       {
		   IntVect iv = bit()+BASISV(idir);
		   if (ebisBox.isRegular(iv) && !m_includeRegular) continue;
		   if (ebisBox.isIrregular(iv) && !m_includeIrregular) continue;
		   if (ebisBox.isCovered(iv) && !m_includeCovered) continue;
		   VolIndex vof(iv,0);
		   PetscInt idx = (*m_VofGlobalMap[lev])[datInd](vof,0) - m_PetscOffsetIdx; //shift for the global vector

		   const bool isDirichlet = domainBC->isDirichletDom(vof, VolIndex(bit(),0),phi);
		   for (int icomp = 0; icomp<num_comp; icomp++)
		   {
		       const FaceIndex face = FaceIndex(vof,VolIndex(bit(),0),idir);
		       arr[idx] += loFaceFlux(iv)*idx1[idir]*ebisBox.areaFrac(face);
		       if (isDirichlet)
		       {
			   arr[idx] += -2.*phi(vof,icomp)*idx2[idir]*(*bcoef)[datInd][idir](face,0)*ebisBox.areaFrac(face);
		       }
		       idx++;
		   }
	       }
	    }
	    
	    if (hasHi == 1)
	    {
		Box hbox=hiBox;
		hbox.shift(idir,1);
		BaseFab<Real> hiFaceFlux(hiBox,num_comp);

		domainBC->getFaceFlux(hiFaceFlux,phiFab,origin,vectDx,idir,Side::Hi,datInd,0.,a_homogeneous);
	       
	       for(BoxIterator bit(hbox); bit.ok(); bit.next())
	       {
		   IntVect iv = bit()-BASISV(idir);
		   if (ebisBox.isRegular(iv) && !m_includeRegular) continue;
		   if (ebisBox.isIrregular(iv) && !m_includeIrregular) continue;
		   if (ebisBox.isCovered(iv) && !m_includeCovered) continue;
		   VolIndex vof(iv,0);

		   PetscInt idx = (*m_VofGlobalMap[lev])[datInd](vof,0) - m_PetscOffsetIdx; //shift for the global vector

		   
		   const bool isDirichlet = domainBC->isDirichletDom(vof, VolIndex(bit(),0),phi);
		   for (int icomp = 0; icomp<num_comp; icomp++)
		   {
		       const FaceIndex face = FaceIndex(vof,VolIndex(bit(),0),idir);
		       arr[idx] += -hiFaceFlux(iv)*idx1[idir]*ebisBox.areaFrac(face);
		       if (isDirichlet)
		       {
		       	   arr[idx] += -2.0*phi(vof,icomp)*idx2[idir]*(*bcoef)[datInd][idir](face,0)*ebisBox.areaFrac(face);
		       }
		       idx++;
		   }
	       }
            }
	}
    }
    
    ierr =  VecRestoreArray(a_petsc_vec,&arr);CHKERRQ(ierr);

    PetscFunctionReturn(0);
}
	  
PetscErrorCode
EBPetscGridConductivity::getMatAIJSizes(PetscInt* a_d_nnz,
				   PetscInt* a_o_nnz,
				   const PetscInt a_nComp)
{
    CH_TIME("EBPetscGridConductivity::getMatAIJSizes");
    //const int numGhostCells = 4;
	
    const int lev = 0;
    //PetscErrorCode ierr;
    
    const EBLevelGrid& eblg = *m_eblgsPtr[lev];
    const DisjointBoxLayout &dbl = eblg.getDBL();
    //const ProblemDomain& domain = eblg.getDomain();
    LayoutData<BaseIVFAB<VoFStencil> >const* ebbcFluxStencil;
    m_op->getEBBCFluxStencil(ebbcFluxStencil);

    LevelData<BaseEBCellFAB<int> > procIDGlobal;
    BaseEBCellFactory<int> baseEBFact(eblg.getEBISL());
    procIDGlobal.define(dbl, 1, m_numGhostCells*IntVect::Unit, baseEBFact);
    DataIterator dit = dbl.dataIterator();
    const int nbox=dit.size();
    for (int mybox=0;mybox<nbox; mybox++)
    {
	const DataIndex& datInd = dit[mybox];
	procIDGlobal[datInd].setVal(procID());
    }
    const Interval interval(0, 0);
    procIDGlobal.exchange(interval);
    
    
    const int n = getPetscVecBaseSize()*a_nComp;
    std::fill_n(a_d_nnz, n, 0);
    std::fill_n(a_o_nnz, n, 0);


    for (int mybox=0;mybox<nbox; mybox++)
    {
	const DataIndex& datInd = dit[mybox];
	const EBISBox& ebisBox = eblg.getEBISL()[datInd];
	
	VoFIterator& vofit = (*m_vofIterator[lev])[datInd];
	for(vofit.reset(); vofit.ok(); ++vofit)
	{
	    const VolIndex& vof = vofit();
	    const IntVect& iv = vof.gridIndex();
	    const PetscInt row_idx = (*m_VofGlobalMap[lev])[datInd](vof,0);
	    if(row_idx<0) continue;

	    PetscInt loc_idx = (*m_VofGlobalMap[lev])[datInd](vof,0) - m_PetscOffsetIdx; //shift for the global vector

	    a_d_nnz[loc_idx]++;

	    VoFStencil opStencil;
	    m_op->getDivFStencil(opStencil, vof, datInd);
	    
	    if (ebisBox.isIrregular(iv) && ebbcFluxStencil !=NULL)
	    {    		
		opStencil += (*ebbcFluxStencil)[datInd](vof,0);
	    }
	    
	    const int64_t size = opStencil.size();
	    for (PetscInt ivof = 0; ivof < size; ivof++)
	    {
		const VolIndex& curVoF = opStencil.vof(ivof);
		if (curVoF==vof) continue;
		const int procid = procIDGlobal[datInd](curVoF,0);
		if (procid==procID()) a_d_nnz[loc_idx]++;
		else a_o_nnz[loc_idx]++;
	    }
	
	    // pout()<<iv<<"-"<<row_idx<<"\t d_nnz="<<a_d_nnz[loc_idx]<<"\t o_nnz="<<a_o_nnz[loc_idx]<<"\t isirr= "<<ebisBox.isIrregular(iv)<<endl;
	}//vofIterator
    }//ditIterator
    
    PetscFunctionReturn(0);
}

#include "NamespaceFooter.H"
#endif
