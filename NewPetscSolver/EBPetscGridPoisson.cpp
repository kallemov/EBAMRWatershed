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
#include "EBDebugOut.H"
#include "EBPetscGridPoisson.H"
#include "BaseEBCellFactory.H"

#include "NamespaceHeader.H"

void EBPetscGridPoisson::setOp(LinearOp<LevelData<EBCellFAB> > *a_op)
{
    m_op = static_cast<EBAMRPoissonOp*>(a_op);
}

PetscErrorCode
EBPetscGridPoisson::formPetscMat(Mat& a_petsc_mat,
				 const LevelData<EBCellFAB>& a_phi)
{
    const int lev = 0;
    PetscErrorCode ierr;

    //set to zero
    // ierr = MatZeroEntries(a_petsc_mat);

    //get number of components of ebcellfab 
    const PetscInt num_comp = a_phi.nComp();
    const EBLevelGrid& eblg = *m_eblgsPtr[lev];
    const ProblemDomain& domain = eblg.getDomain();
    Real alpha,beta;
    m_op->getAlphaBeta(alpha,beta);
    const RealVect vectDx = m_op->dx()*RealVect::Unit;
    const RefCountedPtr<BaseDomainBC>& domainBC = m_op->getDomainBC();

    //get stencil and alphaDiagWeights for irregular cells
    LayoutData<BaseIVFAB<VoFStencil> > const* opIrrVoFStencil;
    LayoutData<BaseIVFAB<Real> > const*       alphaDiagWeight;
    m_op->getVoFStencil(opIrrVoFStencil);
    m_op->getAlphaDiagWeight(alphaDiagWeight);
    
    DataIterator dit = a_phi.dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
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

	    std::vector<PetscScalar> vals;
	    std::vector<PetscInt> cols;
	    
	    if (ebisBox.isRegular(iv))
	    {    
		const RealVect idx2      = 1./(vectDx*vectDx) * beta;
		const RealVect diagValueVec = -2.0 * idx2;
		const Real diagValue = diagValueVec.sum() + alpha;

		//setting standard stencil for Laplace operator

		//off-diagonal
		for( int idir = 0 ; idir < CH_SPACEDIM ; ++idir)
		{
		    const IntVect lIV = iv-BASISV(idir);
		    const IntVect rIV = iv+BASISV(idir);
		    PetscInt col_idx_left = -1;
		    // const Box& gridWithGhosts = ebisBox.getRegion();
		    // IntVectSet ivsWithGhosts(gridWithGhosts);
		    if (domain.contains(lIV))
		    {
			col_idx_left  = (*m_VofGlobalMap[lev])[datInd](VolIndex(lIV,0),0);
		    }
		    PetscInt col_idx_right = -1;
		    if (domain.contains(rIV))
		    {
			col_idx_right = (*m_VofGlobalMap[lev])[datInd](VolIndex(rIV,0),0);
		    }
		    if (col_idx_left>=0)
		    {
			for (PetscInt icomp=0; icomp<num_comp;icomp++)
			{
			    {
				vals.push_back(idx2[idir]);
				cols.push_back(col_idx_left);
			    }
			    col_idx_left++;
			}
		    }
		    if (col_idx_right>=0)
		    {
			for (PetscInt icomp=0; icomp<num_comp;icomp++)
			{
			    {
				vals.push_back(idx2[idir]);
				cols.push_back(col_idx_right);
			    }
			    col_idx_right++;
			}
		    }
		    

		    //now handle domainBC
		    //the value for phi of Dirichlet or flux value for Neumann
		    //go directly to rhs in addBCToRHS function
		    //matrix needs adjustments too
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
			    if (faces[iface].isBoundary())
			    {
				const bool isDirichlet = domainBC->isDirichletDom(vof, VolIndex(conIV,0),phiFab);
				if (isDirichlet)
				{
				    //adjust diagonal value
				    //for (PetscInt icomp=0; icomp<num_comp;icomp++)
				    {
					const Real bcValue = -idx2[idir];
					//MatSetValue(a_petsc_mat,row_idx+icomp,row_idx+icomp, bcValue ,ADD_VALUES);
					vals.push_back(bcValue);
					cols.push_back(row_idx);

				    }
				}else
				{
				    //Neumann adjustment
				    //for (PetscInt icomp=0; icomp<num_comp;icomp++)
				    {
					const Real bcValue = idx2[idir];
					//MatSetValue(a_petsc_mat,row_idx+icomp,row_idx+icomp, bcValue ,ADD_VALUES);
					vals.push_back(bcValue);
					cols.push_back(row_idx);
				    }
				}
			    }
			}//iface
		    }//sideIterator for domainBC
		}//idir
		
		//diagonal part
		//for (PetscInt icomp=0; icomp<num_comp;icomp++)
		{
		    //diagonal value
		    //ierr = MatSetValue(a_petsc_mat,row_idx+icomp,row_idx+icomp,diagValue,ADD_VALUES); CHKERRQ(ierr);
		    vals.push_back(diagValue);
		    cols.push_back(row_idx);
		}
	    }
	    else //now irregular
	    {
		const VoFStencil& opStencil      = (*opIrrVoFStencil)[datInd](vof,0);
		const Real diagWeight = alpha* alphaDiag(vof,0);

		//this is something strange to get matching with EBPetscLinearSolver
		// {
		//     const RealVect idx2      = 1./(vectDx*vectDx) * beta;
		//     const RealVect diagValueVec = -2.0 * idx2;
		//     const Real diagValue = diagValueVec.sum() + alpha;
		//     Real nil=  1e-16 * (abs(diagValue)+abs(idx2[0])) * ( diagValue < 0 ? -1 : 1 );
		//     //for (PetscInt icomp=0; icomp<num_comp;icomp++)
		//     {
		// 	// MatSetValue(a_petsc_mat,row_idx+icomp,row_idx+icomp, nil, ADD_VALUES);
		// 	vals.push_back(nil);
		// 	cols.push_back(row_idx);
		//     }
		// }
		
		const int size = opStencil.size();

		for (PetscInt ivof = 0; ivof < size; ivof++)
		{
		    const VolIndex& curVoF = opStencil.vof(ivof);
		    PetscInt col_idx =(*m_VofGlobalMap[lev])[datInd](curVoF,0);
		    // if (vof.gridIndex()==IntVect(58,16))
		    // {
		    // 	pout()<<setprecision(16);
		    // 	pout()<<"new curVoF="<< curVoF;
		    // }	
		    
		    if(col_idx>=0)
		    {
			const Real offDiag = opStencil.weight(ivof) * beta;
			//if (vof.gridIndex()==IntVect(58,16)) pout()<<"\t value="<< offDiag;

			for (PetscInt icomp=0; icomp<num_comp;icomp++)
			{
			    //for (PetscInt icomp2=0; icomp2<num_comp;icomp2++)
			    {

				//ierr = MatSetValue(a_petsc_mat,row_idx+icomp2,col_idx, offDiag, ADD_VALUES); CHKERRQ(ierr);
				vals.push_back(offDiag);
				cols.push_back(col_idx);
			    }
			    col_idx++;
			}
		    }
		    //if (vof.gridIndex()==IntVect(58,16)) pout()<< endl;
		}


/*
		//now handle domainBC
		VoFStencil domainBCStencil = VoFStencil();
		m_op->getDomainFluxStencil(domainBCStencil,
					 vof,
					 0,
					 datInd);
				
				
		const int size2 = domainBCStencil.size();
		if (size2>0)
		{
		    const Real idx2      = 1./(vectDx[0]*vectDx[0]) * beta;
		    domainBCStencil *= beta;
		    //domainBCStencil *= ebisBox.areaFracScaling(vof);

		    //adjust diagonal value
		    for (PetscInt icomp=0; icomp<num_comp;icomp++)
		    {
			MatSetValue(a_petsc_mat,row_idx+icomp,row_idx+icomp,idx2,ADD_VALUES);
		    }
		    for (PetscInt ivof = 0; ivof < size2; ivof++)
		    {
			const VolIndex& curVoF = domainBCStencil.vof(ivof);
			PetscInt col_idx =(*m_VofGlobalMap[lev])[datInd](curVoF,0);
			if(col_idx>=0)
			{
			    const Real offDiag = domainBCStencil.weight(ivof);
			    for (PetscInt icomp=0; icomp<num_comp;icomp++)
			    {
				for (PetscInt icomp2=0; icomp2<num_comp;icomp2++)
				{
				    
				    MatSetValue(a_petsc_mat,row_idx+icomp2,col_idx, offDiag, ADD_VALUES);
				}
				col_idx++;
			    }
			}
		    }
		}
*/
		//	diagonal part
		//for (PetscInt icomp=0; icomp<num_comp;icomp++)
		{
		    // if (vof.gridIndex()==IntVect(58,16))
		    // {
		    // 	pout()<<"new Diagvalue="<< diagWeight<<endl;
		    // }
			
		    //ierr = MatSetValue(a_petsc_mat,row_idx+icomp, row_idx+icomp, diagWeight, ADD_VALUES); CHKERRQ(ierr);
		    vals.push_back(diagWeight);
		    cols.push_back(row_idx);
		}
	    }
	    for (PetscInt icomp=0; icomp<num_comp;icomp++)
	    {
		const PetscInt row = row_idx+icomp;
#pragma omp critical
		{
		    ierr = MatSetValues(a_petsc_mat,1, &row, cols.size(), &(cols[0]), &(vals[0]), ADD_VALUES); 
		}
	    }
	}//vofIterator
    }//ditIterator
    
    
    ierr = MatAssemblyBegin(a_petsc_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
    ierr = MatAssemblyEnd(a_petsc_mat,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);

    // const int num = numProc();
    // char iter_str[80];
    // sprintf(iter_str, "A%03d.m",num);
    // PetscViewer viewer;
    // PetscViewerASCIIOpen( PETSC_COMM_WORLD, iter_str, &viewer);
    // PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
    // MatView( a_petsc_mat, viewer );
    // PetscViewerDestroy(&viewer);

    PetscFunctionReturn(0);
}

PetscErrorCode
EBPetscGridPoisson::addBCToRHS(Vec a_petsc_vec,
			       const LevelData<EBCellFAB>& a_phi,
			       const bool a_homogeneous)
{
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
    const RealVect vectDx = m_op->dx()*RealVect::Unit;
    const RefCountedPtr<BaseDomainBC>& domainBC = m_op->getDomainBC();
    const RealVect idx1      = 1./vectDx * beta;
    const RealVect idx2      = 1./(vectDx*vectDx) * beta;

    DataIterator dit = a_phi.dataIterator();
    int nbox=dit.size();
#pragma omp parallel for
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
		       arr[idx] += loFaceFlux(iv)*idx1[idir];
		       if (isDirichlet)
		       {
			   arr[idx] += -2.*phi(vof,icomp)*idx2[idir];
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
		       arr[idx] += -hiFaceFlux(iv)*idx1[idir];
		       if (isDirichlet)
		       {
			   arr[idx] += -2.0*phi(vof,icomp)*idx2[idir];
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
EBPetscGridPoisson::getMatAIJSizes(PetscInt* a_d_nnz,
				   PetscInt* a_o_nnz,
				   const PetscInt a_nComp)
{
    //const int numGhostCells = 4;
	
    const int lev = 0;
    //PetscErrorCode ierr;
    
    const EBLevelGrid& eblg = *m_eblgsPtr[lev];
    const DisjointBoxLayout &dbl = eblg.getDBL();
    const ProblemDomain& domain = eblg.getDomain();
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

    //get stencil and alphaDiagWeights for irregular cells
    LayoutData<BaseIVFAB<VoFStencil> > const* opIrrVoFStencil;
    m_op->getVoFStencil(opIrrVoFStencil);
    
    
#pragma omp parallel 
 {
    PetscInt* private_d_nnz;
    PetscInt* private_o_nnz;
    
#ifdef _OPENMP
    const int numthreads=omp_get_num_threads();
    if(numthreads>1)
    {
	private_d_nnz = new PetscInt[n];
	private_o_nnz = new PetscInt[n];
	std::fill_n(private_d_nnz, n, 0);
	std::fill_n(private_o_nnz, n, 0);
    }
    else
    {
	private_d_nnz = a_d_nnz;
	private_o_nnz = a_o_nnz;
    }
	
#else
    private_d_nnz = a_d_nnz;
    private_o_nnz = a_o_nnz;
#endif
    
#pragma omp for
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

	    private_d_nnz[loc_idx]++;
	    
	    if (ebisBox.isRegular(iv))
	    {
		//off-diagonal
		for( int idir = 0 ; idir < CH_SPACEDIM ; ++idir)
		{
		    const IntVect lIV = iv-BASISV(idir);
		    const IntVect rIV = iv+BASISV(idir);
		    {
			VolIndex curVoF(lIV,0);
			if (domain.contains(lIV))
			{
			    const int procid = procIDGlobal[datInd](curVoF,0);
			    if (procid==procID()) private_d_nnz[loc_idx]++;
			    else a_o_nnz[loc_idx]++;
			}
		    }
		    {
			VolIndex curVoF(rIV,0);
			if (domain.contains(rIV))
			{
			    const int procid = procIDGlobal[datInd](curVoF,0);
			    if (procid==procID()) private_d_nnz[loc_idx]++;
			    else a_o_nnz[loc_idx]++;
			}
		    }
		    
		}//idir
	    }
	    else //now irregular
	    {
		const VoFStencil& opStencil      = (*opIrrVoFStencil)[datInd](vof,0);

		const int64_t size = opStencil.size();
		for (PetscInt ivof = 0; ivof < size; ivof++)
		{
		    const VolIndex& curVoF = opStencil.vof(ivof);
		    if (curVoF==vof) continue;
		    const int procid = procIDGlobal[datInd](curVoF,0);
		    if (procid==procID()) private_d_nnz[loc_idx]++;
		    else private_o_nnz[loc_idx]++;
		}
	    }
	    // pout()<<iv<<"-"<<row_idx<<"\t d_nnz="<<a_d_nnz[loc_idx]<<"\t o_nnz="<<a_o_nnz[loc_idx]<<"\t isirr= "<<ebisBox.isIrregular(iv)<<endl;
	}//vofIterator
    }//ditIterator

#ifdef _OPENMP
    if(numthreads>1)
    {
#pragma omp critical
	for (int i=0;i<n; i++)
	{
	    a_d_nnz += private_d_nnz[i];
	    a_o_nnz += private_o_nnz[i];	
	}
	delete[] private_d_nnz;
	delete[] private_o_nnz;
    }
#endif
    
  }//omp parallel   
    PetscFunctionReturn(0);
}

#include "NamespaceFooter.H"
#endif
