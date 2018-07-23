
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

#include "ParmParse.H"
#include "EBPetscMap.H"
#include "EBLevelDataOps.H"
#include "BaseEBCellFactory.H"

#include "NamespaceHeader.H"

EBPetscMap::EBPetscMap()
  :m_PetscOffsetIdx(0),
   m_isDefined(false),
   m_needsRemapping(true),
   m_finestLevel(-1)
{
#ifdef CH_MPI
    m_petsc_comm = Chombo_MPI::comm;
#else
    m_petsc_comm = PETSC_COMM_WORLD;
#endif
}

EBPetscMap::~EBPetscMap()
{
  clean();
}

void 
EBPetscMap::clean()
{
    m_isDefined = false;
    m_needsRemapping =true;
    m_finestLevel=-1;
    
    m_eblgsPtr.clear();
    m_refRatios.clear();
    m_PetscOffsetIdx = 0;
    
    for (int ilev=0;ilev<m_VofGlobalMap.size();ilev++)
    {
	if(m_VofGlobalMap[ilev]) delete m_VofGlobalMap[ilev];
    }
    m_VofGlobalMap.clear();
    
    for (int ilev=0;ilev<m_vofIterator.size();ilev++)
    {
      if(m_vofIterator[ilev]) delete m_vofIterator[ilev];
    }
    m_vofIterator.clear();
}
// define oject
void 
EBPetscMap::defineGrids( const int a_finestLev, 
			 const Vector<const EBLevelGrid*> &a_eblg,
			 const Vector<int> &a_refratios, 
			 const bool a_includeIrregular,
			 const bool a_includeRegular,
			 const bool a_includeCovered)
{
  CH_TIME("EBPetscMap::defineGrids");
  ParmParse pp;

  //clean it before redefining
  if (m_isDefined) clean();
  
  m_finestLevel=a_finestLev;
  m_numGhostCells = 4;
  pp.query("num_ghost", m_numGhostCells);

  int numLevels= a_finestLev+1;
  m_PetscVecLocalSize.resize(numLevels,0);

  m_VofGlobalMap.resize(numLevels,NULL);
  m_vofIterator.resize(numLevels, NULL);
  
  for (int ilev = 0; ilev <numLevels; ilev++)
  {
      m_VofGlobalMap[ilev]    = new LevelData<BaseEBCellFAB<int64_t> >();
      m_vofIterator[ilev]     = new LayoutData<VoFIterator>();
  }

  m_eblgsPtr      = a_eblg;
  m_refRatios     = a_refratios;

  
  m_isDefined = true;
  m_needsRemapping =true;

  for (int ilev = 0; ilev <numLevels; ilev++)
  {
      BaseEBCellFactory<int64_t> ebcellfactInt(m_eblgsPtr[ilev]->getEBISL());
      m_VofGlobalMap[ilev]->define(m_eblgsPtr[ilev]->getDBL(), 1,  m_numGhostCells*IntVect::Unit, ebcellfactInt);
      m_vofIterator[ilev] ->define(m_eblgsPtr[ilev]->getDBL());
  }
  resizeChomboPetscMapping(a_includeIrregular, a_includeRegular, a_includeCovered);
}

void 
EBPetscMap::defineGrids( const int a_finestLev, 
			 const Vector<EBLevelGrid> &a_eblg,
			 const Vector<int> &a_refratios, 
			 const bool a_includeIrregular,
			 const bool a_includeRegular,
			 const bool a_includeCovered)
{
    Vector<const EBLevelGrid*> eblgPtrs;
    for (int ilev = 0; ilev <a_eblg.size(); ilev++)
    {
	eblgPtrs.push_back(&a_eblg[ilev]);
    }
    
    defineGrids(a_finestLev, 
		eblgPtrs,
		a_refratios, 
		a_includeIrregular,
		a_includeRegular,
		a_includeCovered);
}


void 
EBPetscMap::resizeChomboPetscMapping(const bool a_includeIrregular,
				     const bool a_includeRegular,
				     const bool a_includeCovered)
{
  CH_assert(m_isDefined);

  //store cell types to include in petsc vector - affects puChomboInPetsc and puPetscInChombo
  m_includeIrregular = a_includeIrregular;
  m_includeRegular = a_includeRegular;
  m_includeCovered   =  a_includeCovered;

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      PetscInt count=0;
      m_PetscVecLocalSize[ilev]=0;
      const EBLevelGrid* eblevelgrid = m_eblgsPtr[ilev];
      const EBISLayout &ebislayout = eblevelgrid->getEBISL();
      DataIterator dit = m_eblgsPtr[ilev]->getDBL().dataIterator();
      int nbox=dit.size();
#pragma omp parallel for reduction(+:count)
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  const EBISBox& ebbox    = ebislayout[datInd]; 
	  const EBGraph& ebgraph  = ebbox.getEBGraph();
	  VoFIterator& vofit = (*m_vofIterator[ilev])[datInd];
	  const Box& grid = m_eblgsPtr[ilev]->getDBL().get(datInd);
	  
	  //cache vofIterator
	  if (a_includeRegular)
	  {
	      IntVectSet ivsTot(grid);
	      //subtract finer level from ivs
	      if(ilev  < m_finestLevel)
	      {
		  //subtract off finer level regions from ivs
		  const DisjointBoxLayout& finer =  m_eblgsPtr[ilev+1]->getDBL();
		  for(LayoutIterator lit = finer.layoutIterator(); lit.ok(); ++lit)
		  {
		      Box coarsenedFine = finer[lit()];
		      coarsenedFine.coarsen(m_refRatios[ilev]);
		      ivsTot -= coarsenedFine;
		  }
	      }
	      vofit.define(ivsTot, ebgraph);
	  }
	  else
	  {
	      IntVectSet ivsIrreg = ebbox.getIrregIVS(grid);
	      
	      //subtract off finer level regions from ivs
	      if(ilev  < m_finestLevel)
	      {
		  //subtract off finer level regions from ivs
		  const DisjointBoxLayout& finer =  m_eblgsPtr[ilev+1]->getDBL();
		  for(LayoutIterator lit = finer.layoutIterator(); lit.ok(); ++lit)
		  {
		      Box coarsenedFine = finer[lit()];
		      coarsenedFine.coarsen(m_refRatios[ilev]);
		      ivsIrreg -= coarsenedFine;
		  }
	      }
	      vofit.define(ivsIrreg, ebgraph);
	  }
	  
	  for(vofit.reset(); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      const IntVect& iv = vof.gridIndex();
	      //check flags whether include irregular and covered cells into pets vector
	      if (ebbox.isRegular(iv) && !a_includeRegular) continue;
	      if (ebbox.isIrregular(iv) && !a_includeIrregular) continue;
	      //if (ebbox.isCovered(iv) && !a_includeCovered) continue;
	      count++;
	      
	  }
      }
      if (count)
      {
	  m_PetscVecLocalSize[ilev] = count;
      }
  }
}

void 
EBPetscMap::setChomboPetscGlobalMap(Vec& a_vec,
				    ISLocalToGlobalMapping& a_petsc_mapping,
				    const int a_numComp)
{
  PetscErrorCode ierr;
  //setup a global mapping
  int64_t localSize    = getPetscVecBaseSize();
  int64_t globalSize = localSize;

#ifdef CH_MPI
  // const unsigned numProcs = numProc();
  // const unsigned rank = procID();
  // int64_t* allsizes = new int64_t[numProcs];
  // MPI_Allgather(&localSize, 1 , MPI_LONG_LONG_INT, allsizes, 1 , MPI_LONG_LONG_INT, Chombo_MPI::comm);  
  MPI_Allreduce(&localSize, &globalSize, 1, MPI_LONG_LONG_INT, MPI_SUM,  Chombo_MPI::comm);

  // int64_t testGlobalSize=0;
  // //let's test it here
  // for (int i=0;i<numProcs;i++)
  // {
  //     testGlobalSize += allsizes[i];
  // }
  
  // if (testGlobalSize !=globalSize)
  // {
  //     MayDay::Error("mismatch between global size and gathered sizes!!!");
  // }
#endif

//  pout()<<"globalsize="<<globalSize<<"\t localsize="<<localSize<<endl;
  
  if (a_vec) VecDestroy(&a_vec);
  VecCreateMPI(m_petsc_comm, a_numComp*localSize, a_numComp*globalSize, &a_vec);

  //save local to global mapping  
  PetscInt start = 0;
  PetscInt end   = a_numComp*localSize;
  ierr = VecGetOwnershipRange(a_vec,&start,&end); CHKERRV(ierr);
// #ifdef CH_MPI
//   for (int i=0;i<rank;i++)
//   {
//       start += allsizes[i]*a_numComp;
//   }
//   end = start + localSize*a_numComp;
//   delete [] allsizes;
// #endif
  
  CH_assert(a_numComp*localSize == (end-start));

  //setting local offset 
  m_PetscOffsetIdx = start;
  
  int64_t countGlobal = start;
  int64_t localIdx    = 0;
  PetscInt *indxsGlobal=NULL;
  if (localSize) indxsGlobal = new PetscInt[a_numComp*localSize];
  
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      const EBLevelGrid* eblevelgrid = m_eblgsPtr[ilev];
      const EBISLayout &ebislayout = eblevelgrid->getEBISL();
      DataIterator dit = m_eblgsPtr[ilev]->getDBL().dataIterator();
      int nbox=dit.size();
      //     #pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  (*m_VofGlobalMap[ilev])[datInd].setVal(-1);

	  const EBISBox& ebbox = ebislayout[datInd]; 
	  
	  VoFIterator& vofit = (*m_vofIterator[ilev])[datInd];
	  for(vofit.reset(); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      const IntVect& iv = vof.gridIndex();
	      //check flags whether include irregular and covered cells into pets vector
	      if (ebbox.isRegular(iv) && !m_includeRegular) continue;
	      if (ebbox.isIrregular(iv) && !m_includeIrregular) continue;
	      if (ebbox.isCovered(iv) && !m_includeCovered) continue;

	      (*m_VofGlobalMap[ilev])[datInd](vof,0) = countGlobal;
	      
	      for (PetscInt icomp=0; icomp<a_numComp;icomp++)
	      {
		  indxsGlobal[localIdx++] = countGlobal;
		  countGlobal++;
	      }
	  }
      }
      m_VofGlobalMap[ilev]->exchange();
  }
  ierr = ISLocalToGlobalMappingCreate(m_petsc_comm,1, a_numComp*localSize, indxsGlobal,PETSC_COPY_VALUES, &a_petsc_mapping); CHKERRV(ierr);
  VecSetLocalToGlobalMapping(a_vec, a_petsc_mapping);
  if(indxsGlobal) delete[] indxsGlobal;
}


//
// putChomboInPetsc
//
PetscErrorCode 
EBPetscMap::putChomboInPetsc(const Vector<LevelData<EBCellFAB>*>& a_ch, Vec& a_b)
{  
  CH_TIME("EBPetscMap::putChomboInPetsc");
  PetscErrorCode ierr;
  PetscScalar *arr;
  ierr = VecGetArray(a_b, &arr);CHKERRQ(ierr);

  //get number of components of ebcellfab 
  const PetscInt num_comp = (*a_ch[0]).nComp();
  const PetscInt localVecSize = num_comp*getPetscVecBaseSize();
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      const DisjointBoxLayout &dbl = a_ch[ilev]->getBoxes();
      DataIterator dit = dbl.dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  const EBCellFAB &FabEB = (*a_ch[ilev])[datInd];
	  VoFIterator& vofit = (*m_vofIterator[ilev])[datInd];
	  for(vofit.reset(); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      const PetscInt globalIdx = (*m_VofGlobalMap[ilev])[datInd](vof,0);
	      if (globalIdx < 0) continue;
	      PetscInt idx = globalIdx - m_PetscOffsetIdx; //shift for the global vector
	      CH_assert(idx>=0 && (idx+num_comp-1)<localVecSize);
              for (PetscInt icomp=0; icomp<num_comp;icomp++)
	      {
		  arr[idx] = FabEB(vof,icomp);
		  idx++;
	      }
	  }
      }
  }

  ierr =  VecRestoreArray(a_b,&arr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


PetscErrorCode 
EBPetscMap::putChomboInPetsc(const Vector<LevelData<BaseIVFAB<Real> >*>& a_ch, Vec& a_b)
{  
  CH_TIME("EBPetscMap::putChomboInPetsc");
  PetscErrorCode ierr;
  PetscScalar *arr;
  ierr = VecGetArray(a_b, &arr);CHKERRQ(ierr);

  //get number of components of ebcellfab 
  const PetscInt num_comp = (*a_ch[0]).nComp();
  const PetscInt localVecSize = num_comp*getPetscVecBaseSize();
  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      const DisjointBoxLayout &dbl = a_ch[ilev]->getBoxes();
      DataIterator dit = dbl.dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  const BaseIVFAB<Real>&FabIV = (*a_ch[ilev])[datInd];
	  
	  VoFIterator& vofit = (*m_vofIterator[ilev])[datInd];
	  for(vofit.reset(); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      const PetscInt globalIdx = (*m_VofGlobalMap[ilev])[datInd](vof,0);
	      if (globalIdx < 0) continue;
	      PetscInt idx = globalIdx - m_PetscOffsetIdx; //shift for the global vector
	      CH_assert(idx>=0 && (idx+num_comp-1)<localVecSize);
	      for (PetscInt icomp=0; icomp<num_comp;icomp++)
	      {
		  arr[idx] = FabIV(vof,icomp);
		  idx++;
	      }
	  }
      }
  }
  ierr =  VecRestoreArray(a_b,&arr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

//
// putPetscInChombo
PetscErrorCode 
EBPetscMap::putPetscInChombo(Vec a_x, Vector<LevelData<EBCellFAB>*>& a_ch) 
{
  CH_TIME("EBPetscMap::putPetscInChombo");

  PetscScalar *arr;
  PetscErrorCode ierr;
  ierr = VecGetArrayRead(a_x,(const PetscScalar**)&arr);CHKERRQ(ierr);

  //get number of components of ebcellfab 
  const PetscInt num_comp = (*a_ch[0]).nComp();
  const PetscInt localVecSize = num_comp*getPetscVecBaseSize();

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      const DisjointBoxLayout &dbl = a_ch[ilev]->getBoxes();
      DataIterator dit = dbl.dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  EBCellFAB &FabEB = (*a_ch[ilev])[datInd];
	  VoFIterator& vofit = (*m_vofIterator[ilev])[datInd];
	  for(vofit.reset(); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      const PetscInt globalIdx = (*m_VofGlobalMap[ilev])[datInd](vof,0);
	      if (globalIdx < 0) continue;
	      PetscInt idx = globalIdx - m_PetscOffsetIdx; //shift for the global vector
	      CH_assert(idx>=0 && (idx+num_comp-1)<localVecSize);
              for (PetscInt icomp=0; icomp<num_comp;icomp++)
	      {
		  FabEB(vof,icomp) = arr[idx];
		  idx++;
	      }
	  }
      }
      a_ch[ilev]->exchange();
  }
  ierr =  VecRestoreArrayRead(a_x,(const PetscScalar**)&arr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode 
EBPetscMap::putPetscInChombo(Vec a_x, Vector<LevelData<BaseIVFAB<Real> >*>& a_ch) 
{
  CH_TIME("EBPetscMap::putPetscInChombo");

  PetscScalar *arr;
  PetscErrorCode ierr;
  ierr = VecGetArrayRead(a_x,(const PetscScalar**)&arr);CHKERRQ(ierr);

  //get number of components of ebcellfab 
  const PetscInt num_comp = (*a_ch[0]).nComp();
  const PetscInt localVecSize = num_comp*getPetscVecBaseSize();

  for (int ilev = 0; ilev <= m_finestLevel; ilev++)
  {
      //#pragma omp parallel for
      const DisjointBoxLayout &dbl = a_ch[ilev]->getBoxes();
      DataIterator dit = dbl.dataIterator();
      int nbox=dit.size();
#pragma omp parallel for
      for (int mybox=0;mybox<nbox; mybox++)
      {
	  const DataIndex& datInd = dit[mybox];
	  BaseIVFAB<Real>&FabIV = (*a_ch[ilev])[datInd];
	  VoFIterator& vofit = (*m_vofIterator[ilev])[datInd];
	  for(vofit.reset(); vofit.ok(); ++vofit)
	  {
	      const VolIndex& vof = vofit();
	      const PetscInt globalIdx = (*m_VofGlobalMap[ilev])[datInd](vof,0);
	      if (globalIdx < 0) continue;
	      PetscInt idx = globalIdx - m_PetscOffsetIdx; //shift for the global vector
	      CH_assert(idx>=0 && (idx+num_comp-1)<localVecSize);
              for (PetscInt icomp=0; icomp<num_comp;icomp++)
	      {
		  FabIV(vof,icomp) = arr[idx];
		  idx++;
	      }
	  }
      }
      a_ch[ilev]->exchange();
  }
  ierr =  VecRestoreArrayRead(a_x,(const PetscScalar**)&arr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscInt 
EBPetscMap::getPetscVecBaseSize(const int a_lev) const
{
  PetscInt size=0;
  if (a_lev<0)
    {
      for (int ilev = 0; ilev <= m_finestLevel; ilev++)
	{
	  size +=m_PetscVecLocalSize[ilev];
	}
    }
  else
    {
      size =m_PetscVecLocalSize[a_lev];
    }
  return size;
}

#include "NamespaceFooter.H"
#endif
