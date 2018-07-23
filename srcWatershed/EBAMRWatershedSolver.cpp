#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "ParmParse.H"
//#include "EBDebugOut.H"
#include "EBLevelDataOps.H"
#include "EBAMRDataOps.H"
#include "EBAMRWatershedSolver.H"
#include "petscsnes.h"
#include "petscksp.h"

#include "UsingNamespace.H"

    
EBAMRWatershedSolver::EBAMRWatershedSolver(EBPetscCompGridRichards *a_richardsCompGrid, WatershedIBC *a_ibc)
  : m_richardsCompGrid(a_richardsCompGrid),
    m_ibc(a_ibc),
    m_isMatrixFree(true),
    m_isFD(false),
    m_isInitialized(false),
    m_petsc_mat(NULL), 
    m_petsc_x(PETSC_NULL), 
    m_petsc_b(PETSC_NULL), 
    m_petsc_snes(PETSC_NULL),
    m_petsc_ksp(PETSC_NULL),
    m_isDefined(false),
    m_remove_const_nullspace(false)
{
#ifdef CH_MPI
    m_petsc_comm = Chombo_MPI::comm;
#else
    m_petsc_comm = PETSC_COMM_WORLD;
#endif
    
    ParmParse pp;
    
    m_snes_typeID = 1;
    pp.query("SNESType", m_snes_typeID);
    
    m_isMatrixFree = true;
    pp.query("MatrixFree", m_isMatrixFree);
    
    m_isFD = false; 
    pp.query("JacobianFD", m_isFD);
    if (m_isFD && m_isMatrixFree) 
    {
	pp.query("MFFD_type", m_matmffd_type);
	m_mffd_ds_umin=1.0e-6;
	pp.query("MFFD_DS_umin", m_mffd_ds_umin);
	m_mffd_wp_normu=false;
	pp.query("MFFD_WP_normU", m_mffd_wp_normu);
    }
    m_fd_coloring = true;
    if (m_isFD) 
    {
	pp.query("JacobianFD_coloring", m_fd_coloring);
    }    
    string PetscOptionsFile;
    pp.query("petsc_options_file",PetscOptionsFile);
#if PETSC_VERSION_GE(3,7,0)
    PetscOptionsInsertFile(m_petsc_comm, PETSC_NULL,PetscOptionsFile.c_str(), PETSC_TRUE);
#else
    PetscOptionsInsertFile(m_petsc_comm, PetscOptionsFile.c_str(), PETSC_TRUE);
#endif
    
    pp.query("remove_const_nullspace", m_remove_const_nullspace);

    m_presolver = false; 
    pp.query("use_presolver", m_presolver);
    if (m_presolver) 
    {
	pp.get("num_presolver_iterations", m_presolver_iterations);
    }

} 

EBAMRWatershedSolver::~EBAMRWatershedSolver()
{
    m_isDefined=false;
    deallocateSolver();
}

void 
EBAMRWatershedSolver::deallocateSolver()
{
    m_isInitialized=false;
    
    if (m_petsc_mat) MatDestroy(&m_petsc_mat);
    m_petsc_mat = PETSC_NULL;
    
    //  if ( m_petsc_ksp ) KSPDestroy(&m_petsc_ksp);
    //  m_petsc_ksp = NULL;
    
    if ( m_petsc_snes ) SNESDestroy(&m_petsc_snes);
    m_petsc_snes = PETSC_NULL;
    
    if ( m_petsc_x ) VecDestroy(&m_petsc_x);
    m_petsc_x = PETSC_NULL;
}

PetscErrorCode
EBAMRWatershedSolver::formVecFunctionWatershed(SNES snes,Vec x,Vec f,void * ctx)
{
    CH_TIME("EBAMRWatershedSolver::formVecFunctionWatershed");
    EBAMRWatershedSolver* solver = static_cast<EBAMRWatershedSolver*>(ctx);

    // PetscReal norm;
   // VecNorm(x, NORM_2, &norm);
   // pout()<<"xnorm="<<norm<<endl;

   // VecView(x, 	PETSC_VIEWER_STDOUT_WORLD);

    if (solver->m_include_SUBsurface)	
    {
	EBAMRDataOps::setToZero(solver->m_subsur_mfree);
	solver->m_richardsCompGrid->putPetscInChombo(x,solver->m_subsur_mfree);
	//solver->m_richardsCompGrid->averageDown(solver->m_subsur_mfree);
	solver->m_richardsCompGrid->averageDownHydrolic(solver->m_subsur_mfree);
    }
    else
    {
	solver->m_richardsCompGrid->putPetscInChombo(x,solver->m_sur_mfree);
	//solver->m_richardsCompGrid->floorIrregular(solver->m_sur_mfree);
    }

    // static int count=0;
    // char fileChar[1000];
    // sprintf(fileChar, "vectorXX.step.%d.hdf5", count);
    // writeEBAMRname(&solver->m_subsur_mfree, fileChar);
    
    if (!solver->m_include_SUBsurface)	
    {
	solver->m_richardsCompGrid->addAMRIrregExchangeFlux(solver->m_sur_mfree, solver->m_ibc->getBoundarySource(), solver->m_dt);
    }
    else
    {
	if (solver->m_ibc->isSurfaceSolverOn())
	{
	    solver->m_richardsCompGrid->addAMRIrregExchangeFlux(solver->m_subsur_mfree, solver->m_ibc->getBoundarySource(), solver->m_dt);
	}
	
    }
    // sprintf(fileChar, "vectorQQ.step.%d.hdf5", count);
    //  writeIVLevel(&(*(solver->m_richardsCompGrid->m_spreadFlux)[0]));//, fileChar);
    
    if (solver->m_include_SUBsurface)	
    {
	solver->m_richardsCompGrid->getNewtonKrylovVectorOp(solver->m_Fsubsur_mfree, solver->m_subsur_mfree, solver->m_vectDx, solver->m_dt);
    }
    else
    {
	solver->m_richardsCompGrid->getNewtonKrylovVectorOp(solver->m_Fsur_mfree, solver->m_sur_mfree, solver->m_vectDx, solver->m_dt);
    }
    
    // sprintf(fileChar, "vectorFF.step.%d.hdf5", count);
    // writeEBAMRname(&solver->m_Fsubsur_mfree, fileChar);
    // count++;

    //wrapping petsc vector back
    if (solver->m_include_SUBsurface)	
    {
	solver->m_richardsCompGrid->putChomboInPetsc(solver->m_Fsubsur_mfree,f);
    }
    else
    {
	solver->m_richardsCompGrid->putChomboInPetsc(solver->m_Fsur_mfree,f);
    }

    
    // VecNorm(f, NORM_INFINITY, &norm);
    // pout()<<"fnorm="<<norm<<endl;
    // VecView(f, 	PETSC_VIEWER_STDOUT_WORLD);

    PetscFunctionReturn(0);
}

PetscErrorCode
EBAMRWatershedSolver::formJacobianFunctionWatershed(SNES snes,Vec f,Mat J,Mat P,void* ctx)
{
    CH_TIME("EBAMRWatershedSolver::formJacobianFunctionWatershed");
    EBAMRWatershedSolver* solver = static_cast<EBAMRWatershedSolver*>(ctx);
    
    // Vec x;
    // VecDuplicate(solver->m_petsc_x, &x);
    // VecCopy(solver->m_petsc_x, x);
    
    solver->m_richardsCompGrid->putPetscInChombo(solver->m_petsc_x, solver->m_subsur_mfree);
    //  VecDestroy(&x);

    std::vector<PetscInt>  rows;
    std::vector<std::vector<PetscInt> > collumns; 
    std::vector<std::vector<PetscScalar> > data; 
    solver->m_richardsCompGrid->getJacobianVectors(rows,
						   collumns,
						   data,
						   solver->m_subsur_mfree, 
						   solver->m_dt);
    
    PetscInt n;
    VecGetLocalSize(f, &n);
    CH_assert(n==rows.size());
    
    for (int i=0;i<rows.size();i++)
    {
	const PetscInt num_col =collumns[i].size();
	{
	    //insert stencil terms
	    MatSetValues(J, 1, &(rows[i]), num_col, &(collumns[i][0]), &(data[i][0]), INSERT_VALUES);
	}
    }
    
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);
    

    PetscFunctionReturn(0);
}
  
PetscErrorCode
EBAMRWatershedSolver::formColoringJacobian(SNES snes,Vec f,Mat J,Mat P,void* ctx)
{
    CH_TIME("EBAMRWatershedSolver::formColoringJacobianFunctionWatershed");
    EBAMRWatershedSolver* solver = static_cast<EBAMRWatershedSolver*>(ctx);
    
    std::vector<PetscInt>  rows;
    std::vector<std::vector<PetscInt> > collumns;
    std::vector<PetscInt>  num_diag, num_offdiag;
    rows.clear();
    collumns.clear();
    num_diag.clear();
    num_offdiag.clear();

    solver->m_richardsCompGrid->getJacobianColoringVectors(rows, collumns, num_diag, num_offdiag, solver->m_include_SUBsurface);

    PetscInt n;
    VecGetLocalSize(f, &n);
    CH_assert(n==rows.size());
    
    for (int i=0;i<rows.size();i++)
    {
	const PetscInt num_col =collumns[i].size();
	{
	    //insert stencil terms
	    std::vector<PetscScalar> dummyVals(num_col, 1.0);
	    MatSetValues(J, 1, &(rows[i]), num_col, &(collumns[i][0]),  &(dummyVals[0]), INSERT_VALUES);
	}
    }
    MatAssemblyBegin(J, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(J, MAT_FINAL_ASSEMBLY);

  PetscFunctionReturn(0);
}


PetscErrorCode
EBAMRWatershedSolver::MatVecMult(Mat A, Vec x, Vec y)
{
    CH_TIME("EBAMRWatershedSolver::MatVecMult_Watershed");
    
    void* ctx;
    MatShellGetContext(A, &ctx);
    
    //  SNESShellGetContext(A, &ctx);
    
    //  EBAMRWatershedSolver *solver = static_cast<EBAMRWatershedSolver*>(ctx);



    PetscFunctionReturn(0);
}


PetscErrorCode 
EBAMRWatershedSolver::PCApply(PC pc, Vec x, Vec y)
{
  //currently using none
  PetscErrorCode ierr;
  ierr = VecCopy(x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

bool
EBAMRWatershedSolver::solve(Vector<Vector<LevelData<EBCellFAB>*>*>& a_SUBsurfaceComposite,
			    const Vector<RealVect>& a_vectDx,
			     const Real a_dt,
			    int a_maxLev /* =-1 */, 
			    int a_ibase /* =0 */ )
{
//    PetscErrorCode ierr;
    m_include_SUBsurface = true;

    Vector<LevelData<EBCellFAB>*>& psi_subsurface   = *a_SUBsurfaceComposite[0];
    Vector<LevelData<EBCellFAB>*>& rhs_subsurface   = *a_SUBsurfaceComposite[1];
    
    //saves reference to h_current for subsurface vector evaluation 
    //  m_richardsCompGrid->setCurrentPressure(*a_SUBsurfaceComposite[2]);
    
 
    CH_TIME("EBPetscAMRSolver::solve");
    
    m_vectDx = a_vectDx;
    m_dt = a_dt;
    
    //this will pass current ebcellfab as temporary vars  
    //this is fine since we convert rhs to petsc Vec before using them.
    m_subsur_mfree  = psi_subsurface;
    m_Fsubsur_mfree = rhs_subsurface;


    if (!m_isInitialized)
    {
	m_richardsCompGrid->setChomboPetscGlobalMap(m_petsc_x, m_petsc_mapping);
    }

    //put data into global petsc vector
    m_richardsCompGrid->putChomboInPetsc(psi_subsurface, m_petsc_x);


    VecDuplicate(m_petsc_x, &m_petsc_b);
    m_richardsCompGrid->putChomboInPetsc(rhs_subsurface,m_petsc_b);
    
    // Vec testvec = PETSC_NULL;
    // VecDuplicate(m_petsc_x, &testvec);
    // m_richardsCompGrid->fillRandom(psi_subsurface);
    // m_richardsCompGrid->putChomboInPetsc(psi_subsurface,testvec);
    // m_richardsCompGrid->putPetscInChombo(testvec,rhs_subsurface);
    // EBAMRDataOps::incr(rhs_subsurface,psi_subsurface,-1.0);
    // EBAMRDataOps::setCoveredVal(rhs_subsurface, 0.0);
    // EBAMRDataOps::exchangeAll(rhs_subsurface);
    // Real maxval=-1.0e123, minval=1.e123;
    // for (int ilev=0;ilev<rhs_subsurface.size();ilev++)
    // {
    // 	EBLevelDataOps::getMaxMin(maxval, minval,*rhs_subsurface[ilev], 0);
    // 	pout()<<"level="<<ilev<<"\t Test maxmin "<<maxval<<"\t"<<minval<<endl;
    // }
    // VecDestroy(&testvec);
 
    //VecView(m_petsc_x, 	PETSC_VIEWER_STDOUT_WORLD);
    
    bool converged = solveSystem();

    m_richardsCompGrid->putPetscInChombo(m_petsc_x, psi_subsurface);

    EBAMRDataOps::exchangeAll(psi_subsurface);
//    m_richardsCompGrid->averageDown(psi_subsurface);
    m_richardsCompGrid->averageDownHydrolic(psi_subsurface);

    VecDestroy(&m_petsc_b);
    m_petsc_b = PETSC_NULL;
    
    //deallocateSolver();

    return converged;
}

bool
EBAMRWatershedSolver::solveSurfaceFlow(Vector<Vector<LevelData<BaseIVFAB<Real> >*>*>& a_surfaceComposite,
				       const Vector<RealVect>& a_vectDx,
				       const Real a_dt,
				       int a_maxLev /* =-1 */, 
				       int a_ibase /* =0 */ )
{
//    PetscErrorCode ierr;
    m_include_SUBsurface = false;
    Vector<LevelData<BaseIVFAB<Real> >*>& psiIrr   = *a_surfaceComposite[0];
    Vector<LevelData<BaseIVFAB<Real> >*>& rhsIrr   = *a_surfaceComposite[1];
    
    CH_TIME("EBPetscAMRSolver::solve");
    
    m_vectDx = a_vectDx;
    m_dt = a_dt;
    
    //this will pass current ebcellfab as temporary vars  
    //this is fine since we convert rhs to petsc Vec before using them.
    m_sur_mfree  = psiIrr;
    m_Fsur_mfree = rhsIrr;


    if (!m_isInitialized)
    {
	m_richardsCompGrid->setChomboPetscGlobalMap(m_petsc_x, m_petsc_mapping);
    }

    //put data into global petsc vector
    m_richardsCompGrid->putChomboInPetsc(psiIrr, m_petsc_x);


    VecDuplicate(m_petsc_x, &m_petsc_b);
    m_richardsCompGrid->putChomboInPetsc(rhsIrr,m_petsc_b);
    
    // Vec testvec = PETSC_NULL;
    // VecDuplicate(m_petsc_x, &testvec);
    // m_richardsCompGrid->fillRandom(h_subsurface);
    // m_richardsCompGrid->putChomboInPetsc(h_subsurface,testvec);
    // m_richardsCompGrid->putPetscInChombo(testvec,rhs_subsurface);
    // EBAMRDataOps::incr(h_subsurface, rhs_subsurface,-1.0);
    // EBAMRDataOps::exchangeAll(h_subsurface);
    // EBAMRDataOps::setCoveredVal(h_subsurface, 0.0);
    // Real maxval=0, minval=0;
    // EBAMRDataOps::getMaxMin(maxval, minval,h_subsurface, 0);
    // pout()<<"maxmin"<<maxval<<"\t"<<minval<<endl;
    // VecDestroy(&testvec);
    bool converged = solveSystem();

    m_richardsCompGrid->putPetscInChombo(m_petsc_x, psiIrr);

    for (int ilev=0;ilev<psiIrr.size();ilev++)
    {
	psiIrr[ilev]->exchange();
    }

    VecDestroy(&m_petsc_b);
    m_petsc_b = PETSC_NULL;
    
    return converged;
}



bool 
EBAMRWatershedSolver::initializeSolver()
{
    //CH_TIMERS("EBAMRWatershedSolver::initilizeSolver()");
    
    //create operators, mat, options
    initializeSNES();  
    
    m_isInitialized=true;

    return true;
}

PetscErrorCode
EBAMRWatershedSolver::initializeSNES()
{
  PetscErrorCode ierr;

  // Create the SNES solver.
  ierr = SNESCreate( m_petsc_comm, &m_petsc_snes ); CHKERRQ(ierr);
  resetSNESOperators();
  resetSNESOptions();
  resetSNESPC();
  
  // Set the SNES options from the PETSc options database.
  if (m_options_prefix != "")
  {
      ierr = SNESSetOptionsPrefix(m_petsc_snes, m_options_prefix.c_str());
      CHKERRQ(ierr);
  }
  
  //overrides inputfile options if provided runtime
  ierr = SNESSetFromOptions(m_petsc_snes); CHKERRQ(ierr);
  
  //check PC for matrix free method
  if (m_isMatrixFree)
  {
      static const size_t len = 255;
      char pc_type_str[len];
      PetscBool flg;

#if PETSC_VERSION_GE(3,7,0)
      ierr = PetscOptionsGetString(PETSC_NULL, m_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg); CHKERRQ(ierr);
#else
      ierr = PetscOptionsGetString(m_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg); CHKERRQ(ierr);
#endif
      
      const std::string pctype(pc_type_str);
      if (!(pctype == "none" || pctype == "shell")) 
      {
	  pout()<<"EBAMRWatershedSolver::initializeSNES() WARNING!: inconsistent PC type for MatrixFree method. Setting pc_type to default none.\n";
	  PC petsc_pc;
	  ierr = KSPGetPC(m_petsc_ksp, &petsc_pc); CHKERRQ(ierr);
	  ierr = PCSetType(petsc_pc, PCNONE); CHKERRQ(ierr);
      }
  }
 
  PetscFunctionReturn(0);
} // initializeSNES

PetscErrorCode
EBAMRWatershedSolver::resetSNESOptions()
{
  PetscErrorCode ierr;
    if (!m_petsc_snes)   PetscFunctionReturn(0);
    switch (m_snes_typeID)
    {
      case 0:  m_snes_type = "shell";break;
      case 1:  m_snes_type = "newtonls";break;
      case 2:  m_snes_type = "newtontr";break;
      case 3:  m_snes_type = "nrichardson";break;
      case 4:  m_snes_type = "ncg";break;
      case 5:  m_snes_type = "qn";break;
      case 6:  m_snes_type = "fas";break;
      case 7:  m_snes_type = "nasm";break;
      case 8:  m_snes_type = "composite";break;
      case 9:  m_snes_type = "ngs";break;
      case 10: m_snes_type = "anderson";break;
      case 11: m_snes_type = "ngmres";break;
      case 12: m_snes_type = "test";break;
      default: m_snes_type = "ngmres";
    }
    pout()<< "EBAMRWatershedSolver::resetSNESOptions() SNESType is "<<m_snes_type<<endl;
    if (m_snes_typeID)
    {
	const SNESType snes_type = m_snes_type.c_str();
	ierr = SNESSetType(m_petsc_snes, snes_type); CHKERRQ(ierr);
    }

  PetscFunctionReturn(0);
} // resetKSPOptions
 
PetscErrorCode
EBAMRWatershedSolver::resetSNESOperators()
{
  pout()<< "EBAMRWatershedSolver::resetSNESOperators() MatrixFree="<<m_isMatrixFree<<"    JacobianFD="<<m_isFD<<endl;
  PetscErrorCode ierr;
  // Create and configure the MatShell object.
  if (m_petsc_mat)
    {
      MatDestroy(&m_petsc_mat);
      m_petsc_mat = NULL;
    }
  PetscInt n;
  VecGetLocalSize(m_petsc_b, &n);

  //set F vector function
  ierr = SNESSetFunction( m_petsc_snes, m_petsc_r, (PetscErrorCode (*)(SNES,Vec,Vec,void*)) &EBAMRWatershedSolver::formVecFunctionWatershed, (void*)this); CHKERRQ(ierr);

  //check if matrix free methos
  if (m_isMatrixFree)
    {
      ierr = MatCreateSNESMF(m_petsc_snes, &m_petsc_mat); CHKERRQ(ierr);
      if(m_isFD) //use finite difference to approx Jacobian
	{
	  //Jacobian based on finite differencing
	  ierr = SNESSetJacobian(m_petsc_snes, m_petsc_mat, m_petsc_mat, MatMFFDComputeJacobian, (void*) this); CHKERRQ(ierr);
	  
	  // two options MATMFFD_WP and MATMFFD_DS
	  if (m_matmffd_type==1)
	    {
	  ierr = MatMFFDSetType(m_petsc_mat, MATMFFD_DS); 
	  MatMFFDDSSetUmin(m_petsc_mat, m_mffd_ds_umin); 
	    }
	  else if (m_matmffd_type==2)
	    {
	  ierr = MatMFFDSetType(m_petsc_mat, MATMFFD_WP); CHKERRQ(ierr);
	  ierr = MatMFFDWPSetComputeNormU(m_petsc_mat, (PetscBool) m_mffd_wp_normu); CHKERRQ(ierr);
	    }
	  else
	    MayDay::Error("EBAMRWatershedSolver::resetSNESOperators() nonsupported MFFD type");
	}
      else //use shell functions for Jacobian
	{
	  //ierr = MatCreateSNESMF(m_petsc_snes, &m_petsc_mat); 
	  //shell Jacobian
	  ierr = MatCreateShell(m_petsc_comm, n, n, PETSC_DETERMINE, PETSC_DETERMINE, (void*)(this), &m_petsc_mat); CHKERRQ(ierr);
	  ierr = MatShellSetOperation(m_petsc_mat, MATOP_MULT,(void(*)())&EBAMRWatershedSolver::MatVecMult); // reinterpret_cast<void (*)(void)>(EBAMRWatershedSolver::MatVecMult)); 
	  CHKERRQ(ierr);
	  ierr = SNESSetJacobian( m_petsc_snes, m_petsc_mat, m_petsc_mat, (PetscErrorCode(*) (SNES snes,Vec ,Mat ,Mat ,void* )) &EBAMRWatershedSolver::formJacobianFunctionWatershed, (void*)this ); CHKERRQ(ierr);
	}
    }
  else //otherwise create and allocate matrix
  {
      //pout()<<"EBAMRWatershedSolver::resetSNESOperators() creating matrix"<<endl;
      // ierr = MatCreate(m_petsc_comm, &m_petsc_mat); CHKERRQ(ierr);
      // ierr = MatSetType(m_petsc_mat, MATMPIAIJ); CHKERRQ(ierr);
      // ierr = MatSetSizes(m_petsc_mat, n, n, PETSC_DETERMINE, PETSC_DETERMINE); CHKERRQ(ierr); 
      // ierr = MatSetFromOptions(m_petsc_mat); CHKERRQ(ierr);
      //ierr = MatSetUp(m_petsc_mat); CHKERRQ(ierr);
      ierr = MatCreateAIJ(m_petsc_comm, n, n, PETSC_DETERMINE, PETSC_DETERMINE, (n<36?n:36), NULL, (n<36?n:36), NULL, &m_petsc_mat); CHKERRQ(ierr); 

      ierr =MatSetLocalToGlobalMapping(m_petsc_mat, m_petsc_mapping, m_petsc_mapping); CHKERRQ(ierr);

      if(m_isFD) //use finite difference to approx Jacobian
      {

	  //Jacobian based on finite differencing
	  if (m_fd_coloring)
	  {	    
	      ISColoring iscoloring;
	      MatFDColoring fdcoloring;
	      MatColoring coloring;

	      //pout()<<"EBAMRWatershedSolver::resetSNESOperators() setup coloring of Jacobian"<<endl;
	      EBAMRWatershedSolver::formColoringJacobian(m_petsc_snes,m_petsc_x, m_petsc_mat, m_petsc_mat,(void*) this); 
 
	      ierr = MatColoringCreate(m_petsc_mat,&coloring); CHKERRQ(ierr);
	      ierr = MatColoringSetType(coloring,MATCOLORINGSL); CHKERRQ(ierr);
	      ierr = MatColoringSetFromOptions(coloring); CHKERRQ(ierr);
	      ierr = MatColoringApply(coloring,&iscoloring); CHKERRQ(ierr);
	      ierr = MatColoringDestroy(&coloring); CHKERRQ(ierr);
	      ierr = MatFDColoringCreate(m_petsc_mat,iscoloring,&fdcoloring); CHKERRQ(ierr);
	      ierr = MatFDColoringSetFunction(fdcoloring,(PetscErrorCode (*)(void)) &EBAMRWatershedSolver::formVecFunctionWatershed,(void*) this); CHKERRQ(ierr);
	      ierr = MatFDColoringSetFromOptions(fdcoloring); CHKERRQ(ierr);
	      ierr = MatFDColoringSetUp(m_petsc_mat,iscoloring,fdcoloring); CHKERRQ(ierr);
	      ierr = ISColoringDestroy(&iscoloring); CHKERRQ(ierr);

	      //pout()<<"EBAMRWatershedSolver::resetSNESOperators() set Jacobian"<<endl;
	      ierr = SNESSetJacobian(m_petsc_snes,m_petsc_mat,m_petsc_mat,SNESComputeJacobianDefaultColor,fdcoloring); CHKERRQ(ierr);
	  }
	  else 
	  {
	      // full jacobian very slow
	      ierr = SNESSetJacobian(m_petsc_snes, m_petsc_mat, NULL, SNESComputeJacobianDefault, (void*) this); CHKERRQ(ierr);
	  }	  	  
      }
      else
      {
	  ierr = SNESSetJacobian( m_petsc_snes, m_petsc_mat, m_petsc_mat, (PetscErrorCode(*) (SNES snes,Vec ,Mat ,Mat ,void* )) &EBAMRWatershedSolver::formJacobianFunctionWatershed, (void*)this ); CHKERRQ(ierr);
	  
	  // MayDay::Error("EBAMRWatershedSolver::resetSNESOperators() shell Jacobian is not supported yet ");
      }
  }
  
  
  pout()<<"EBAMRWatershedSolver::resetSNESOperators() setup KSP"<<endl;

  // Reset the configuration of the PETSc KSP object.
  ierr = SNESGetKSP( m_petsc_snes, &m_petsc_ksp ); CHKERRQ(ierr);
  ierr = KSPSetOperators(m_petsc_ksp, m_petsc_mat, m_petsc_mat); CHKERRQ(ierr);
  ierr = KSPSetReusePreconditioner(m_petsc_ksp, PETSC_TRUE); CHKERRQ(ierr);
  // ierr = KSPSetType(m_petsc_ksp, KSPFGMRES);

  //remove const nullspace 
  if (m_remove_const_nullspace)
    {
      MatNullSpace nullsp;
      ierr = MatNullSpaceCreate(m_petsc_comm, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
      CHKERRQ(ierr);
      ierr = MatSetNullSpace(m_petsc_mat, nullsp); 
      CHKERRQ(ierr);
      ierr = MatSetTransposeNullSpace(m_petsc_mat, nullsp);
      CHKERRQ(ierr);
      ierr = MatNullSpaceDestroy(&nullsp);
      CHKERRQ(ierr);
    }

  //ierr = KSPSetFromOptions(m_petsc_ksp); CHKERRQ(ierr);

  PetscFunctionReturn(0);
} // resetSNESOperators

PetscErrorCode
EBAMRWatershedSolver::resetSNESPC()
{
    if (!m_petsc_ksp)   PetscFunctionReturn(0);

    /*
    // Determine the preconditioner type to use.
    static const size_t len = 255;
    char pc_type_str[len];
    PetscBool flg;
    PetscOptionsGetString(m_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);
    std::string pc_type = m_pc_type;
    if (flg)
      {
        pc_type = std::string(pc_type_str);
      }
    else
      {
	//default pc_type is none
	pc_type="none";
      }

    PC petsc_pc;
    KSPGetPC(m_petsc_ksp, &petsc_pc);
    
    if (pc_type == "none")
      {
	PCSetType(petsc_pc, PCNONE);
      }
    else if (pc_type == "shell")
      {
	PCSetType(petsc_pc, PCSHELL);
	PCShellSetContext(petsc_pc, static_cast<void*>(this));
	//PCShellSetApply(petsc_pc, EBAMRWatershedSolver::PCApply_W);
      }
    */
    PetscFunctionReturn(0);
} // resetSNESPC

bool 
EBAMRWatershedSolver::solveSystem()
{
    //   pout()<<"checkpoint"<<endl;
    PetscErrorCode ierr;
    VecDuplicate(m_petsc_x, &m_petsc_r);
    
    if (!m_isInitialized) initializeSolver();
    

    if(m_presolver && m_presolver_iterations)
    {
	ierr = SNESSetOptionsPrefix(m_petsc_snes, "pre_");
	
	PetscReal atol, rtol, stol;
	PetscInt maxit, maxf;
	ierr = SNESGetTolerances(m_petsc_snes, &atol, &rtol, &stol, &maxit, &maxf);
	CHKERRQ(ierr);
	
	pout()<<"EBAMRWatershedSolver::presolver step.\n";
	
	PetscInt numIt = m_presolver_iterations;
	ierr = SNESSetTolerances(m_petsc_snes, atol, rtol, stol, numIt,maxf);
	ierr = SNESSetType(m_petsc_snes, SNESNRICHARDSON);
	CHKERRQ(ierr);
	ierr = SNESSolve(m_petsc_snes, m_petsc_b, m_petsc_x);
	//CHKERRQ(ierr);
	ierr = SNESSetTolerances(m_petsc_snes, atol, rtol, stol, maxit, maxf);
	//CHKERRQ(ierr);
	resetSNESOptions();
	ierr = SNESSetOptionsPrefix(m_petsc_snes, m_options_prefix.c_str());
	CHKERRQ(ierr);
    }
  
    // Solve the system using a PETSc KSP object.
    ierr = SNESSolve(m_petsc_snes, m_petsc_b, m_petsc_x);
    //CHKERRQ(ierr);
    ierr = SNESGetIterationNumber(m_petsc_snes, &m_current_iterations);
    //CHKERRQ(ierr);
    // SNESGetResidualNorm(m_petsc_snes, &m_current_residual_norm);
    // MatView(m_petsc_mat,PETSC_VIEWER_STDOUT_WORLD);	
    // Determine the convergence reason.
    SNESConvergedReason reason;
    SNESGetConvergedReason(m_petsc_snes, &reason);
    //   pout()<<"reason="<<reason<<endl;
    
    const bool converged = (static_cast<int>(reason) > 0);
    //if (m_enable_logging) 
    reportSNESConvergedReason(reason);
    
    // Deallocate the solver, when necessary.
    VecDestroy(&m_petsc_r);
    m_petsc_r = NULL;
    
    //if (deallocate_after_solve) deallocateSolver();
    
    return converged;
}

void
EBAMRWatershedSolver::reportSNESConvergedReason(const SNESConvergedReason& reason) const
{
    switch (static_cast<int>(reason))
    {
    case SNES_CONVERGED_FNORM_RELATIVE:
        pout() << ": converged: ||F|| < rtol --- residual norm is less than specified relative tolerance.\n";
        break;
    case SNES_CONVERGED_FNORM_ABS:
        pout() << ": converged: ||F|| <= atol --- residual norm is less than specified absolute tolerance.\n";
        break;
    case SNES_CONVERGED_SNORM_RELATIVE:
        pout() << "converged: || delta x || < stol || x || Newton computed step size small.\n";
        break;
    case SNES_CONVERGED_ITS:
        pout() << ": converged: maximum number of iterations reached.\n";
        break;
    case SNES_CONVERGED_TR_DELTA:
        pout() << ": converged: TR Delta.\n";
        break;
    case SNES_DIVERGED_FUNCTION_DOMAIN:
        pout() << ": diverged: the new x location passed the function is not in the domain of F.\n";
        break;
    case SNES_DIVERGED_FUNCTION_COUNT:
        pout() << ": diverged: function has been called more times then the final argument to SNESSetTolerances().\n";
        break;
    case SNES_DIVERGED_MAX_IT:
        pout() << ": diverged: reached maximum number of iterations before any convergence criteria were satisfied.\n";
        break;
    case SNES_DIVERGED_LINEAR_SOLVE:
        pout() << ": diverged: the linear solver failed.\n";
        break;
    case SNES_DIVERGED_FNORM_NAN:
        pout() << ":diverged: norm is nan.\n";
        break;
    case SNES_DIVERGED_LINE_SEARCH:
        pout() << ": diverged: line search is failed.\n";
        break;
    case SNES_DIVERGED_INNER:
        pout() << ": diverged: inner solve failed.\n";
        break;
    case SNES_DIVERGED_LOCAL_MIN:
        pout() << ": diverged: || J^T b || is small, implies converged to local minimum of F().\n";
        break;
    case SNES_CONVERGED_ITERATING:
        pout() << ": iterating: SNESSolve() is still running.\n";
        break;
    default:
        pout() << ": unknown completion code " << static_cast<int>(reason) << " reported.\n";
        break;
    }
    return;
} // reportKSPConvergedReason
