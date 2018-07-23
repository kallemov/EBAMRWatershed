#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEWPETSCSOLVERI_H_
#define _NEWPETSCSOLVERI_H_

#ifdef CH_USE_PETSC
//#include "advisor-annotate.h"

#include "ParmParse.H"
#include "EBLevelDataOps.H"
#include "EBAMRDataOps.H"
#include "NewPetscSolver.H"
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "EBPetscMap.H"
#include "NamespaceHeader.H"

template <class T>
NewPetscSolver<T>::NewPetscSolver()
    :m_homogeneous(false),
     m_isMatrixFree(false),
     m_isInitialized(false),
     m_verbosity(1),
     m_petsc_mat(NULL), 
     m_petsc_x(PETSC_NULL), 
     m_petsc_b(PETSC_NULL), 
     m_petsc_ksp(PETSC_NULL),
     m_isDefined(false),
     m_remove_const_nullspace(false),
     //m_op(NULL),
     m_petscGridMapping(NULL),
     m_d_nnz(PETSC_NULL),
     m_o_nnz(PETSC_NULL)
{}


template <class T>
void
NewPetscSolver<T>::define( LinearOp<LevelData<EBCellFAB> > *a_op,
	      bool a_homogeneous)
{
    //destroy();
    if (!m_isDefined)
    {
	if (!m_petscGridMapping) m_petscGridMapping = new T();
	m_petscGridMapping->setOp(a_op);
	m_homogeneous = a_homogeneous;
	define(m_petscGridMapping->m_op->getEBLG());
    }
}

template <class T>
void 
NewPetscSolver<T>::define(const EBLevelGrid& a_eblg)
{
#ifdef CH_MPI
    m_petsc_comm = Chombo_MPI::comm;
#else
    m_petsc_comm = PETSC_COMM_WORLD;
#endif


    m_eblg = a_eblg;
    m_petscGridMapping->defineGrids(0, Vector<const EBLevelGrid*>(1, &m_eblg), Vector<int>(), 
				    true,
				    true,
				    false);
    
    ParmParse pp;
    
    m_isMatrixFree = false;
    pp.query("MatrixFree", m_isMatrixFree);
    
    string PetscOptionsFile;
    pp.query("petsc_options_file",PetscOptionsFile);
#if PETSC_VERSION_GE(3,7,0)
    if (PetscOptionsFile!="")
	PetscOptionsInsertFile(m_petsc_comm, PETSC_NULL,PetscOptionsFile.c_str(), PETSC_TRUE);
#else
    if (PetscOptionsFile!="")
	PetscOptionsInsertFile(m_petsc_comm, PetscOptionsFile.c_str(), PETSC_TRUE);
#endif
    m_isDefined = true;
}

template <class T>
NewPetscSolver<T>::~NewPetscSolver()
{
    destroy();
}

template <class T>
void NewPetscSolver<T>::destroy()
{
    if (m_petscGridMapping) delete m_petscGridMapping;
    m_petscGridMapping=NULL;
    m_isDefined = false;
    
    if (!m_isInitialized) return;
    m_isInitialized=false;
    
    if (m_petsc_mat) MatDestroy(&m_petsc_mat);
    m_petsc_mat = PETSC_NULL;
    
    if ( m_petsc_ksp ) KSPDestroy(&m_petsc_ksp);
    m_petsc_ksp = PETSC_NULL;
    
    if ( m_petsc_x ) VecDestroy(&m_petsc_x);
    m_petsc_x = PETSC_NULL;
    if ( m_petsc_b ) VecDestroy(&m_petsc_b);
    m_petsc_b = PETSC_NULL;
}

template <class T>
PetscErrorCode 
NewPetscSolver<T>::MatVecMult(Mat A, Vec x, Vec y)
{
    CH_TIME("NewPetscSolver::MatVecMult");
    
    void* ctx;
    MatShellGetContext(A, &ctx);
    NewPetscSolver *solver = static_cast<NewPetscSolver*>(ctx);

    solver->m_petscGridMapping->putPetscInChombo(x, solver->m_phi_mfree);

    solver->m_petscGridMapping->m_op->applyOp(*(solver->m_rhs_mfree[0]), *(solver->m_phi_mfree[0]), solver->m_homogeneous);
    
    solver-> m_petscGridMapping->putChomboInPetsc(solver->m_rhs_mfree, y);
    
    PetscFunctionReturn(0);
}

template <class T>
PetscErrorCode NewPetscSolver<T>::PCApply(PC pc, Vec x, Vec y)
{
  //currently using none
  PetscErrorCode ierr;
  ierr = VecCopy(x,y);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

template <class T>
bool NewPetscSolver<T>::initializeSolver(const PetscInt a_nComp)
{
    CH_TIMERS("NewPetscSolver::initilizeSolver()");
    CH_TIMER("getMatSizes", t1);
    CH_TIMER("initializeKSP", t2);
    CH_START(t1);
    allocateTemporaries(a_nComp);
    
    m_petscGridMapping->getMatAIJSizes(m_d_nnz, m_o_nnz, a_nComp);
    CH_STOP(t1);

    CH_START(t2);
    //create operators, mat, options
    initializeKSP();  
    CH_STOP(t2);
      
    deallocateTemporaries();
    m_isInitialized=true;

    return true;
}

template <class T>
PetscErrorCode NewPetscSolver<T>::initializeKSP()
{
  PetscErrorCode ierr;

  ierr = KSPCreate(m_petsc_comm, &m_petsc_ksp); CHKERRQ(ierr);

  resetKSPOperators();
  //resetKSPOptions();
  resetKSPPC();
  
  if (m_options_prefix != "")
  {
      ierr = KSPSetOptionsPrefix(m_petsc_ksp, m_options_prefix.c_str());
      CHKERRQ(ierr);
  }
  
  //overrides inputfile options if provided runtime
  ierr = KSPSetFromOptions(m_petsc_ksp); CHKERRQ(ierr);
  
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
	  if (m_verbosity>0)
	      pout()<<"NewPetscSolver::initializeKSP() WARNING!: inconsistent PC type for MatrixFree method. Setting pc_type to default none.\n";
	  PC petsc_pc;
	  ierr = KSPGetPC(m_petsc_ksp, &petsc_pc); CHKERRQ(ierr);
	  ierr = PCSetType(petsc_pc, PCNONE); CHKERRQ(ierr);
      }
  }else
  {
      
  }
 
  PetscFunctionReturn(0);
} // initializeKSP

template <class T>
PetscErrorCode NewPetscSolver<T>::resetKSPOperators()
{
    if (m_verbosity>2)
	pout()<< "NewPetscSolver::resetKSPOperators() MatrixFree="<<m_isMatrixFree<<endl;
  PetscErrorCode ierr;
  // Create and configure the MatShell object.
  if (m_petsc_mat)
    {
      MatDestroy(&m_petsc_mat);
      m_petsc_mat = NULL;
    }
  PetscInt n;
  VecGetLocalSize(m_petsc_b, &n);

  if (m_isMatrixFree)
  {
      ierr = MatCreateShell(m_petsc_comm, n, n, PETSC_DETERMINE, PETSC_DETERMINE, static_cast<void*>(this), &m_petsc_mat); CHKERRQ(ierr);
      
      ierr = MatShellSetOperation(m_petsc_mat, MATOP_MULT, (void(*)(void))MatVecMult); CHKERRQ(ierr); 
  }
  else
  {
      ierr = MatCreate(m_petsc_comm,&m_petsc_mat);CHKERRQ(ierr);
      ierr = MatSetSizes(m_petsc_mat, n, n, PETSC_DETERMINE, PETSC_DETERMINE);CHKERRQ(ierr);
      ierr = MatSetBlockSize(m_petsc_mat,1);CHKERRQ(ierr);
      ierr = MatSetType(m_petsc_mat,MATAIJ);CHKERRQ(ierr);
      ierr = MatSeqAIJSetPreallocation(m_petsc_mat,n, m_d_nnz);CHKERRQ(ierr);
      ierr = MatMPIAIJSetPreallocation(m_petsc_mat,0, m_d_nnz, 0, m_o_nnz);CHKERRQ(ierr);
      ierr = MatSetFromOptions( m_petsc_mat ); CHKERRQ(ierr);
      //ierr = MatSetUp(m_petsc_mat);

#ifdef CH_MPI
      MPI_Barrier(Chombo_MPI::comm);
#endif
      //exit(0);
      ierr = MatSetLocalToGlobalMapping(m_petsc_mat, m_petsc_mapping, m_petsc_mapping); CHKERRQ(ierr);
  }
  
  if (m_verbosity>3)
      pout()<<"NewPetscSolver::resetKSPOperators() setup KSP"<<endl;

  // Reset the configuration of the PETSc KSP object.
  ierr = KSPSetOperators(m_petsc_ksp, m_petsc_mat, m_petsc_mat); CHKERRQ(ierr);
  ierr = KSPSetReusePreconditioner(m_petsc_ksp, PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPSetType(m_petsc_ksp, KSPGMRES);

  //remove const nullspace 
  if (m_remove_const_nullspace)
  {
      MatNullSpace nullsp;
      ierr = MatNullSpaceCreate(m_petsc_comm, PETSC_TRUE, 0, PETSC_NULL, &nullsp);
      CHKERRQ(ierr);
      ierr = MatSetNullSpace(m_petsc_mat, nullsp); 
      CHKERRQ(ierr);
      ierr = MatNullSpaceDestroy(&nullsp);
      CHKERRQ(ierr);
  }

  PetscFunctionReturn(0);
} // resetKSPOperators

template <class T>
PetscErrorCode NewPetscSolver<T>::resetKSPPC()
{
    if (!m_petsc_ksp)   PetscFunctionReturn(0);

    
    // Determine the preconditioner type to use.
    static const size_t len = 255;
    char pc_type_str[len];
    PetscBool flg;
#if PETSC_VERSION_GE(3,7,0)
    PetscOptionsGetString(PETSC_NULL, m_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);
#else
    PetscOptionsGetString(m_options_prefix.c_str(), "-pc_type", pc_type_str, len, &flg);
#endif
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
	//PCShellSetApply(petsc_pc, NewPetscSolver::PCApply);
    }
    
    PetscFunctionReturn(0);
} // resetKSPPC

template <class T>
bool NewPetscSolver<T>::solveSystem()
{

//ANNOTATE_DISABLE_COLLECTION_PUSH;
    int ierr;
    PetscInt numIterations;
    Real residual;
    
    // const int num = numProc();
    // char iter_str[80];
    // sprintf(iter_str, "A%03d.m",num);
    // PetscViewer viewer;
    // PetscViewerASCIIOpen( PETSC_COMM_WORLD, iter_str, &viewer);
    // PetscViewerSetFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
    // //MatView( a_petsc_mat, viewer );
    // VecView( m_petsc_x, viewer );
    // PetscViewerDestroy(&viewer);

    
    ierr = KSPSolve(m_petsc_ksp, m_petsc_b, m_petsc_x); CHKERRQ(ierr);
    ierr = KSPGetIterationNumber(m_petsc_ksp, &numIterations); CHKERRQ(ierr);
    ierr = KSPGetResidualNorm(m_petsc_ksp, &residual); CHKERRQ(ierr);
 	// Determine the convergence reason.
	KSPConvergedReason reason;
	ierr = KSPGetConvergedReason(m_petsc_ksp, &reason); CHKERRQ(ierr);
	const bool converged = (static_cast<int>(reason) > 0);
   if (m_verbosity>3)
    {
	reportKSPConvergedReason(reason);
    }
//ANNOTATE_DISABLE_COLLECTION_POP;
    return converged;
}// solveSystem


template <class T>
void NewPetscSolver<T>::reportKSPConvergedReason(const KSPConvergedReason& reason) const
{
    
    switch (static_cast<int>(reason))
    {
    case KSP_CONVERGED_RTOL:
	pout() <<": converged: |Ax-b| <= rtol*|b| --- residual norm is less than specified relative tolerance.\n";
	break;
    case KSP_CONVERGED_ATOL:
	pout() <<": converged: |Ax-b| <= atol --- residual norm is less than specified absolute tolerance.\n";
	break;
    case KSP_CONVERGED_ITS:
	pout() <<": converged: maximum number of iterations reached.\n";
	break;
    case KSP_CONVERGED_STEP_LENGTH:
	pout() <<": converged: step size less than specified tolerance.\n";
	break;
    case KSP_DIVERGED_NULL:
	pout() <<": diverged: null.\n";
	break;
    case KSP_DIVERGED_ITS:
	pout() <<": diverged: reached maximum number of iterations before any convergence criteria were satisfied.\n";
	break;
    case KSP_DIVERGED_DTOL:
	pout() <<": diverged: |Ax-b| >= dtol*|b| --- residual is greater than specified divergence tolerance.\n";
	break;
    case KSP_DIVERGED_BREAKDOWN:
	pout() <<": diverged: breakdown in the Krylov method.\n";
	break;
    case KSP_DIVERGED_BREAKDOWN_BICG:
	pout() <<": diverged: breakdown in the bi-congugate gradient method.\n";
	break;
    case KSP_DIVERGED_NONSYMMETRIC:
	pout() <<": diverged: it appears the operator or preconditioner is not symmetric, but this Krylov method (KSPCG, KSPMINRES, KSPCR) requires symmetry\n";
	break;
    case KSP_DIVERGED_INDEFINITE_PC:
	pout() <<": diverged: it appears the preconditioner is indefinite (has both positive and negative eigenvalues), but this Krylov method (KSPCG) requires it to be positive definite.\n";
	break;
    case KSP_CONVERGED_ITERATING:
	pout() <<": iterating: KSPSolve() is still running.\n";
	break;
    default:
	pout() <<": unknown completion code " << static_cast<int>(reason) << " reported.\n";
	break;
    }
    
    return;
}// reportKSPConvergedReason

template <class T>
void NewPetscSolver<T>::setOptionsPrefix(const std::string a_prefix)
{
    m_options_prefix = a_prefix;
}

template <class T>
void NewPetscSolver<T>::solve(LevelData<EBCellFAB>& a_phi,
			   const LevelData<EBCellFAB>& a_rhs)
{
    CH_TIMERS("NewPetscSolver::solve");
    CH_TIMER("formPetscVectors", t1);
    CH_TIMER("formPetscMatrix", t2);
    CH_TIMER("solvePetscSystem", t3);
    CH_TIMER("getSolution", t4);
  
    if (!m_isDefined)
    {
	MayDay::Error("NewPetscSolver:: Solver is not properly defined.");
    }
    
    CH_START(t1);
    const PetscInt nComp = a_phi.nComp();
    
    Vector<LevelData<EBCellFAB>*> phiAMR, rhsAMR;
    phiAMR.push_back(&a_phi);
    rhsAMR.push_back(const_cast<LevelData<EBCellFAB>*>(&a_rhs));
    if (m_isMatrixFree)
    {
	m_phi_mfree = phiAMR;
	m_rhs_mfree.resize(1, NULL);
	IntVect ghosts = a_rhs.ghostVect();
	const DisjointBoxLayout &dbl = a_phi.disjointBoxLayout();
	const EBCellFactory ebcellfact(m_eblg.getEBISL() );
	LevelData<EBCellFAB>* mult_rhs = new LevelData<EBCellFAB>(dbl, 1, ghosts, ebcellfact);
	m_rhs_mfree[0] = mult_rhs;
    }

    if (!m_isInitialized)
    {
	m_petscGridMapping->setChomboPetscGlobalMap(m_petsc_b, m_petsc_mapping, nComp);
    }
    //put data into global petsc vector
    m_petscGridMapping->putChomboInPetsc(rhsAMR,m_petsc_b);


    if (!m_isInitialized)
    {
	PetscErrorCode ierr;
        ierr=VecDuplicate(m_petsc_b, &m_petsc_x);CHKERRV(ierr);
    }
    
    // if (!m_homogeneous)
    {
	//we need adjust rhs for domainBC
	m_petscGridMapping->addBCToRHS(m_petsc_b, a_phi, m_homogeneous);
    }
    CH_STOP(t1);
     
    CH_START(t2);
    if (!m_isInitialized)
    {
	initializeSolver(nComp);
	if (!m_isMatrixFree)
	{
	    m_petscGridMapping->formPetscMat(m_petsc_mat, a_phi);
	}
    }
    CH_STOP(t2);

    CH_START(t3);
    solveSystem();
    CH_STOP(t3);
    
    CH_START(t4);
    m_petscGridMapping->putPetscInChombo(m_petsc_x, phiAMR);
    a_phi.exchange();
    CH_STOP(t4);
    
    if (m_isMatrixFree) delete m_rhs_mfree[0];

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
}

template <class T>
void NewPetscSolver<T>::applyOp(LevelData<EBCellFAB>& a_lhs,
			     const LevelData<EBCellFAB>& a_phi,
			     const bool a_homogeneous)
{
    if (!m_isDefined)
    {
	MayDay::Error("NewPetscSolver:: Solver is not properly defined.");
    }
    const PetscInt nComp = a_phi.nComp();
    if (m_isMatrixFree)
    {
	m_petscGridMapping->m_op->applyOp(a_lhs, a_phi, a_homogeneous);
	return;	
    }

    Vector<LevelData<EBCellFAB>*> lhsAMR, phiAMR;
    lhsAMR.push_back(&a_lhs);
    phiAMR.push_back(const_cast<LevelData<EBCellFAB>*>(&a_phi));

    if (!m_isInitialized)
    {
	m_petscGridMapping->setChomboPetscGlobalMap(m_petsc_x, m_petsc_mapping, nComp);
    }

    //put data into global petsc vector
    m_petscGridMapping->putChomboInPetsc(phiAMR, m_petsc_x);
    if (!m_isInitialized)
    {
      VecDuplicate(m_petsc_x, &m_petsc_b);
    }
    
    if (!m_isInitialized)
    {
	initializeSolver(nComp);
	if (!m_isMatrixFree)
	{
	    m_petscGridMapping->formPetscMat(m_petsc_mat, a_phi);
	}
    }
    
    MatMult( m_petsc_mat, m_petsc_x, m_petsc_b); 

    m_petscGridMapping->putPetscInChombo(m_petsc_b, lhsAMR);

    a_lhs.exchange();
    
    return;
}

template <class T>
void NewPetscSolver<T>::residual(LevelData<EBCellFAB>& a_lhs,
			      const LevelData<EBCellFAB>& a_phi,
			      const LevelData<EBCellFAB>& a_rhs,
			      const bool a_homogeneous)
{
    if (!m_isDefined)
    {
	MayDay::Error("NewPetscSolver:: Solver is not properly defined.");
    }

    if (m_isMatrixFree)
    {
	m_petscGridMapping->m_op->residual(a_lhs, a_phi, a_rhs, a_homogeneous);
	return;	
    }
    const PetscInt nComp = a_phi.nComp();
    Vector<LevelData<EBCellFAB>*> lhsAMR, phiAMR, rhsAMR;
    lhsAMR.push_back(&a_lhs);
    phiAMR.push_back(const_cast<LevelData<EBCellFAB>*>(&a_phi));
    rhsAMR.push_back(const_cast<LevelData<EBCellFAB>*>(&a_rhs));

    if (!m_isInitialized)
    {
	m_petscGridMapping->setChomboPetscGlobalMap(m_petsc_x, m_petsc_mapping);
    }

    //put data into global petsc vector
     m_petscGridMapping->putChomboInPetsc(phiAMR, m_petsc_x);

    if (!m_isInitialized)
    {
      VecDuplicate(m_petsc_x, &m_petsc_b);
    }

    Vec lhs_petsc;
    VecDuplicate(m_petsc_x, &lhs_petsc);
    
    m_petscGridMapping->putChomboInPetsc(rhsAMR,m_petsc_b);
    if (!a_homogeneous)
    {
	//we need adjust rhs for domainBC
	m_petscGridMapping->addBCToRHS(m_petsc_b, a_phi, a_homogeneous);
    }
    
    if (!m_isInitialized)
    {
	initializeSolver(nComp);
	if (!m_isMatrixFree)
	{
	    m_petscGridMapping->formPetscMat(m_petsc_mat, a_phi);
	}
    }
    
    MatResidual( m_petsc_mat, m_petsc_b, m_petsc_x, lhs_petsc);  

    m_petscGridMapping->putPetscInChombo(lhs_petsc, lhsAMR);

    a_lhs.exchange();
    
    VecDestroy(&lhs_petsc);
    
    return;
}
	
template <class T>
void  NewPetscSolver<T>::allocateTemporaries(const PetscInt a_nComp)
{
    if (!m_isMatrixFree)
    {
	const PetscInt n = m_petscGridMapping->getPetscVecBaseSize()*a_nComp;
	m_d_nnz = new PetscInt[n];
	m_o_nnz = new PetscInt[n];
    }
    return;
}

template <class T>
void  NewPetscSolver<T>::deallocateTemporaries()
{
    if (!m_isMatrixFree)
    {
	delete[] m_d_nnz;
	m_d_nnz= PETSC_NULL;
	
	delete[] m_o_nnz;
	m_o_nnz= PETSC_NULL;
    }
}
#include "NamespaceFooter.H"
#endif
#endif
