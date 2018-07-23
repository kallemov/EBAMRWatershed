#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _NEWPETSCSOLVER_H_
#define _NEWPETSCSOLVER_H_
#ifdef CH_USE_PETSC
#include "LinearSolver.H"
#include "petsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscviewer.h"
#include "NamespaceHeader.H"
template <class T>
class NewPetscSolver: public LinearSolver<LevelData<EBCellFAB> > 
{
public:
    NewPetscSolver();
    
    ~NewPetscSolver();
    
    virtual void setHomogeneous(bool a_homogeneous)
    {
	m_homogeneous = a_homogeneous;
    }

    virtual void define( LinearOp<LevelData<EBCellFAB> > *a_op,
			 bool a_homogeneous = true );
    
    virtual void define(const EBLevelGrid& a_eblg);
    
    virtual void solve( LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs );
        
    bool initializeSolver(const PetscInt a_numComp=1);
    
    void destroy();
    
    /*!
     * \brief Set the options prefix used by this PETSc solver object.
     */
    void setOptionsPrefix(const std::string options_prefix);
    
    /*!
     * \brief Compute the matrix vector product \f$y=Ax\f$.
     */
    static PetscErrorCode MatVecMult(Mat A, Vec x, Vec y);
    
    /*!
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply(PC pc, Vec x, Vec y);

    bool solveSystem();
    
    void applyOp(LevelData<EBCellFAB>& a_lhs,
		 const LevelData<EBCellFAB>& a_phi,
		 const bool a_homogeneous = true);
    
    void residual(LevelData<EBCellFAB>& a_lhs,
		  const LevelData<EBCellFAB>& a_phi,
		  const LevelData<EBCellFAB>& a_rhs,
		  const bool a_homogeneous = true);
	
    
    
private:
    bool m_homogeneous;
    
    bool m_isMatrixFree;

    /*!
     * \brief Report the SNESConvergedReason.
     */
    void reportKSPConvergedReason(const KSPConvergedReason& reason) const;
    
    /*!
     * \brief Routine to setup SNES object.
     */
    PetscErrorCode initializeKSP();
    
    /*!
     * \brief Reset the values of the convergence tolerances for the PETSc SNES
     * object.
     */
    PetscErrorCode resetKSPOptions();
    
    /*!
     * \brief Reset the SNES operators to correspond to the supplied
     * LinearOperator.
     */
    PetscErrorCode resetKSPOperators();
    
    /*!
     * \brief Reset the SNES PC to correspond to the supplied preconditioner.
     */
    PetscErrorCode resetKSPPC();
    

    void allocateTemporaries(const PetscInt a_nComp);
    
    void deallocateTemporaries();
    
    // Solver stuff
    bool m_isInitialized;
    int m_verbosity;
    Mat m_petsc_mat;//, m_petsc_precond_mat;
    
    Vec m_petsc_x, m_petsc_b; 
    KSP m_petsc_ksp; // linear solver context
    MPI_Comm m_petsc_comm;
    std::string m_options_prefix;
    std::string m_ksp_type, m_pc_type;
    int m_max_iterations;
    PetscInt m_current_iterations;
    double m_abs_residual_tol, m_rel_residual_tol;
    double m_current_residual_norm;
    bool m_enable_logging;
    bool m_isDefined;
    bool m_remove_const_nullspace;
    
    EBLevelGrid m_eblg;
    //EBAMRPoissonOp* m_op;
    T* m_petscGridMapping;
    //mapping Local to Global
    ISLocalToGlobalMapping m_petsc_mapping;
    
    //temporary data holders for matrix vector operations
    Vector<LevelData<EBCellFAB>*> m_phi_mfree;
    Vector<LevelData<EBCellFAB>*> m_rhs_mfree;
    
    //temporaries for MatSetSizes
    PetscInt *m_d_nnz, *m_o_nnz;
};
#include "NamespaceFooter.H"
#include "NewPetscSolverI.H"
#endif
#endif
