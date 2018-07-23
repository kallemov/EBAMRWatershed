#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBAMRWATERSHEDSOLVER_H_
#define _EBAMRWATERSHEDSOLVER_H_

#include <petscsnes.h>
#include "WatershedIBC.H"
#include "EBPetscCompGridRichards.H"
#include "LevelData.H"
#include "EBCellFAB.H"

#include "UsingNamespace.H"

class EBAMRWatershedSolver
{
public:
    
    /*!
     * \brief constructor of nonlinear solver class 
     * for coupled surface water and subsurface solver 
     *
     */
    EBAMRWatershedSolver(EBPetscCompGridRichards *a_richardsCompGrid, WatershedIBC *a_ibc);
    
    /*!
     * \brief destructor
     *
     */
    ~EBAMRWatershedSolver();
    
    /*!
     * \brief Remove all hierarchy dependent data allocated by
     * initializeSolverState().
     *
     * \note It is safe to call deallocateSolverState() when the solver state is
     * already deallocated.
     */
    void deallocateSolver();
    
    
    bool solve(Vector<Vector<LevelData<EBCellFAB>*>*>& a_subSurfaceComposite,
	       const Vector<RealVect>& a_vectDx,
	       const Real a_dt,
	       int a_maxLev  =-1, int a_ibase  =0);

    bool solveSurfaceFlow(Vector<Vector<LevelData<BaseIVFAB<Real> >*>*>& a_surfaceComposite,
			  const Vector<RealVect>& a_vectDx,
			  const Real a_dt,
			  int a_maxLev =-1, int a_ibase  =0);
    
    static PetscErrorCode formVecFunctionWatershed(SNES snes,Vec x,Vec f,void * foo);
    static PetscErrorCode formJacobianFunctionWatershed(SNES snes,Vec f,Mat J,Mat P,void* ctx);
    static PetscErrorCode formColoringJacobian(SNES snes,Vec f,Mat J,Mat P,void* ctx);
    
    /*!
     * \brief Set the SNES type.
     */
    void setSNESType(const std::string& ksp_type);
    
    /*!
     * \brief Set the options prefix used by this PETSc solver object.
     */
    void setOptionsPrefix(const std::string& options_prefix);
    
    /*!
     * \brief Get the PETSc KSP object.
     */
    const KSP& getPETScKSP() const;
    
    
    /*!
     * \brief Compute the matrix vector product \f$y=Ax\f$.
     */
    static PetscErrorCode MatVecMult(Mat A, Vec x, Vec y);
    
    /*!
     * \brief Apply the preconditioner to \a x and store the result in \a y.
     */
    static PetscErrorCode PCApply(PC pc, Vec x, Vec y);
    
    
    /*!
     * \brief Compute hierarchy dependent data 
     *
     * \param x solution vector
     */
    bool initializeSolver();
    
    /*!
     * \brief Solve the non-linear system of equations
     *
     * \param x solution multivector
     * \param b right-hand-side multivector
     * \return \p true if the solver converged to the specified tolerances, \p
     * false otherwise
     */
    bool solveSystem();
    
    
private:
    EBPetscCompGridRichards *m_richardsCompGrid;
    WatershedIBC *m_ibc;
    //temporary data holders for matrix vector operations
    Vector<LevelData<EBCellFAB>*> m_subsur_mfree;
    Vector<LevelData<EBCellFAB>*> m_Fsubsur_mfree;

    Vector<LevelData<BaseIVFAB<Real> >*> m_sur_mfree;
    Vector<LevelData<BaseIVFAB<Real> >*> m_Fsur_mfree;
    
    bool m_include_SUBsurface;
    
    bool m_isMatrixFree;
    bool m_isFD;
    bool m_fd_coloring;
    int m_matmffd_type;

    //parameters for MFFD types (see snes docs)
    Real m_mffd_ds_umin;
    bool m_mffd_wp_normu;
    
    /*!
     * \brief Report the SNESConvergedReason.
     */
    void reportSNESConvergedReason(const SNESConvergedReason& reason) const;
    
    /*!
     * \brief Routine to setup SNES object.
     */
    PetscErrorCode initializeSNES();
    
    /*!
     * \brief Reset the values of the convergence tolerances for the PETSc SNES
     * object.
     */
    PetscErrorCode resetSNESOptions();
    
    /*!
     * \brief Reset the SNES operators to correspond to the supplied
     * LinearOperator.
     */
    PetscErrorCode resetSNESOperators();
    
    /*!
     * \brief Reset the SNES PC to correspond to the supplied preconditioner.
     */
    PetscErrorCode resetSNESPC();
    
    /*!
     * \brief Set SNES monitoring routine for the SNES.
     */
    PetscErrorCode monitorSNES(SNES snes, int it, PetscReal rnorm, void* mctx);
    
    
    // Solver stuff
    bool m_isInitialized;
    Mat m_petsc_mat, m_petsc_precond_mat;
    
    Vec m_petsc_x, m_petsc_b, m_petsc_r; 
    SNES m_petsc_snes; //solver context
    KSP m_petsc_ksp; // linear solver context
    MPI_Comm m_petsc_comm;
    std::string m_options_prefix;
    std::string m_snes_type, m_ksp_type, m_pc_type;
    int m_max_iterations;
    PetscInt m_current_iterations;
    double m_abs_residual_tol, m_rel_residual_tol;
    double m_current_residual_norm;
    bool m_presolver;
    int m_presolver_iterations;
    bool m_enable_logging;
    bool m_isDefined;
    bool m_remove_const_nullspace;
    
    //grid res (can be got from CompGrid)
    //leave it as additional copy for now.
    Vector<RealVect>            m_vectDx;
    //current time step
    Real m_dt;
    
    //parameter for snes type
    int m_snes_typeID;
    
    //mapping Local to Global
    ISLocalToGlobalMapping m_petsc_mapping;
};

#endif

