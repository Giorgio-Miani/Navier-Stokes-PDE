#ifndef NAVIER_STOKES_HPP
#define NAVIER_STOKES_HPP

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

using namespace dealii;

// Class implementing a solver for the Navier-Stokes problem.
class NavierStokes
{
public:
  // Physical dimension (1D, 2D, 3D)
  static constexpr unsigned int dim = 3;

  // Function for the forcing term.
  class ForcingTerm : public Function<dim>
  {
  public:
    virtual void
    vector_value(const Point<dim> & /*p*/,
                 Vector<double> &values) const override
    {
      for (unsigned int i = 0; i < dim; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  // Function for the initial condition.
  class FunctionU0 : public Function<dim>
  {
  public:
    virtual void
    vector_value(const Point<dim> & /*p*/,
                 Vector<double> &values) const override
    {
      values[0] = 0.0;
      values[1] = 0.0;
      values[2] = 0.0;
    }

    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int /*component*/ = 0) const override
    {
      return 0.0;
    }
  };

  // Function for inlet velocity.
  class InletVelocity : public Function<dim>
  {
  public:
    InletVelocity()
      : Function<dim>(dim + 1)
    {}

    virtual void
    vector_value(const Point<dim> &p, Vector<double> &values) const override
    {
      values[0] = 16 * Um * p[1] * p[2] *(H - p[1]) * (H - p[2]) / (std::pow(H, 4));
      for (unsigned int i = 1; i < dim + 1; ++i)
        values[i] = 0.0;
    }

    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      if (component == 0){
          return 16 * Um * p[1] * p[2] *(H - p[1]) * (H - p[2]) / (std::pow(H, 4));;
      }
      else{
        return 0.0;
      }
    }

  protected:
    double Um = 0.45;
    double H = 0.41;
  };

  // Identity preconditioner.
  class PreconditionIdentity
  {
  public:
    // Application of the preconditioner: we just copy the input vector (src)
    // into the output vector (dst).
    void
    vmult(TrilinosWrappers::MPI::BlockVector &      dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      dst = src;
    }

  protected:
  };

  // Block-diagonal preconditioner.
  class PreconditionBlockDiagonal
  {
  public:
    // Initialize the preconditioner, given the velocity stiffness matrix, the
    // pressure mass matrix.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass      = &pressure_mass_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector &      dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl                           solver_control_velocity(1000,
                                            1e-2 * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_cg_velocity(
        solver_control_velocity);
      solver_cg_velocity.solve(*velocity_stiffness,
                               dst.block(0),
                               src.block(0),
                               preconditioner_velocity);

      SolverControl                           solver_control_pressure(1000,
                                            1e-2 * src.block(1).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
        solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               src.block(1),
                               preconditioner_pressure);
    }

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;
  };

  // Block-triangular preconditioner.
  class PreconditionBlockTriangular
  {
  public:
    // Initialize the preconditioner, given the velocity stiffness matrix, the
    // pressure mass matrix.
    void
    initialize(const TrilinosWrappers::SparseMatrix &velocity_stiffness_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_,
               const TrilinosWrappers::SparseMatrix &B_)
    {
      velocity_stiffness = &velocity_stiffness_;
      pressure_mass      = &pressure_mass_;
      B                  = &B_;

      preconditioner_velocity.initialize(velocity_stiffness_);
      preconditioner_pressure.initialize(pressure_mass_);
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector &      dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {
      SolverControl                           solver_control_velocity(1000,
                                            1e-2 * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_cg_velocity(
        solver_control_velocity);
      solver_cg_velocity.solve(*velocity_stiffness,
                               dst.block(0),
                               src.block(0),
                               preconditioner_velocity);

      tmp.reinit(src.block(1));
      B->vmult(tmp, dst.block(0));
      tmp.sadd(-1.0, src.block(1));

      SolverControl                           solver_control_pressure(1000,
                                            1e-2 * src.block(1).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_cg_pressure(
        solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               dst.block(1),
                               tmp,
                               preconditioner_pressure);
    }

  protected:
    // Velocity stiffness matrix.
    const TrilinosWrappers::SparseMatrix *velocity_stiffness;

    // Preconditioner used for the velocity block.
    TrilinosWrappers::PreconditionILU preconditioner_velocity;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;

    // B matrix.
    const TrilinosWrappers::SparseMatrix *B;

    // Temporary vector.
    mutable TrilinosWrappers::MPI::Vector tmp;
  };


  // PCD preconditioner.
  class PreconditionBlockPCD
  {
  public:
    // Initialize the preconditioner, given the F matrix, the
    // pressure mass matrix, the B matrix, the Ap matrix and
    // the Fp matrix.
    void
    initialize(const TrilinosWrappers::SparseMatrix &F_,
               const TrilinosWrappers::SparseMatrix &pressure_mass_,
               const TrilinosWrappers::SparseMatrix &B_,
               const TrilinosWrappers::SparseMatrix &Bt_,
               const TrilinosWrappers::SparseMatrix &Ap_,
               const TrilinosWrappers::SparseMatrix &Fp_)
    {
      F                  = &F_;
      pressure_mass      = &pressure_mass_;
      B                  = &B_;
      Bt                  = &Bt_;
      Ap                 = &Ap_;
      Fp                 = &Fp_;

      preconditioner_F.initialize(F_);
      preconditioner_pressure.initialize(pressure_mass_);
      preconditioner_Ap.initialize(Ap_);  
    }

    // Application of the preconditioner.
    void
    vmult(TrilinosWrappers::MPI::BlockVector &      dst,
          const TrilinosWrappers::MPI::BlockVector &src) const
    {

      // block 0

      tmp1.reinit(src.block(0));
      tmp2.reinit(src.block(0));

      SolverControl solver_control_pressure(1000, 1e-2 * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_cg_pressure(solver_control_pressure);
      solver_cg_pressure.solve(*pressure_mass,
                               tmp1,
                               src.block(1),
                               preconditioner_pressure);
  
      Fp->vmult(tmp2, tmp1);

      SolverControl solver_control_Ap(1000, 1e-2 * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_cg_Ap(solver_control_Ap);
      solver_cg_Ap.solve(*Ap,
                          tmp1,
                          tmp2,
                          preconditioner_Ap);


      Bt->vmult(tmp2, tmp1);

      tmp2.sadd(1.0, src.block(0));

      SolverControl solver_control_F(1000, 1e-2 * src.block(0).l2_norm());
      SolverGMRES<TrilinosWrappers::MPI::Vector> solver_cg_F(solver_control_F);
      solver_cg_F.solve (*F,
                          dst.block(0),
                          tmp2,
                          preconditioner_F);

      // block 1

      dst.block(1).sadd(-1.0, tmp1); // non sicuro

    }

  protected:
    // F matrix.
    const TrilinosWrappers::SparseMatrix *F;

    // Preconditioner used for the F block.
    TrilinosWrappers::PreconditionILU preconditioner_F;

    // Pressure mass matrix.
    const TrilinosWrappers::SparseMatrix *pressure_mass;

    // Preconditioner used for the pressure block.
    TrilinosWrappers::PreconditionILU preconditioner_pressure;

    // B matrix.
    const TrilinosWrappers::SparseMatrix *B;

        // B matrix.
    const TrilinosWrappers::SparseMatrix *Bt;

    // Ap matrix.
    const TrilinosWrappers::SparseMatrix *Ap;

    // Preconditioner used for the Ap block.
    TrilinosWrappers::PreconditionILU preconditioner_Ap;

    // Fp matrix.
    const TrilinosWrappers::SparseMatrix *Fp;

    // Temporary vectors.
    mutable TrilinosWrappers::MPI::Vector tmp1;
    mutable TrilinosWrappers::MPI::Vector tmp2;
  };

  // Constructor.
  NavierStokes(const unsigned int &N_,
               const unsigned int &degree_velocity_,
               const unsigned int &degree_pressure_,
               const double &      T_,
               const double &      deltat_,
               const double &      theta_)
    : mpi_size(Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD))
    , mpi_rank(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    , pcout(std::cout, mpi_rank == 0)
    , T(T_)
    , N(N_)
    , degree_velocity(degree_velocity_)
    , degree_pressure(degree_pressure_)
    , deltat(deltat_)
    , theta(theta_)
    , mesh(MPI_COMM_WORLD)
  {}

  // Setup system.
  void
  setup();

  // Solve problem.
  void
  solve();

protected:
  // Assemble recurrent matrices. 
  // We also assemble the pressure mass matrix (needed for the
  // preconditioner).
  void
  assemble_matrices();

  // Assemble system.
  void
  assemble_system();

  // Solve system at current time.
  void
  solve_time_step();

  // Output results.
  void
  output(const unsigned int &time_step, const double &time);

  // MPI parallel. /////////////////////////////////////////////////////////////

  // Number of MPI processes.
  const unsigned int mpi_size;

  // This MPI process.
  const unsigned int mpi_rank;

  // Parallel output stream.
  ConditionalOStream pcout;

  // Problem definition. ///////////////////////////////////////////////////////

  // Kinematic viscosity [m2/s].
  const double nu = 0.01;

  // Outlet pressure [Pa].
  const double p_out = 0.0;

  // Forcing term.
  ForcingTerm forcing_term;

  // Initial condition.
  FunctionU0 u_0;

  // Inlet velocity.
  InletVelocity inlet_velocity;

  // Current time.
  double time;

  // Final time.
  const double T;

  // Discretization. ///////////////////////////////////////////////////////////

  // Mesh refinement.
  const unsigned int N;

  // Polynomial degree used for velocity.
  const unsigned int degree_velocity;

  // Polynomial degree used for pressure.
  const unsigned int degree_pressure;

  // Time step.
  const double deltat;

  // Theta parameter of the theta method.
  const double theta;

  // Mesh.
  parallel::fullydistributed::Triangulation<dim> mesh;

  // Finite element space.
  std::unique_ptr<FiniteElement<dim>> fe;

  // Quadrature formula.
  std::unique_ptr<Quadrature<dim>> quadrature;

  // Quadrature formula for face integrals.
  std::unique_ptr<Quadrature<dim - 1>> quadrature_face;

  // DoF handler.
  DoFHandler<dim> dof_handler;

  // DoFs owned by current process.
  IndexSet locally_owned_dofs;

  // DoFs owned by current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_owned_dofs;

  // DoFs relevant to the current process (including ghost DoFs).
  IndexSet locally_relevant_dofs;

  // DoFs relevant to current process in the velocity and pressure blocks.
  std::vector<IndexSet> block_relevant_dofs;

  // System matrix.
  TrilinosWrappers::BlockSparseMatrix system_matrix;

  // Pressure mass matrix, needed for preconditioning. We use a block matrix for
  // convenience, but in practice we only look at the pressure-pressure block.
  TrilinosWrappers::BlockSparseMatrix pressure_mass;

  
  //preconditioner
  TrilinosWrappers::BlockSparseMatrix Fp_matrix;

  TrilinosWrappers::BlockSparseMatrix inverse_diagonal_mass_matrix;

  // Matrix on the right-hand side (M / deltat - theta * A).
  TrilinosWrappers::BlockSparseMatrix lhs_matrix;
  
  // Matrix on the right-hand side (M / deltat - (1 - theta) * A).
  TrilinosWrappers::BlockSparseMatrix rhs_matrix;

  // Right-hand side vector in the linear system.
  TrilinosWrappers::MPI::BlockVector system_rhs;

  // System solution (without ghost elements).
  TrilinosWrappers::MPI::BlockVector solution_owned;

  // System solution (including ghost elements).
  TrilinosWrappers::MPI::BlockVector solution;
};

#endif