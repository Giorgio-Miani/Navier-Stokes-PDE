#include "Navier-Stokes_TD.hpp"

void
NavierStokes::setup()
{
  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;
    GridIn<dim> grid_in;

    grid_in.attach_triangulation(mesh_serial);
    const std::string mesh_file_name = "../mesh/mesh-0.1.msh";
    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space.
  {
    pcout << "Initializing the finite element space" << std::endl;

    const FE_SimplexP<dim> fe_scalar_velocity(degree_velocity);
    const FE_SimplexP<dim> fe_scalar_pressure(degree_pressure);
    fe = std::make_unique<FESystem<dim>>(fe_scalar_velocity,
                                         dim,
                                         fe_scalar_pressure,
                                         1);

    pcout << "  Velocity degree:           = " << fe_scalar_velocity.degree
          << std::endl;
    pcout << "  Pressure degree:           = " << fe_scalar_pressure.degree
          << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;

    quadrature_face = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);

    pcout << "  Quadrature points per face = " << quadrature_face->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We want to reorder DoFs so that all velocity DoFs come first, and then
    // all pressure DoFs.
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);

    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    // Besides the locally owned and locally relevant indices for the whole
    // system (velocity and pressure), we will also need those for the
    // individual velocity and pressure blocks.
    std::vector<types::global_dof_index> dofs_per_block =
      DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
    const unsigned int n_u = dofs_per_block[0];
    const unsigned int n_p = dofs_per_block[1];

    block_owned_dofs.resize(2);
    block_relevant_dofs.resize(2);
    block_owned_dofs[0]    = locally_owned_dofs.get_view(0, n_u);
    block_owned_dofs[1]    = locally_owned_dofs.get_view(n_u, n_u + n_p);
    block_relevant_dofs[0] = locally_relevant_dofs.get_view(0, n_u);
    block_relevant_dofs[1] = locally_relevant_dofs.get_view(n_u, n_u + n_p);

    pcout << "  Number of DoFs: " << std::endl;
    pcout << "    velocity = " << n_u << std::endl;
    pcout << "    pressure = " << n_p << std::endl;
    pcout << "    total    = " << n_u + n_p << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // Velocity DoFs interact with other velocity DoFs (the weak formulation has
    // terms involving u times v), and pressure DoFs interact with velocity DoFs
    // (there are terms involving p times v or u times q). However, pressure
    // DoFs do not interact with other pressure DoFs (there are no terms
    // involving p times q). We build a table to store this information, so that
    // the sparsity pattern can be built accordingly.
    Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
    for (unsigned int c = 0; c < dim + 1; ++c)
      {
        for (unsigned int d = 0; d < dim + 1; ++d)
          {
            if (c == dim && d == dim) // pressure-pressure term
              coupling[c][d] = DoFTools::none;
            else // other combinations
              coupling[c][d] = DoFTools::always;
          }
      }

    TrilinosWrappers::BlockSparsityPattern sparsity(block_owned_dofs,
                                                    MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity);
    sparsity.compress();

    //Sparsity pattern for velocity mass matrix
    for (unsigned int c = 0; c < dim + 1; ++c)
      {
        for (unsigned int d = 0; d < dim + 1; ++d)
          {
            if (c == dim || d == dim) // pressure-pressure term
              coupling[c][d] = DoFTools::none;
            else // other combinations
              coupling[c][d] = DoFTools::always;
          }
      }

    TrilinosWrappers::BlockSparsityPattern sparsity_velocity_mass(block_owned_dofs,
                                                    MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, coupling, sparsity_velocity_mass);
    sparsity_velocity_mass.compress();

    // We also build a sparsity pattern for the pressure mass matrix.
    for (unsigned int c = 0; c < dim + 1; ++c)
      {
        for (unsigned int d = 0; d < dim + 1; ++d)
          {
            if (c == dim && d == dim) // pressure-pressure term
              coupling[c][d] = DoFTools::always;
            else // other combinations
              coupling[c][d] = DoFTools::none;
          }
      }
    TrilinosWrappers::BlockSparsityPattern sparsity_pressure_mass(
      block_owned_dofs, MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler,
                                    coupling,
                                    sparsity_pressure_mass);
    sparsity_pressure_mass.compress();

    system_matrix.reinit(sparsity);
    rhs_matrix.reinit(sparsity);
    lhs_matrix.reinit(sparsity);
    pressure_mass.reinit(sparsity_pressure_mass);
    //Initialize preconditioner Fp, not user if the sparsity pattern is right
    inverse_diagonal_mass_matrix.reinit(sparsity_velocity_mass);
    Fp_matrix.reinit(sparsity_pressure_mass);
    intermediate_fp_matrix.reinit(sparsity_pressure_mass);
    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(block_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(block_owned_dofs, MPI_COMM_WORLD);
    solution.reinit(block_owned_dofs, block_relevant_dofs, MPI_COMM_WORLD);

    pcout << "-----------------------------------------------" << std::endl;
    pcout << "Reynolds number: " << Re << std::endl;
  }

  // Initialize drag and lift vectors
  {
    drag_coefficients.reserve(static_cast<int>(T/deltat + 1.0));
    lift_coefficients.reserve(static_cast<int>(T/deltat + 1.0));
  }
}

void
NavierStokes::assemble_matrices()
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the recurrent matrices" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);

  FullMatrix<double> cell_lhs_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_rhs_matrix(dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_pressure_mass_matrix(dofs_per_cell, dofs_per_cell);         //Per precondizionatore
  FullMatrix<double> cell_inverse_diagonal_mass_matrix(dofs_per_cell, dofs_per_cell); //Per precondizionatore
  FullMatrix<double> cell_Fp_matrix(dofs_per_cell, dofs_per_cell);                    //Per precondizionatore

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  lhs_matrix    = 0.0;
  rhs_matrix    = 0.0;
  pressure_mass = 0.0;
  inverse_diagonal_mass_matrix = 0.0;
  intermediate_fp_matrix = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_lhs_matrix           = 0.0;
      cell_rhs_matrix           = 0.0;
      cell_pressure_mass_matrix = 0.0;
      cell_inverse_diagonal_mass_matrix = 0.0;

      for (unsigned int q = 0; q < n_q; ++q)
        {
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                { 
                  // M/deltaT
                  cell_lhs_matrix(i, j) += 
                    scalar_product(fe_values[velocity].value(i, q),
                                   fe_values[velocity].value(j, q)) /
                    deltat * fe_values.JxW(q);

                  // A*theta                 
                  cell_lhs_matrix(i, j) +=
                    nu * theta * 
                    scalar_product(fe_values[velocity].gradient(i, q),
                                   fe_values[velocity].gradient(j, q)) *
                    fe_values.JxW(q);

                  // Pressure term in the momentum equation.
                  cell_lhs_matrix(i, j) -= fe_values[velocity].divergence(i, q) *
                                           fe_values[pressure].value(j, q) *
                                           fe_values.JxW(q);

                  // Pressure term in the continuity equation.
                  cell_lhs_matrix(i, j) -= fe_values[velocity].divergence(j, q) *
                                           fe_values[pressure].value(i, q) *
                                           fe_values.JxW(q);

                  // M/deltaT
                  cell_rhs_matrix(i, j) += 
                    scalar_product(fe_values[velocity].value(i, q),
                                   fe_values[velocity].value(j, q)) /
                    deltat * fe_values.JxW(q);

                  // A*(1-theta)               
                  cell_rhs_matrix(i, j) +=
                    nu * (1.0 - theta) * 
                    scalar_product(fe_values[velocity].gradient(i, q),
                                   fe_values[velocity].gradient(j, q)) *
                    fe_values.JxW(q);

                  // Pressure mass matrix.
                  cell_pressure_mass_matrix(i, j) +=
                    fe_values[pressure].value(i, q) *
                    fe_values[pressure].value(j, q) / nu * fe_values.JxW(q);

                  ////////////////////////// PCD Preconditioner //////////////////////////
                  // M/deltaT
                  cell_Fp_matrix(i, j) += fe_values[pressure].value(i, q) *
                                          fe_values[pressure].value(j, q) /
                                          deltat * fe_values.JxW(q);
                  
                  // A*theta                 
                  cell_Fp_matrix(i, j) += nu * theta * 
                                          scalar_product(fe_values[pressure].gradient(i, q),
                                                         fe_values[pressure].gradient(j, q)) *
                                          fe_values.JxW(q);     

                  cell_inverse_diagonal_mass_matrix(i, j) += scalar_product(fe_values[velocity].value(i, q),
                                                             fe_values[velocity].value(j, q)) * fe_values.JxW(q);                             
                }
            }
        }

      cell->get_dof_indices(dof_indices);

      lhs_matrix.add(dof_indices, cell_lhs_matrix);
      rhs_matrix.add(dof_indices, cell_rhs_matrix);
      pressure_mass.add(dof_indices, cell_pressure_mass_matrix);
      inverse_diagonal_mass_matrix.add(dof_indices, cell_inverse_diagonal_mass_matrix);
      intermediate_fp_matrix.add(dof_indices, cell_Fp_matrix);
    }

  lhs_matrix.compress(VectorOperation::add);
  rhs_matrix.compress(VectorOperation::add);
  pressure_mass.compress(VectorOperation::add);
  inverse_diagonal_mass_matrix.compress(VectorOperation::add);
  intermediate_fp_matrix.compress(VectorOperation::add);
}

void
NavierStokes::assemble_system()
{
  pcout << "===============================================" << std::endl;
  pcout << "Assembling the system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();
  // const unsigned int n_q_face      = quadrature_face->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                          update_quadrature_points | update_JxW_values);
  
  // FEFaceValues<dim> fe_face_values(*fe,
  //                                  *quadrature_face,
  //                                  update_values | update_normal_vectors |
  //                                  update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  FullMatrix<double> cell_Fp_matrix(dofs_per_cell, dofs_per_cell); //Per precondizionatore


  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  system_matrix = 0.0;
  system_rhs    = 0.0;
  Fp_matrix = 0.0;

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  std::vector<Tensor<1, dim>> velocity_loc(n_q);
  std::vector<Tensor<2, dim>> velocity_gradient_loc(n_q);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (!cell->is_locally_owned())
        continue;

      fe_values.reinit(cell);

      cell_matrix = 0.0;
      cell_rhs    = 0.0;
      cell_Fp_matrix = 0.0;

      fe_values[velocity].get_function_values(solution, velocity_loc);
      fe_values[velocity].get_function_gradients(solution, velocity_gradient_loc);

      for (unsigned int q = 0; q < n_q; ++q)
        {

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                { 
                  //Non-linear term
                  // cell_matrix(i, j) += 0.5 * velocity_loc[q] * fe_values[velocity].gradient(j, q) * fe_values[velocity].value(i, q) *
                  //                      fe_values.JxW(q);
                  // cell_matrix(i, j) += 0.5 * present_velocity_divergence * fe_values[velocity].value(j, q) * fe_values[velocity].value(i, q) *
                  //                      fe_values.JxW(q);

                  cell_matrix(i, j) += rho * velocity_loc[q] * fe_values[velocity].gradient(j, q) * fe_values[velocity].value(i, q) *
                                       fe_values.JxW(q);

                  ////////////////////////// PCD Preconditioner //////////////////////////
                  cell_Fp_matrix(i, j) += rho * velocity_loc[q] * fe_values[pressure].gradient(j, q) * fe_values[pressure].value(i, q) *
                        fe_values.JxW(q);                                  
                }

              // Compute F(t_n)
              forcing_term.set_time(time - deltat);
              Vector<double> forcing_term_pre(dim);
              forcing_term.vector_value(fe_values.quadrature_point(q), forcing_term_pre);
              Tensor<1, dim> forcing_term_tensor_pre;
              for (unsigned int d = 0; d < dim; ++d)
                forcing_term_tensor_pre[d] = forcing_term_pre[d];

              // Compute F(t_n+1)
              forcing_term.set_time(time);
              Vector<double> forcing_term_loc(dim);
              forcing_term.vector_value(fe_values.quadrature_point(q), forcing_term_loc);
              Tensor<1, dim> forcing_term_tensor;
              for (unsigned int d = 0; d < dim; ++d)
                forcing_term_tensor[d] = forcing_term_loc[d];

              // Forcing term.
              // theta*F(t_n+1)
              cell_rhs(i) += theta * scalar_product(forcing_term_tensor, fe_values[velocity].value(i, q)) *
                             fe_values.JxW(q);

              // (1-theta)*F(t_n)
              cell_rhs(i) += (1 - theta) * scalar_product(forcing_term_tensor_pre, fe_values[velocity].value(i, q)) *
                             fe_values.JxW(q);
            }
        }

      // // Boundary integral for Neumann BCs.
      // if (cell->at_boundary())
      //   {
      //     for (unsigned int f = 0; f < cell->n_faces(); ++f)
      //       {
      //         if (cell->face(f)->at_boundary() &&
      //             cell->face(f)->boundary_id() == 0)
      //           {
      //             fe_face_values.reinit(cell, f);

      //             for (unsigned int q = 0; q < n_q_face; ++q)
      //               {
      //                 for (unsigned int i = 0; i < dofs_per_cell; ++i)
      //                   {
      //                     cell_rhs(i) +=
      //                       -p_out *
      //                       scalar_product(fe_face_values.normal_vector(q),
      //                                      fe_face_values[velocity].value(i,
      //                                                                     q)) *
      //                       fe_face_values.JxW(q);
      //                   }
      //               }
      //           }
      //       }
      //   }

      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);
      Fp_matrix.add(dof_indices, cell_Fp_matrix);
    }

  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);
  Fp_matrix.compress(VectorOperation::add);

  // rhs = (M/deltaT + A*(1-theta))*u_n + (1-theta)*F(t_n) + theta*F(t_n+1) 
  rhs_matrix.vmult_add(system_rhs, solution_owned);
  system_matrix.add(1.0, lhs_matrix);
  Fp_matrix.add(1.0, intermediate_fp_matrix);

  // Dirichlet boundary conditions.
  {
    std::map<types::global_dof_index, double>           boundary_values;
    std::map<types::boundary_id, const Function<dim> *> boundary_functions;
    Functions::ZeroFunction<dim> zero_function(dim + 1);

    // Inflow
    inlet_velocity.set_time(time);
    boundary_functions[0] = &inlet_velocity;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                               {true, true, true, false}));
    boundary_functions.clear();

    // Cylinder
    boundary_functions[6] = &zero_function;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                               {true, true, true, false}));
    boundary_functions.clear();

    // Up/Down
    boundary_functions[2] = &zero_function;
    boundary_functions[4] = &zero_function;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                               {false, true, false, false})); 
    boundary_functions.clear();

    // Left/Right
    boundary_functions[3] = &zero_function;
    boundary_functions[5] = &zero_function;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values,
                                             ComponentMask(
                                               {false, false, true, false})); 

    // Boundary conditions on all walls for tunnel test case to test drag and lift coeff
    // boundary_functions[2] = &zero_function;
    // boundary_functions[4] = &zero_function;
    // boundary_functions[3] = &zero_function;
    // boundary_functions[5] = &zero_function;    
    // VectorTools::interpolate_boundary_values(dof_handler,
    //                                          boundary_functions,
    //                                          boundary_values,
    //                                          ComponentMask(
    //                                            {true, true, true, false}));

    boundary_functions.clear();

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, false);
  }
}

void
NavierStokes::solve_time_step()
{
  pcout << "===============================================" << std::endl;

  SolverControl solver_control(10000, 1e-6 * system_rhs.l2_norm());
  SolverGMRES<TrilinosWrappers::MPI::BlockVector> solver(solver_control);

  ///////////////////// PCD PRECONDITIONER /////////////////////
  // TrilinosWrappers::MPI::Vector diagonal;
  //     {
  //     for (unsigned int i = 0; i < inverse_diagonal_mass_matrix.block(0,0).m(); ++i)
  //       {
  //         auto value = inverse_diagonal_mass_matrix.block(0,0).diag_element(i);
  //         diagonal(i) = (1.0/value);
  //       }
  //     }
  // system_matrix.block(1,0).mmult(inverse_diagonal_mass_matrix.block(0,0), system_matrix.block(0,1), diagonal);
  // PreconditionBlockPCD preconditioner;
  // preconditioner.initialize(system_matrix.block(0,0), pressure_mass.block(1,1), system_matrix.block(0, 1), inverse_diagonal_mass_matrix.block(0,0), Fp_matrix.block(1,1));  

  ///////////////////// DIAGONAL PRECONDITIONER /////////////////////
  PreconditionBlockDiagonal preconditioner;
  preconditioner.initialize(system_matrix.block(0, 0),
                            pressure_mass.block(1, 1));

  ///////////////////// TRIANGULAR PRECONDITIONER /////////////////////
  // PreconditionBlockTriangular preconditioner;
  // preconditioner.initialize(system_matrix.block(0, 0),
  //                           pressure_mass.block(1, 1),
  //                           system_matrix.block(1, 0));

  pcout << "Solving the linear system" << std::endl;
  solver.solve(system_matrix, solution_owned, system_rhs, preconditioner);
  pcout << "  " << solver_control.last_step() << " GMRES iterations"
        << std::endl;

  solution = solution_owned;
}

void
NavierStokes::calculate_coefficients()
{
  pcout << "===============================================" << std::endl;
  pcout << "Calculating coefficients" << std::endl;

  const unsigned int n_q_face = quadrature_face->size();

  FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_face,
                                   update_values | update_gradients | update_normal_vectors |
                                   update_JxW_values);
  
  double f_d = 0.0;
  double f_l = 0.0; 

  FEValuesExtractors::Vector velocity(0);
  FEValuesExtractors::Scalar pressure(dim);

  std::vector<double> pressure_loc(n_q_face);
  std::vector<Tensor<2, dim>> velocity_gradient_loc(n_q_face);

  Tensor<2, dim> fluid_pressure;
  Tensor<1, dim> normal_vector;
  Tensor<1, dim> forces;

  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    if (!cell->is_locally_owned())
      continue;

    if(cell->at_boundary()){
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
      {
        if (cell->face(f)->at_boundary() &&
            cell->face(f)->boundary_id() == 6)
        {
          fe_face_values.reinit(cell, f);
          fe_face_values[pressure].get_function_values(solution, pressure_loc);
          fe_face_values[velocity].get_function_gradients(solution, velocity_gradient_loc);

          for (unsigned int q = 0; q < n_q_face; ++q)
          {

            normal_vector = -fe_face_values.normal_vector(q);

            fluid_pressure[0][0] = pressure_loc[q];
            fluid_pressure[1][1] = pressure_loc[q];

            forces = (rho * nu * velocity_gradient_loc[q] - fluid_pressure) * normal_vector * fe_face_values.JxW(q);
            f_d += forces[0];
            f_l += forces[1];
          }
        }
      }
    }
  }

  if(mpi_rank == 0)
  { 
    double temp_f_d, temp_f_l;
    for(unsigned int i = 1; i < mpi_size; ++i){
      MPI_Recv(&temp_f_d, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&temp_f_l, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 

      f_d += temp_f_d;
      f_l += temp_f_l;
    }
    drag_coefficients.push_back(multiplicative_const * f_d);
    lift_coefficients.push_back(multiplicative_const * f_l);
  }
  else{
    MPI_Send(&f_d, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&f_l, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  }
}

void
NavierStokes::output_coefficients()
{
  const unsigned int mpi_rank = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  pcout << "===============================================" << std::endl;
  
  if(mpi_rank == 0)
  {
    pcout << "Creating drag_coefficient.csv" << std::endl;
    std::ofstream drag_coefficient_file("drag_coefficient.csv");
    drag_coefficient_file << "T,Cd" << std::endl;

    pcout << "Creating lift_coefficient.csv" << std::endl;
    std::ofstream lift_coefficient_file("lift_coefficient.csv");
    lift_coefficient_file << "T,Cl" << std::endl;

    unsigned int head = 0;
    for (double t = 0.0; t <= T; t = t + deltat)
    {
      drag_coefficient_file << t << "," << drag_coefficients[head] << std::endl;
      lift_coefficient_file << t << "," << lift_coefficients[head] << std::endl;
      head++;
    }
  }

  pcout << "Coefficients written to file" << std::endl;
  pcout << "===============================================" << std::endl;
}

void
NavierStokes::solve()
{ 
  pcout << "===============================================" << std::endl;
  TimerOutput timer (pcout, TimerOutput::never, TimerOutput::wall_times);
  
  time = 0.0;

  timer.enter_subsection ("Matrix assembly");
  assemble_matrices();
  timer.leave_subsection();

  // Apply the initial condition.
  {
    pcout << "Applying the initial condition" << std::endl;

    u_0.set_time(time);
    VectorTools::interpolate(dof_handler, u_0, solution_owned);
    solution = solution_owned;

    timer.enter_subsection ("Coefficients calculation");
    calculate_coefficients();
    timer.leave_subsection();
    timer.enter_subsection ("Output");
    output(0, 0.0);
    timer.leave_subsection();
    pcout << "-----------------------------------------------" << std::endl;
  }

  unsigned int time_step = 0;
  double intpart, fractpart;

  while (time < T)
    {
      time += deltat;
      ++time_step;
      fractpart = std::modf(time, &intpart);

      pcout << "n = " << std::setw(3) << time_step << ", t = " << std::setw(5)
            << time << ":\n" << std::flush;

      timer.enter_subsection ("System assembly");
      assemble_system();
      timer.leave_subsection();
      timer.enter_subsection ("System solving");
      solve_time_step();
      timer.leave_subsection();
      timer.enter_subsection ("Coefficients calculation");
      calculate_coefficients();
      timer.leave_subsection();
      timer.enter_subsection ("Output");
      if(fractpart == 0.0 || fractpart == 0.5) output(time_step, time);
      timer.leave_subsection();
    }
  timer.enter_subsection ("Output coefficients");
  output_coefficients();
  timer.leave_subsection();

  timer.print_wall_time_statistics(MPI_COMM_WORLD);
}

void
NavierStokes::output(const unsigned int &time_step, const double &time)
{
  pcout << "===============================================" << std::endl;

  DataOut<dim> data_out;

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
  data_component_interpretation.push_back(
    DataComponentInterpretation::component_is_scalar);
  std::vector<std::string> names = {"velocity",
                                    "velocity",
                                    "velocity",
                                    "pressure"};

  data_out.add_data_vector(dof_handler,
                           solution,
                           names,
                           data_component_interpretation);

  // std::vector<unsigned int> partition_int(mesh.n_active_cells());
  // GridTools::get_subdomain_association(mesh, partition_int);
  // const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  // data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  std::string output_file_name = "output-" + std::to_string(time_step);

  DataOutBase::DataOutFilter data_filter(
    DataOutBase::DataOutFilterFlags(/*filter_duplicate_vertices = */ true,
                                    /*xdmf_hdf5_output = */ true));
  data_out.write_filtered_data(data_filter);
  data_out.write_hdf5_parallel(data_filter,
                               output_file_name + ".h5",
                               MPI_COMM_WORLD);

  std::vector<XDMFEntry> xdmf_entries({data_out.create_xdmf_entry(
    data_filter, "output-" + std::to_string(time_step) + ".h5", time, MPI_COMM_WORLD)});
  data_out.write_xdmf_file(xdmf_entries,
                           output_file_name + ".xdmf",
                           MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;
  pcout << "===============================================" << std::endl;
}