#include "Navier-Stokes_TD.hpp"

// Main function.
int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv);  

  const unsigned int degree_velocity = 2;
  const unsigned int degree_pressure = 1;
  const double T     = 6.0;
  const double theta = 1.0;
  const std::vector<double> deltat_vector = {0.125};

  for (const auto &deltat : deltat_vector)
  {   
    NavierStokes problem(degree_velocity, degree_pressure, T, deltat, theta);
    problem.setup();
    problem.solve();
  }

  return 0;
}