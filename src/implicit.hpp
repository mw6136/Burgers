#include "grid.hpp"

#ifndef IMPLICIT_HPP_
#define IMPLICIT_HPP_

using namespace arma;

void
printArray (mat A, int nx)
{

  for (int i = 0; i < nx; ++i)
    {
      for (int j = 0; j < nx; ++j)
        {
          std::cout << A (i, j) << ", ";
        }
      std::cout << "\n";
    }
}

mat
fillJacobian (mat dfdu, vec uk, int nx1, double dx, double dt, double v,
              std::string x_disc)
{

  for (int j = 0; j < nx1; ++j)
    {

      if (x_disc == "center")
        {
          if (j == 0)
            {
              dfdu (j, j)
                  = 1
                    + dt / dx * ((uk (j + 1) - uk (nx1 - 1)) / 2 + 2 * v / dx);
              dfdu (j, j + 1) = dt / dx * (uk (j) / 2 - v / dx);
            }
          else if (j == nx1 - 1)
            {
              dfdu (j, j)
                  = 1 + dt / dx * ((uk (0) - uk (j - 1)) / 2 + 2 * v / dx);
              dfdu (j, j - 1) = -dt / dx * (uk (j) / 2 + v / dx);
            }
          else
            {
              dfdu (j, j)
                  = 1 + dt / dx * ((uk (j + 1) - uk (j - 1)) / 2 + 2 * v / dx);
              dfdu (j, j + 1) = dt / dx * (uk (j) / 2 - v / dx);
              dfdu (j, j - 1) = -dt / dx * (uk (j) / 2 + v / dx);
            }
        }
      else if (x_disc == "upwind")
        {
          if (uk (j) > 0)
            {
              if (j == 0)
                {
                  dfdu (j, j) = 1 + dt * (uk (j) - uk (j)) / (2 * dx)
                                - dt * v * (-2 / (dx * dx));
                  dfdu (j, j + 1) = -dt * v * (1 / (dx * dx));
                }
              else if (j == nx1 - 1)
                {
                  dfdu (j, j) = 1 + dt * (uk (j) - uk (j - 1)) / (2 * dx)
                                - dt * v * (-2 / (dx * dx));
                  dfdu (j, j - 1)
                      = -dt * (uk (j)) / (2 * dx) - dt * v * (1 / (dx * dx));
                }
              else
                {
                  dfdu (j, j) = 1 + dt * (uk (j) - uk (j - 1)) / (2 * dx)
                                - dt * v * (-2 / (dx * dx));
                  dfdu (j, j + 1) = -dt * v * (1 / (dx * dx));
                  dfdu (j, j - 1)
                      = -dt * (uk (j)) / (2 * dx) - dt * v * (1 / (dx * dx));
                }
            }
          else
            {
              if (j == 0)
                {
                  dfdu (j, j) = 1 + dt * (uk (j + 1) - uk (j)) / (2 * dx)
                                - dt * v * (-2 / (dx * dx));
                  dfdu (j, j + 1)
                      = dt * (uk (j)) / (2 * dx) - dt * v * (1 / (dx * dx));
                }
              else if (j == nx1 - 1)
                {
                  dfdu (j, j) = 1 + dt * (uk (j) - uk (j)) / (2 * dx)
                                - dt * v * (-2 / (dx * dx));
                  dfdu (j, j - 1) = -dt * v * (1 / (dx * dx));
                }
              else
                {
                  dfdu (j, j) = 1 + dt * (uk (j + 1) - uk (j)) / (2 * dx)
                                - dt * v * (-2 / (dx * dx));
                  dfdu (j, j + 1)
                      = dt * (uk (j)) / (2 * dx) - dt * v * (1 / (dx * dx));
                  dfdu (j, j - 1) = -dt * v * (1 / (dx * dx));
                }
            }
        }
    }
  return dfdu;
}

vec
NewtonsMethod (Mesh *pmesh, double dt, double dx, double v, std::string x_disc,
               double tol = 1e-4)
{
  /*
  Performs Newton's method for implicit euler time integration
  */

  // Initialize empty
  int nx1 = pmesh->m_nx1;

  mat dfdu (nx1, nx1);
  vec fuk (nx1);

  vec uk (nx1);
  vec ukp1 (nx1);

  vec umabsu (nx1);
  vec upabsu (nx1);

  mat umabsu_diag (nx1, nx1);
  mat upabsu_diag (nx1, nx1);

  double error = 1e3;

  // Fill initial guesses for uk and ukp1
  for (int i = 0; i < nx1; ++i)
    {
      uk (i) = pmesh->u (i);
      ukp1 (i) = pmesh->u (i);
    }

  // Loop over iterations to get u_np1
  while (error > tol)
    {

      mat udiag = diagmat (uk);

      if (x_disc == "upwind")
        {
          // Get diagonal matrices for upwind derivative operator
          for (int i = 0; i < nx1; ++i)
            {
              umabsu (i) = (uk (i) - abs (uk (i))) / 2;
              upabsu (i) = (uk (i) + abs (uk (i))) / 2;
            }

          umabsu_diag = diagmat (umabsu);
          upabsu_diag = diagmat (upabsu);
        }

      // Get the expression for F(Uk)
      if (x_disc == "center")
        {
          fuk = uk - pmesh->u
                + dt
                      * (1 / (2 * dx) * udiag * pmesh->ddx_ctr * uk
                         - v / (dx * dx) * pmesh->ddx2_op * uk);
        }
      else if (x_disc == "upwind")
        {
          fuk = uk - pmesh->u
                + dt
                      * (1 / (2 * dx)
                             * (upabsu_diag * pmesh->ddx_bwd * uk
                                + umabsu_diag * pmesh->ddx_fwd * uk)
                         - v / (dx * dx) * pmesh->ddx2_op * uk);
        }

      // Fill Jacobian
      dfdu = fillJacobian (dfdu, uk, nx1, dx, dt, v, x_disc);

      // Calculate u_k+1
      ukp1 = uk - inv (dfdu) * fuk;

      // Check for convergence
      error = norm ((ukp1 - uk), 2);

      // Update for next iteration
      uk = ukp1;
    }

  return uk;
}

#endif