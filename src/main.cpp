#include <armadillo>
#include <iomanip>
#include <iostream>

#include "grid.hpp"
#include "implicit.hpp"

using namespace arma;

void
ProblemGenerator (Mesh *pmesh, int nx1, double xmin, double xmax, double dx)
{

  for (int i = 0; i < nx1; ++i)
    {

      pmesh->x1s (i) = xmin + i * dx;

      pmesh->u (i)
          = sin (pmesh->x1s (i)) * exp (-pow ((pmesh->x1s (i) - M_PI), 2.0));
    }
}

int
main (int argc, char *argv[])
{

  int nx1 = atoi (argv[4]);
  double v = atof (argv[3]);

  double xmax = 2.0 * M_PI;
  double xmin = 0.0;

  double dx = (xmax - xmin) / nx1;

  double dt;
  // if (v == 0) {
  //     dt = 0.005;
  // } else {
  //     dt = 0.25 * (dx * dx) / v;
  // }
  // std::cout<< dt<< std::endl;
  dt = 0.05;

  double tmax = 40.0 * M_PI;

  int maxiters = 1000000;
  int savefreq = 400;

  std::string time_int, spat_disc;

  if (atoi (argv[1]) == 0)
    {
      time_int = "explicit";
    }
  else
    {
      time_int = "implicit";
    }

  if (atoi (argv[2]) == 0)
    {
      spat_disc = "center";
    }
  else
    {
      spat_disc = "upwind";
    }

  std::ostringstream savepref;

  savepref << time_int << "." << spat_disc << "." << v << ".";

  Mesh *pmesh;
  pmesh = new Mesh (nx1, spat_disc);

  vec temp_conv_array (nx1);
  vec temp_diff_array (nx1);

  ProblemGenerator (pmesh, nx1, xmin, xmax, dx);
  pmesh->fillDiscreteMatrix ();

  int iter = 0;
  double t = 0;
  int saveiter = 0;
  for (iter; iter < maxiters; ++iter)
    {

      if (t > tmax)
        {
          break;
        }

      if (time_int == "explicit")
        {
          // Explicit implementation

          if (spat_disc == "center")
            {
              // Construct diagonal velocity matrix for nonlinearity
              mat udiag = diagmat (pmesh->u);

              temp_conv_array
                  = pmesh->u
                    + dt
                          * (-1 / (2 * dx) * udiag * pmesh->ddx_ctr * pmesh->u
                             + v / (dx * dx) * pmesh->ddx2_op * pmesh->u);
            }
          else if (spat_disc == "upwind")
            {
              // Construct diagonal velocity matrix for nonlinearity
              pmesh->fillAbsVal ();

              mat umabsu_diag = diagmat (pmesh->umabsu);
              mat upabsu_diag = diagmat (pmesh->upabsu);

              temp_conv_array
                  = pmesh->u
                    + dt
                          * (-1 / (2 * dx)
                                 * (upabsu_diag * pmesh->ddx_bwd * pmesh->u
                                    + umabsu_diag * pmesh->ddx_fwd * pmesh->u)
                             + v / (dx * dx) * pmesh->ddx2_op * pmesh->u);
            }
        }
      else if (time_int == "implicit")
        {

          // Do Newton's method to get u for next timestep
          temp_conv_array = NewtonsMethod (pmesh, dt, dx, v, spat_disc);
        }

      // Swap arrays
      for (int i = 0; i < nx1; ++i)
        {
          pmesh->u (i) = temp_conv_array (i);
        }

      if (iter % savefreq == 0)
        {
          std::cout << "Iteration : " << iter << ", Time : " << t
                    << ", dt = " << dt << std::endl;
          pmesh->saveOutput (saveiter, t, savepref.str ());
          saveiter += 1;
        }

      t += dt;
    }

  return 0;
}