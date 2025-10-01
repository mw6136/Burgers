#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#ifndef GRID_HPP_
#define GRID_HPP_

using namespace arma;

std::string
ZeroPadNumber (int num)
{
  std::ostringstream ss;
  ss << std::setw (4) << std::setfill ('0') << num;
  std::string result = ss.str ();
  if (result.length () > 7)
    {
      result.erase (0, result.length () - 7);
    }
  return result;
}

std::string
getFname (int iter, std::string prefix)
{

  std::string fname;
  std::stringstream ss;
  std::stringstream ssfname;
  std::string newnum = ZeroPadNumber (iter);

  ss << prefix << newnum << ".csv";

  ss >> fname;

  return fname;
}

class Mesh
{

public:
  int m_nx1;
  std::string m_x_disc_meth;

  // Vector of u and x values
  vec x1s;
  vec u;

  // Second order deriv operator
  mat ddx2_op;

  // First order deriv (forward and backward) operators
  mat ddx_fwd;
  mat ddx_bwd;

  // Centered first order deriv
  mat ddx_ctr;

  // Arrays to store absolute value coefficients
  vec umabsu;
  vec upabsu;

  Mesh (int nx1, std::string x_disc_meth)
      : u (nx1), x1s (nx1), umabsu (nx1), upabsu (nx1), ddx2_op (nx1, nx1),
        ddx_fwd (nx1, nx1), ddx_ctr (nx1, nx1), ddx_bwd (nx1, nx1)
  {
    m_nx1 = nx1;
    m_x_disc_meth = x_disc_meth;
  };
  void fillDiscreteMatrix ();
  void fillAbsVal ();

  void saveOutput (int iter, double t, std::string fpref);
};

void
Mesh::fillAbsVal ()
{
  /*
  Fills absolute value coefficients for upwind scheme
  */
  for (int i = 0; i < m_nx1; ++i)
    {
      umabsu (i) = (u (i) - abs (u (i))) / 2;
      upabsu (i) = (u (i) + abs (u (i))) / 2;
    }
}

void
Mesh::saveOutput (int iter, double t, std::string fpref)
{
  /*
  Saves simulation output to csv file specified in main.cpp
  */

  std::string fname = getFname (iter, fpref);

  std::ofstream file (fname);
  for (int i = 0; i < m_nx1; ++i)
    {
      file << u (i) << ", ";
    }

  file << t;
}

void
Mesh::fillDiscreteMatrix ()
{
  /*
  Fills derivative operators in the mesh class with periodic
  boundary conditions
  */

  for (int i = 0; i < m_nx1; ++i)
    {
      for (int j = 0; j < m_nx1; ++j)
        {

          // Fill second order derivative spatial T matrix
          // Fill first order derivative forward difference
          // Fill first order derivative backward diff
          if (i == j)
            {
              ddx2_op (i, j) = -2;
              if (m_x_disc_meth == "upwind")
                {
                  ddx_fwd (i, j) = -1;
                  ddx_bwd (i, j) = 1;
                }
              else
                {
                  ddx_ctr (i, j) = 0;
                }
            }
          if (i + 1 == j)
            {
              ddx2_op (i, j) = 1;
              if (m_x_disc_meth == "upwind")
                {
                  ddx_fwd (i, j) = 1;
                  ddx_bwd (i, j) = 0;
                }
              else
                {
                  ddx_ctr (i, j) = 1;
                }
            }
          if (i - 1 == j)
            {
              ddx2_op (i, j) = 1;
              if (m_x_disc_meth == "upwind")
                {
                  ddx_fwd (i, j) = 0;
                  ddx_bwd (i, j) = -1;
                }
              else
                {
                  ddx_ctr (i, j) = -1;
                }
            }
        }
    }

  // Periodic Boundary Conditions
  ddx2_op (0, m_nx1 - 1) = 1;
  ddx2_op (m_nx1 - 1, 0) = 1;

  ddx_bwd (0, m_nx1 - 1) = 1;
  ddx_bwd (m_nx1 - 1, 0) = 1;

  ddx_fwd (0, m_nx1 - 1) = 1;
  ddx_fwd (m_nx1 - 1, 0) = 1;

  ddx_ctr (0, m_nx1 - 1) = 1;
  ddx_ctr (m_nx1 - 1, 0) = 1;
}

#endif