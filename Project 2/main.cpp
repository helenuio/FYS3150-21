#include "jacobi_rotational.hpp"
#include <iostream>
#include <armadillo>
#include <cmath>

int main(int argc, char const *argv[]) {
  int n = 80;
  double h = 1.0/(n);
  double hh = h*h;
  double a = -1/hh;
  double d = 2/hh;
  double eps = 1.0E-9;

  // Analytical solutions for eigenvectors and eigenvalues
  arma::vec eigval_analytical = arma::vec(n);
  arma::mat eigvec_analytical(n, n);

  for (int i=1; i <= n; i++)
  {
    eigval_analytical(i-1) = d + 2*a*cos(i*M_PI/(n+1));
    for (int j=1; j <= n; j++)
    {
      eigvec_analytical(i-1,j-1) = sin(i*j*M_PI/(n+1));
    }
  }

  /*
  std::cout << std::setprecision(4) << "analytical eigenvalues" << "\n" << eigval_analytical << std::endl;
  std::cout << std::setprecision(4) << "analytical eigenvectors" << "\n" << normalise(eigvec_analytical)*(-1) << std::endl;
  */

  // Creating matrix A
  Jacobi_rotate solver;
  solver.initialize(a, d, eps, n);
  arma::mat A = solver.create_tridiagonal();

  // Printing eigval, eigvec using Armadillo
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, A);
  /*
  std::cout << "Using eigval from Armadillo" << "\n" << eigval << std::endl;
  std::cout << "Using eigvec from Armadillo" << "\n" << eigvec << std::endl;
  */

  // Numerical solution
  // Prints numerical eigenvalues and eigenvalues and
  // writes the three eigenvectors corresponding to the
  // three lowest eigenvalues to a file
  solver.jacobi_eigensolver(A);

  /*
  // Testing create_tridiagonal-function using a 4x4 matrix "C"
  Jacobi_rotate tester;
  tester.initialize(0, 1, eps, 4);
  arma::mat C = tester.create_tridiagonal();
  C(0,3) = 0.5; C(3,0) = 0.5; C(1,2) = -0.7; C(2,1) = -0.7;
  double expected_value = 0.7;
  double max_C = tester.max_offdiag_symmetric(C);
  //std::cout << C << std::endl;
  if (fabs(max_C - expected_value) < eps) {
    std::cout << max_C << " is indeed " << expected_value << std::endl;
  }
  */
  return 0;
}
