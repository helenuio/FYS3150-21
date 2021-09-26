#include "jacobi_rotational.hpp"
#include <iostream>
#include <armadillo>
#include <cmath>

int main(int argc, char const *argv[]) {
  int n = 6;
  double h = 1.0/(n);
  double hh = h*h;
  double a = -1/hh;
  double d = 2/hh;
  double eps = 1.0E-9;


  arma::vec eigval;
  arma::mat eigvec;

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

  //std::cout << "analytical eigenvectors" << "\n" << normalise(eigvec_analytical)*(-1) << std::endl;
  std::cout << "analytical eigenvalues" << "\n" << eigval_analytical << std::endl;

  Jacobi_rotate solver;
  solver.initialize(a, d, eps, n);
  arma::mat A = solver.create_tridiagonal();
  std::cout << "A before" << "\n" << A << std::endl;

  eig_sym(eigval, eigvec, A);
  //std::cout << A << std::endl;
  std::cout << "Using eig_sym eigval from Armadillo" << "\n" << eigval << std::endl;
  //std::cout << "Using eig_sym eigvec from Armadillo" << "\n" << eigvec << std::endl;
  solver.jacobi_eigensolver(A);


/*
  // Test create_tridiagonal-function using a 4x4 matrix
  // Create a real test later
  Jacobi_rotate tester;
  tester.initialize(0, 1, eps, 4);
  arma::mat C = tester.create_tridiagonal();
  C(0,3) = 0.5; C(3,0) = 0.5; C(1,2) = -0.7; C(2,1) = -0.7;
  double max_C = tester.max_offdiag_symmetric(C);
  std::cout << C << std::endl;
  std::cout << max_C << std::endl;
*/
  return 0;
}
