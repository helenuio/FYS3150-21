#ifndef JACOBI_ROTATIONAL.HPP
#define JACOBI_ROTATIONAL.HPP

#include <armadillo>

class Jacobi_rotate {
private:
  double m_a, m_d, m_eps;
  int m_n, m_k, m_l;
public:
  void initialize(double a, double d, double eps, int n);
  arma::mat create_tridiagonal();
  double max_offdiag_symmetric(const arma::mat& A);
  void rotate(arma::mat& A, arma::mat& R);
  void jacobi_eigensolver(arma::mat& A);
  


};

#endif
