#include "jacobi_rotational.hpp"
#include <iostream>
#include <cmath>

void Jacobi_rotate::initialize(double a, double d, double eps, int n)
{
  m_a = a;
  m_d = d;
  m_eps = eps;
  m_n = n;
  int m_k;
  int m_l;
}

// Create a tridiagonal matrix tridiag(a,d,e) of size n*n
arma::mat Jacobi_rotate::create_tridiagonal()
{
  // Start from identity matrix
  arma::mat A = arma::mat(m_n, m_n, arma::fill::eye);

  // Fill the first row (row index 0), e.g.
  A(0,0) = m_d;
  A(0,1) = m_a;

  // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
  for (int i=1; i < m_n-1; i++)
  {
    A(i+1,i) = m_a;
    A(i,i+1) = m_a;
    A(i,i) = m_d;
  }

  // Fill last row (row index n-1)
  A(m_n-1,m_n-1) = m_d;
  A(m_n-1, m_n-2) = m_a;

  return A;
}

// Find the largest off-diagonal element
double Jacobi_rotate::max_offdiag_symmetric(const arma::mat& A)
{
  double maxval = 0.0;
  for (int i = 0; i < m_n; i++) {
    for (int j = i+1; j < m_n; j++) {
      if (fabs(A(i,j)) > maxval) {
        maxval = fabs(A(i,j));
        m_k = i;
        m_l = j;
      }
    }
  }
  return maxval;
}

// Performs a single Jacobi rotation
void Jacobi_rotate::rotate(arma::mat& A, arma::mat& R)
{

  // Compute sin(theta), cos(theta) and tan(theta)
  double s, c;
  if (A(m_k,m_l) == 0.0) {
    s = 0.0;
    c = 1.0;
  }
  else {
    double t, tau;
    tau = (A(m_l,m_l) - A(m_k,m_k))/(2*A(m_k,m_l));
    if (tau > 0.0) {
      t = -tau + sqrt(1 + (tau*tau));
    }
    else {
      t = -tau - sqrt(1 + (tau*tau));
    }
    c = 1/sqrt(1 + (t*t));
    s = t*c;
  }


  // Modify matrix elements with indices k and l
  double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
  a_kk = A(m_k,m_k);
  a_ll = A(m_l,m_l);
  A(m_k,m_k) = a_kk*c*c - 2*A(m_k,m_l)*c*s + a_ll*s*s;
  A(m_l,m_l) = a_ll*c*c + 2*A(m_k,m_l)*c*s + a_kk*s*s;
  //A(m_k,m_k) = A(m_k,m_k)*c*c - 2*A(m_k,m_l)*c*s + A(m_l,m_l)*s*s;
  //A(m_l,m_l) = A(m_l,m_l)*c*c + 2*A(m_k,m_l)*c*s + A(m_k,m_k)*s*s;
  A(m_k,m_l) = 0;
  A(m_l,m_k) = 0;

  // Update the rest of the elements in the matrix
  for (int i = 0; i < m_n; i++) {
    if (i != m_k && i != m_l) { // for all i != k,l
      a_ik = A(i,m_k); // Keep previous elements A(i,k) and A(i,l)
      a_il = A(i,m_l); // to use in next calculation
      A(i,m_k) = a_ik*c - a_il*s;
      A(m_k,i) = A(i,m_k);
      A(i,m_l) = a_il*c + a_ik*s;
      A(m_l,i) = A(i,m_l);
    }

    // Update the overall rotation matrix R
    r_ik = R(i,m_k);
    r_il = R(i,m_l);
    R(i,m_k) = r_ik*c - r_il*s;
    R(i,m_l) = r_il*c + r_ik*s;
  }
}

void Jacobi_rotate::jacobi_eigensolver(arma::mat& A)
{
  // R starts out as identity matrix
  arma::mat R = arma::mat(m_n, m_n, arma::fill::eye);
  bool converged;
  int iterations = 0;
  int maxiter = m_n*m_n*m_n;
  arma::mat eigenvectors(m_n, m_n);
  arma::vec eigenvalues = arma::vec(m_n);
  double max = max_offdiag_symmetric(A);

  while (fabs(max) > m_eps && iterations < maxiter) {
    //std::cout << "l " << m_l << " k " << m_k << std::endl;
    max = max_offdiag_symmetric(A);
    rotate(A,R);
    iterations++;
  }

  // Set the bool reference "converged" to true if convergence was reached before hitting maxiter
  if (iterations < maxiter) {
    converged = true;
  }
  else {
    converged = false;
  }

  // Write the eigenvalues as entries in the vector "eigenvalues"
  for (int i = 0; i < m_n; i++) {
    eigenvalues(i) = A(i,i);
  }

  std::cout << "numerical eigenvalues" << "\n" << eigenvalues << std::endl;
  std::cout << "Diagonalized A" << "\n" << A << std::endl;
  //std::cout << "numerical eigenvectors" << "\n" << normalise(R) << std::endl;

}
