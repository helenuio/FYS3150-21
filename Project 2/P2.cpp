#include <armadillo>
#include <iostream>

// Create a tridiagonal matrix tridiag(a,d,e) of size n*n,
// from scalar input a, d, and e. That is, create a matrix where
// - all n-1 elements on the subdiagonal have value a
// - all n elements on the diagonal have value d
// - all n-1 elements on the superdiagonal have value e
arma::mat create_tridiagonal(int n, double a, double d)
{
  // Start from identity matrix
  arma::mat A = arma::mat(n, n, arma::fill::eye);

  // Fill the first row (row index 0), e.g.
  A(0,0) = d;
  A(0,1) = a;

  // Loop that fills rows 2 to n-1 (row indices 1 to n-2)
  for (int i=1; i < n-1; i++)
  {
    A(i+1,i) = a;
    A(i,i+1) = a;
    A(i,i) = d;
  }

  // Fill last row (row index n-1)
  A(n-1,n-1) = d;
  A(n-1, n-2) = a;

  return A;
}

// Create a symmetric tridiagonal matrix tridiag(a,d,a) of size n*n
// from scalar input a and d.
arma::mat create_symmetric_tridiagonal(int n, double a, double d)
{
  // Call create_tridiagonal and return the result
  return create_tridiagonal(n, a, d);
}

int main(int argc, char* argv[])
{
  if( argc <= 3)
  {
    std::cout << "argh" << std::endl;
    exit(1); // terminate with error
  }
    else
    {
      int n = atoi(argv[1]);
      double a = atoi(argv[2]);
      double d = atoi(argv[3]);
      arma::mat A = create_symmetric_tridiagonal(n, a, d);

      std::cout << A << std::endl;

      arma::vec eigval;
      arma::mat eigvec;

      eig_sym(eigval, eigvec, A);
      std::cout << eigval << std::endl;
      std::cout << normalise(eigvec, 1) << std::endl;

    }
}
