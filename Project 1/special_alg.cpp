#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

// object for output files
ofstream ofile;

// f(x) -> right side of the equation
inline double f(double x)
{
  return 100.0*exp(-10.0*x);
}

// exact solution
inline double exact(double x)
{
  return 1.0-(1-exp(-10))*x-exp(-10*x);
}

int main(int argc, char* argv[])
{
  int exponent; // 10^exponent
  string filename; // name of file taken from command line to read
  if( argc <= 1)
  {
    cout << "Please write name of file to read from as well as the number of grid points n" << endl;
    exit(1); // terminate with error
  }
    else
    {
      filename = argv[1]; // argv[0] is the name of the program
      exponent = atoi(argv[2]); // convert from string to integer
    }
  clock_t start, finish; // declare start and final time
  start = clock();
  for (int i = 1; i <= exponent; i++) // n=10, n=100, n=3 if exponent = 3
  {
      int n = (int) pow(10.0,i);
      double* x = new double[n+1];        x[0] = 0.0;
      double* u = new double[n+1];
      double* g = new double[n+1];        g[0] = 0.0;
      double* d = new double[n+1];        d[0] = 0.0;

      double h = 1.0/(n);
      double hh = h*h;

      for (int i = 1; i <= n+1; i++)
      {
        x[i] = i*h;
        g[i] = hh*f(x[i]);
      }

    // boundary conditions
    u[0] = u[n+1] = 0.0;
    d[1] = 2;

    // forward substitution
    for (int i = 2; i <= n+1; i++)
    {
      d[i] = d[1] - 1/d[i-1];
      g[i] = g[i] + g[i-1]/d[i-1];
    }

    // backward substitution
    for (int i = n; i >= 1; i--)
    {
      u[i] = (g[i] + u[i+1])/d[i];
    }

    // filename for each exponent
    string fileout = filename;

    // convert exponent to string
    string exponent = to_string(i);

    // filename-i (where i = exponent)
    fileout.append(exponent);

    ofile.open(fileout);

    // write to file for easier plotting
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    //      ofile << "       x:             approx:          exact:       relative error" << endl;
    for (int i = 1; i <= n-1; i++)
    {
      double xval = x[i];
      double RelativeError = fabs((exact(xval)-u[i])/exact(xval));
      ofile << setw(15) << setprecision(8) << xval;
      ofile << setw(15) << setprecision(8) << u[i];
      ofile << setw(15) << setprecision(8) << exact(xval);
      ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
    }

    ofile.close();

    delete [] x; delete [] d;
    delete [] u; delete [] g;

    finish = clock();
    double timeused = (double) (finish - start)/(CLOCKS_PER_SEC );
    cout << setw(5) << n;
    cout << setw(15) << timeused << endl;
  } // closing function looping over n=0,1,2
return 0;
} // closing main function
