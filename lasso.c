// lasso.c - a LASSO solver
//
// by Artem Sokolov

#include <R.h>
#include <math.h>

// Soft threshold
void sthresh( double* x, double a )
{
  if( fabs(*x) < a ) *x = 0;
  else if( *x < 0 ) *x = *x + a;
  else *x = *x - a;
}

// Updates fits and L2-norm penalties
void updateFits( double* X, double* S, int* np, int* jp, double* wj_diffp )
{
  int n = *np;
  int j = *jp;
  double wjd = *wj_diffp;

  for( int i = 0; i < n; ++i )
    S[i] += X[i+j*n] * wjd;
}

// Computes the objective function value
void objVal( double* S, double* z, double* w,
	     double* lambdap, int* np, int* pp, double* res )
{
  int n = *np;
  int p = *pp;
  double lambda = *lambdap;

  // Loss term
  double loss = 0.0;
  for( int i = 0; i < n; ++i )
    {
      double r = (z[i] - S[i]);
      loss += r * r;
    }

  // Regularization term
  double regL1 = 0.0;
  for( int j = 0; j < p; ++j )
    regL1 += fabs( w[j] );

  *res = 0.5 * loss / n + lambda * regL1;
}

// Computes the new value for coordinate *jp
void computeCoord( double* X, double* z, double* lambdap, double* S,
		   int* np, int* pp, int* jp, double* w,
		   double* work_zj, double* res )
{
  // Dereference
  int n = *np; int p = *pp; int j = *jp;
  double lambda = *lambdap;

  // Compute the working space values
  for( int i = 0; i < n; ++i )
    work_zj[i] = S[i] - X[i+j*n] * w[j];

  // Compute the numerator
  double num = 0.0;
  for( int i = 0; i < n; ++i )
    num += X[i+j*n] * (z[i] - work_zj[i]);

  // Normalize the numerator
  num /= n;
  sthresh( &num, lambda );
  if( num == 0.0 ) { *res = 0.0; return; }

  // Compute the denominator
  double denom = 0.0;
  for( int i = 0; i < n; ++i )
    denom += X[i+j*n] * X[i+j*n];

  // Normalize the denominator
  denom /= n;

  *res = num / denom;
}

// Optimizes the a LASSO objective via coordinate descent
void optLASSO( double* X, double* z, double* lambdap,
	       double* S, int* np, int* pp,
	       int* max_iter, double* eps,
	       double* w, double* b, int* bSilentp )
{
  // Dereference
  int n = *np; int p = *pp;
  double lambda = *lambdap;

  // Working storage
  double* work_zj = (double*) R_alloc( n, sizeof( double ) );

  if( !(*bSilentp) )
    Rprintf( "Running base optimization with lambda = %f\n", lambda );

  // Compute the initial objective function value
  double fprev;
  objVal( S, z, w, &lambda, np, pp, &fprev );

  //  Rprintf( "nTop = %d, obj = %f\n", nTop, fprev );

  // Perform coordinate descent
  int iter; double f;
  for( iter = 1; iter <= (*max_iter); ++iter )
    {
      //      Rprintf( "\n=== Iteration %d ===\n", iter );

      // Update the weights
      for( int j = 0; j < p; ++j )
	{
	  //	  Rprintf( "== Coord %d ==\n", j );

	  // Perform the update
	  double wj_old = w[j];
	  computeCoord( X, z, &lambda, S, np, pp, &j, w, work_zj, w+j );

	  //	  Rprintf( "computeCoord returned %f\n", w[j] );

	  // Update fits and L2-norm penalty term accordingly
	  double wj_diff = w[j] - wj_old;
	  if( wj_diff != 0.0 )
	    updateFits( X, S, np, &j, &wj_diff);
	}

      // Update the bias term
      double b_num = 0.0;
      double b_denom = n;
      for( int i = 0; i < n; ++i )
	{
	  double s = S[i] - *b;
	  b_num += (z[i] - s);
	}
      
      double b_old = *b;
      *b = b_num / b_denom;
      double b_diff = *b - b_old;

      // Update the fits accordingly
      if( b_diff != 0 )
	{
	  for( int i = 0; i < n; ++i )
	    S[i] += b_diff;
	}

      // Compute the objective function value and check the stopping criterion
      objVal( S, z, w, &lambda, np, pp, &f );
      //      Rprintf( "f = %f\n", f );
      if( fabs( f - fprev ) / fabs( fprev ) < *eps ) break;
      else fprev = f;
    }
  if( iter > (*max_iter) ) --iter;	// Corner case: loop didn't end via break
  if( !(*bSilentp) )
    Rprintf( "f = %f after iteration %d\n", f, iter );
}

