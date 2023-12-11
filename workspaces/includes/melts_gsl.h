#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permute_vector.h>
//#include <gsl/gsl_permute_matrix.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_linalg.h>

int melts_LU_decomp (gsl_matrix * A, gsl_permutation * p, int *signum);
int melts_LU_solve (const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x);
int melts_LU_svx (const gsl_matrix * LU, const gsl_permutation * p, gsl_vector * x);
int melts_LU_refine (const gsl_matrix * A, const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x, gsl_vector * residual);
