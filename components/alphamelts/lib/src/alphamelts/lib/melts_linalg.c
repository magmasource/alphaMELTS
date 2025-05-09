const char *melts_linalg_ver(void) { return "$Id: melts_linalg_ver.c,v 1.0 2017/03/03 18:24:06 psmith Exp $"; }

/*
**++
**  FACILITY:  Silicate Melts Regression/Crystallization Package
**
**  MODULE DESCRIPTION:
**
**      Wrappers for the Linear Algebra moduls of the GNU Scientific Library
**
**      https://www.gnu.org/software/gsl/
**      https://www.gnu.org/software/gsl/manual/html_node/
**
**      File: MELTS_LINALG.C
**
**  MODIFICATION HISTORY:
**      V1.0-1  March 03, 2017
**
**              Pack arrays and perform LU decomp etc. only for non-singular part
**--
*/

#include "melts_gsl.h"

int melts_LU_decomp (gsl_matrix * A, gsl_permutation * p, int *signum)
{

  /* Assumes that the matrix A is symmetrical. Any zero rows / columns are
     ignored for LU decomposition and later solution. If the packed matrix,
     A2, is 1x1 LU decomposition does nothing so 1/A is stored. If A2 is
     0x0 (singular) then store 0.0 instead. */

  if (A->size1 != A->size2)
    {
      GSL_ERROR ("LU decomposition requires square matrix", GSL_ENOTSQR);
    }
  else if (p->size != A->size1)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = A->size1;
      size_t i, j, N2;
      gsl_vector_view row;

      gsl_permutation_init (p);

      for (j = 0, N2 = N; j < N; j++)
        {
          row = gsl_matrix_row(A, j);
          if (gsl_vector_isnull(&row.vector))
            {
              N2--;
              for (i = j + 1; i < N; i++)
                {
                  row = gsl_matrix_row(A, i);
                  if (gsl_vector_isnull(&row.vector))
                    {
                      N2--;
                    }
                  else
                    {
                      gsl_matrix_swap_rows(A, j, i);
                      gsl_permutation_swap(p, j, i);
                      j = i;
                      break;
                    }
                }
            }
        }

      if (N2 <= 1)
        {
          size_t p0 = gsl_permutation_get(p, 0);
          double A0 = gsl_matrix_get(A, 0, p0); /* columns were not swapped */

          gsl_matrix_set(A, 0, p0, 0.0);
          if (A0 != 0.0) gsl_matrix_set(A, 0, 0, 1.0/A0);

        }
      else if (N2 < N)
        {
          gsl_matrix_view A2 = gsl_matrix_submatrix(A, 0, 0, N2, N2);
          gsl_permutation *p2 = gsl_permutation_alloc(N2);
          size_t i;

          // In version 2.2
          //      gsl_permute_matrix(p, A); /* swap non-zero columns into A2 */
          for (i = 0; i < N; i++)
            {
              gsl_vector_view row = gsl_matrix_row(A, i);
              gsl_permute_vector(p, &row.vector);
            }

          gsl_linalg_LU_decomp(&A2.matrix, p2, signum);

          for (i = 0; i < N2; i++)
            {
              gsl_permutation_swap(p, i, p2->data[i]);
            }
          gsl_permutation_free(p2);

        }
      else
        {
          gsl_linalg_LU_decomp(A, p, signum);
        }

  return GSL_SUCCESS;
  }

}

int
melts_LU_solve (const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x)
{
  if (b->size != x->size)
    {
      GSL_ERROR ("b size must match solution size", GSL_EBADLEN);
    }
  else
    {
      int status;

      /* see gsl_linalg_LU_solve */
      gsl_vector_memcpy (x, b);
      status = melts_LU_svx (LU, p, x);

      return status;
    }
}

int
melts_LU_svx (const gsl_matrix * LU, const gsl_permutation * p, gsl_vector * x)
{
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution/rhs size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = LU->size1;
      size_t j, N2;

      for (j = 0, N2 = 0; j < N; j++)
        {
          gsl_vector_const_view row = gsl_matrix_const_row(LU, j);
          if (!gsl_vector_isnull(&row.vector))
            {
              N2++;
            }
        }

      if (N2 <= 1)
        {
          size_t p0 = gsl_permutation_get(p, 0);
          double x0 = gsl_vector_get(x, p0);

          gsl_vector_set_zero(x);
          gsl_vector_set(x, p0, x0*gsl_matrix_get(LU, 0, 0));

        }
      else if (N2 < N)
        {
          gsl_matrix_const_view LU2 = gsl_matrix_const_submatrix(LU, 0, 0, N2, N2);
          gsl_vector_view x2 = gsl_vector_subvector(x, 0, N2);

          gsl_permute_vector(p, x);

          /* see gsl_linalg_LU_svx */
          gsl_blas_dtrsv (CblasLower, CblasNoTrans, CblasUnit, &LU2.matrix, &x2.vector);
          gsl_blas_dtrsv (CblasUpper, CblasNoTrans, CblasNonUnit, &LU2.matrix, &x2.vector);

          for (j = N2; j < N; j++)
            {
              gsl_vector_set(x, j, 0.0);
            }
            gsl_permute_vector_inverse(p, x);
        }
      else
        {
          gsl_linalg_LU_svx (LU, p, x);
        }

      return GSL_SUCCESS;
  }
}

int
melts_LU_refine (const gsl_matrix * A, const gsl_matrix * LU, const gsl_permutation * p, const gsl_vector * b, gsl_vector * x, gsl_vector * residual)
{
  /* Assumes that the matrix A is symmetrical. Any zero rows / columns are
     ignored for LU decomposition and later solution. If the packed matrix
     in melts_LU_decomp, A2, is 1x1 LU decomposition does nothing so 1/A
     was stored. If A2 was 0x0 (singular) then 0.0 was stored. */
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("matrix a must be square", GSL_ENOTSQR);
    }
  if (LU->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be square", GSL_ENOTSQR);
    }
  else if (A->size1 != LU->size2)
    {
      GSL_ERROR ("LU matrix must be decomposition of a", GSL_ENOTSQR);
    }
  else if (LU->size1 != p->size)
    {
      GSL_ERROR ("permutation length must match matrix size", GSL_EBADLEN);
    }
  else if (LU->size1 != b->size)
    {
      GSL_ERROR ("matrix size must match b size", GSL_EBADLEN);
    }
  else if (LU->size1 != x->size)
    {
      GSL_ERROR ("matrix size must match solution size", GSL_EBADLEN);
    }
  else
    {
      const size_t N = LU->size1;
      size_t j, N2;

      for (j = 0, N2 = 0; j < N; j++)
        {
          gsl_vector_const_view row = gsl_matrix_const_row(LU, j);
          if (!gsl_vector_isnull(&row.vector))
            {
              N2++;
            }
        }

      int status;

      if (N2 <= 1)
        {
          status = GSL_SUCCESS; /* don't do anything */
        }
      else if (N2 < N)
        {
          /* see gsl_linalg_LU_refine */
          gsl_vector_memcpy (residual, b);
          gsl_blas_dgemv (CblasNoTrans, 1.0, A, x, -1.0, residual);

          status = melts_LU_svx (LU, p, residual);
          gsl_blas_daxpy (-1.0, residual, x);
        }
      else
        {
          status = gsl_linalg_LU_refine(A, LU, p, b, x, residual);
        }

      return status;
  }
}
