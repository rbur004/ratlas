/*      SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  SGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) REAL array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*/
static VALUE rb_lapack_xgetrf_mod(VALUE self, VALUE ipiv)
{
  Matrix *m,  *i;
  int error;
  int info;
  char *args[] = { "Number of columns in A"
              "Number of rows in 'A'",
              "Input Array 'A' / Output LU Array",
              "Leading Dimesion of 'A'",
              "Output Permutation Matrix 'IPIV'"
              };
  char error_msg[64];
  
  Data_Get_Struct(self, Matrix, m);
  Data_Get_Struct(ipiv, Matrix, i);
  
///* Subroutine */ int cgetrf_(__CLPK_integer *m, __CLPK_integer *n, __CLPK_complex *a, __CLPK_integer *lda,
//         __CLPK_integer *ipiv, __CLPK_integer *info);
// Fortran           XGETRF( M,  N, A, LDA, IPIV, INFO)
  switch(m->data_type)
  {
  case Single_t: 
    error = sgetrf_(&m->ncols, &m->nrows, m->data, &m->ncols, i->data, &info); 
    break;
  case Double_t:
    error = dgetrf_(&m->ncols, &m->nrows, m->data, &m->ncols, i->data, &info);
    break;
    
  case Complex_t:
    error = cgetrf_(&m->ncols, &m->nrows, m->data, &m->ncols, i->data, &info);
    break;
  case Double_Complex_t:
    error = zgetrf_(&m->ncols, &m->nrows, m->data, &m->ncols, i->data, &info);
    break;
    
  }

  if(info < 0)
  { sprintf(error_msg, "Arg %s has illegal value", args[-info]);
    rb_raise(rb_eRuntimeError, error_msg);
  }
  else if(info > 0)
  { sprintf(error_msg, "Singular Matrix at %d", info);
    rb_raise(rb_eRuntimeError, error_msg);
  }
  
  m->ipiv = i;
  m->matrix_type = LU;
  
  return self;
}



