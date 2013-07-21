#include "rb_lapack.h"

/*
Purpose
*  =======
*
*  XGESV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as
*     A = P * L * U,
*  where P is a permutation matrix, L is unit lower triangular, and U is
*  upper triangular.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*Fortran SUBROUTINE XGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
  *  N       (input) INTEGER
  *          The number of linear equations, i.e., the order of the
  *          matrix A.  N >= 0.
  *
  *  NRHS    (input) INTEGER
  *          The number of right hand sides, i.e., the number of columns
  *          of the matrix B.  NRHS >= 0. You can have multiple solutions
  *          and not just a single vector for B.
  *
  *  A       (input/output) REAL array, dimension (LDA,N)
  *          On entry, the N-by-N coefficient matrix A.
  *          On exit, the factors L and U from the factorization
  *          A = P*L*U; the unit diagonal elements of L are not stored.
  *
  *  LDA     (input) INTEGER
  *          The leading dimension of the array A.  LDA >= max(1,N) for row major.
  *          or LDA >= max(1,M) for column major. (this allows non-contiguous arrays)
  *
  *  IPIV    (output) INTEGER array, dimension (N)
  *          The pivot indices that define the permutation matrix P;
  *          row i of the matrix was interchanged with row IPIV(i).
  *
  *  B       (input/output) REAL array, dimension (LDB,NRHS)
  *          On entry, the N-by-NRHS matrix of right hand side matrix B.
  *          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
  *
  *  LDB     (input) INTEGER
  *          The leading dimension of the array B.  LDB >= max(1,N).
  *
  *  INFO    (output) INTEGER
  *          = 0:  successful exit
  *          < 0:  if INFO = -i, the i-th argument had an illegal value
  *          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
  *                has been completed, but the factor U is exactly
  *                singular, so the solution could not be computed.
* Ruby method
  * A is an instance of the class
  * B is the result argument
  * IPIV is the ipiv argument
  * The other values, N, NRHS, LDA, LDB are inferred from the objects A & B
*/

VALUE rb_lapack_xgesv_mod(VALUE self, VALUE result, VALUE ipiv)
{
  Matrix *m, *r, *i;
  int error;
  int info;
  const char *args[] = { 
			  "Number of equations in 'A'",
              "Number of result columns in 'B'",
              "Input Array 'A'",
              "Leading Dimesion of 'A'",
              "Permutation Matrix 'IPIV'",
              "Result Array 'B'",
              "Leading Dimension of 'B'" };
  //char error_msg[64];
  
  Data_Get_Struct(self, Matrix, m);
  Data_Get_Struct(result, Matrix, r);
  Data_Get_Struct(ipiv, Matrix, i);
  
  //Storage is M rows x N Columns. Matrix is lda x N, which could be a sub-Matrix of MxN
  
  //int sgesv_(__CLPK_integer *n, __CLPK_integer *nrhs, __CLPK_real *a, __CLPK_integer *lda, 
  //        __CLPK_integer *ipiv, __CLPK_real *b, __CLPK_integer *ldb, __CLPK_integer *info);
  
  
  //SUBROUTINE XGESV( N, NRHS,   A,       LDA,        IPIV,     B,        LDB,    INFO )
  switch(m->data_type)
  {
  case Single_t: 
    error = sgesv_( &m->ncols, &r->ncols, m->data, &m->nrows, i->data, r->data, &r->nrows, &info);
    break;
  case Double_t:
    error = dgesv_( &m->ncols, &r->ncols, m->data, &m->nrows, i->data, r->data, &r->nrows, &info);
    break;
    
  case Complex_t:
    error = cgesv_(&m->ncols, &r->ncols, m->data, &m->nrows, i->data, r->data, &r->nrows, &info);
    break;
  case Double_Complex_t:
    error = zgesv_(&m->ncols, &r->ncols, m->data, &m->nrows, i->data, r->data, &r->nrows, &info);
    break;
    
  }

  if(info < 0)
  { //sprintf(error_msg, "Arg %s has illegal value", args[-info]);
    rb_raise(rb_eRuntimeError, "Arg %s has illegal value", args[-info]);
  }
  else if(info > 0)
  { //sprintf(error_msg, "Singular Matrix at %d", info);
    rb_raise(rb_eRuntimeError, "Singular Matrix at %d", info);
  }
  
  m->ipiv = i;
  m->matrix_type = LU;
  
  return result;
}

VALUE rb_lapack_xgesv( VALUE self, VALUE result, VALUE ipiv )
{
  Matrix *dx, *dy;
  VALUE self_copy, result_copy;
    
  Data_Get_Struct(self, Matrix, dx);
  self_copy = rblas_new_instance(dx->class_id, dx->data, dx->nrows, dx->ncols, dx->data_type, dx->matrix_type, dx->cblas_order);
  Data_Get_Struct(result, Matrix, dy);
  result_copy = rblas_new_instance(dy->class_id, dy->data, dy->nrows, dy->ncols, dy->data_type, dy->matrix_type, dy->cblas_order);
  return rb_lapack_xgesv_mod(self_copy, result_copy, ipiv);
}
