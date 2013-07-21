#include "rb_lapack.h"

/*
      SUBROUTINE SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  SGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by SGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) REAL array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by SGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from SGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) REAL array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*/



///* Subroutine */ int cgetrs_(char *trans, __CLPK_integer *n, __CLPK_integer *nrhs, __CLPK_complex *
//       a, __CLPK_integer *lda, __CLPK_integer *ipiv, __CLPK_complex *b, __CLPK_integer *ldb, __CLPK_integer *info);

VALUE rb_lapack_xgetrs_mod(VALUE self,  VALUE result, VALUE ipiv, VALUE transpose )
{
  Matrix *m,  *r, *i;
  int error;
  int info;
  char trans = 'N';
  const char *args[] = { 
              "Transpose argument error",
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
  
  switch(NUM2INT(transpose))
  {
    case CblasNoTrans: trans = 'N'; break;
    case CblasTrans: trans = 'T'; break;
    case CblasConjTrans: trans = 'C'; break;
  }
//   SUBROUTINE SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )

  switch(m->data_type)
  {
  case Single_t: 
    error = sgetrs_(&trans,  &m->ncols, &r->ncols, m->data, &m->nrows, i->data, r->data, &r->nrows, &info);
    break;
  case Double_t:
    error = dgetrs_(&trans,  &m->ncols, &r->ncols, m->data, &m->nrows, i->data, r->data, &r->nrows, &info);
    break;
    
  case Complex_t:
    error = cgetrs_(&trans,  &m->ncols, &r->ncols, m->data, &m->nrows, i->data, r->data, &r->nrows, &info);
    break;
  case Double_Complex_t:
    error = zgetrs_(&trans,  &m->ncols, &r->ncols, m->data, &m->nrows, i->data, r->data, &r->nrows, &info);
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
  
  return result;
}

VALUE rb_lapack_xgetrs( VALUE self, VALUE result, VALUE ipiv, VALUE transpose )
{
  Matrix *dy;
  VALUE  result_copy;
    
  Data_Get_Struct(result, Matrix, dy);
  result_copy = rblas_new_instance(dy->class_id, dy->data, dy->nrows, dy->ncols, dy->data_type, dy->matrix_type, dy->cblas_order);
  return rb_lapack_xgetrs_mod(self, result_copy, ipiv, transpose);
}
