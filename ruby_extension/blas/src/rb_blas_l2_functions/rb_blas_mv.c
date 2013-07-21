
/*Ruby will pass in 
  The Blas class (self) already holds the matrix datatype and dimensions, so these don't need to be passed in.
  the first argument is an array of options argument (or nil) for the UPLO, TRANS and DIAG options
  The second argument is the scalar alpha.
  the third and forth arguments are the vector X and incX
  the fourth, fifth and sixth arguments are optional, and are the scalar beta, the vector Y and incY 

  The function to use is selected from the matrix data_type as an index into mv.
  The Vector Y and its associated arguments are not used for operations on TB, TP and TR matrices.

	Purpose   
    =======
    xGEMV and xGBMV performs one of the matrix-vector operations
       y = alpha*A*x + Beta*y 
    or y = alpha*A.transpose*x + Beta*y 
    or y = alpha*A.conjugate_transpose*x + Beta*y
    where x and y are n element vectors, alpha and beta are scalars, and A is a general matrix.
    
    xHEMV, xHBMV, xHPMV, xSYMV, xSBMV and xSPMV performs one of the matrix-vector operations
       y = alpha*A*x + Beta*y
    where x and y are n element vectors, alpha and beta are scalars, and A is a general matrix.
 
    xTRMV, xTBMV and xTPMV  performs one of the matrix-vector operations   
         x := A*x
    or   x := A.transpose*x  
    or   x := A.conjugate_transpose*x 
    where x is an n element vector and  A is an n by n unit, or non-unit,   
    upper or lower triangular matrix.   
*/

#include "rb_blas.h"

VALUE rb_blas_mv(int argc, VALUE *argv, VALUE self)
{
  return self;
}