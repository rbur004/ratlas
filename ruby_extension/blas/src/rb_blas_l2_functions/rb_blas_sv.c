/*Ruby will pass in 
  The Blas class (self) already holds the matrix datatype and dimensions, so these don't need to be passed in.
  the first argument is an array of options argument (or nil) for the UPLO, TRANS and DIAG options
  the second and third arguments are the vector X and incX

  The function to use is selected from the matrix data_type as an index into sv.
 
    xTRSV, xTBSV and xTPSV  performs one of the matrix-vector operations   
         x := A.inverse*x
    or   x := A.transpose.inverse*x  
    or   x := A.conjugate_transpose.inverse*x 
    where x is an n element vector and  A is an n by n unit, or non-unit,   
    upper or lower triangular matrix.   
*/


#include "rb_blas.h"

VALUE rb_blas_sv(int argc, VALUE *argv, VALUE self)
{
  return self;
}
