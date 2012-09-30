/*Ruby will pass in 
  The Blas class (self) already holds the matrix datatype and dimensions, so these don't need to be passed in.
  the first argument is an array of options argument (or nil) for the UPLO, TRANS and DIAG options
  the second and third arguments are the vector X and incX
  the forth and fifth arguments are the vector Y and incY which aren't used by H and S matrices

  The function to use is selected from the matrix data_type as an index into sv.
 
    xGER, xGERU  performs one of the matrix-vector operations   
         x := alpha*X*Y.transpose + A
    where x and y are n element vectors 

    xHER, xHPR performs one of the matrix-vector operations   
         x := alpha*X*X.conjugate_transpose + A

    xSPR, xSYR performs one of the matrix-vector operations   
         x := alpha*X*X.transpose + A
*/


static VALUE rb_blas_r(int argc, VALUE *argv, VALUE self)
{
  return self;
}