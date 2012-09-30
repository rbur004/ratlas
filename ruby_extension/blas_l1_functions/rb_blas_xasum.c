
/*
      REAL/DOUBLE FUNCTION XASUM(N, X, INCX)
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      MATRIX X(*)
*     ..
*
*  Purpose
*  =======
*
*     takes the sum of the absolute values.
*     uses unrolled loops for increment equal to one.

*   Ruby
*     x.asum(n = x.ncols, incx = 1)
*/

static VALUE rb_blas_xasum(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx;
  int incx;
  int incy;
  int n;
  char error_msg[64];
  VALUE n_value,  incx_value;
  
  rb_scan_args(argc, argv, "02",   &incx_value, &n_value);
  
  Data_Get_Struct(self, Matrix, dx);

  if(incx_value == Qnil)
    incx = 1;
  else
    incx = NUM2INT(incx_value);
  
  if(n_value == Qnil)
    n = dx->nrows;
  else
    n = NUM2INT(n_value);

  if(dx == NULL || dx->ncols != 1)
  { sprintf(error_msg, "Self is not a Vector");
    rb_raise(rb_eRuntimeError, error_msg);
  }
  
  switch(dx->data_type)
  {
  case Single_t: //s
    return rb_float_new(cblas_sasum(n , (float *)dx->data, incx)); 
  case Double_t: //d
    return rb_float_new(cblas_dasum(n , (double *)dx->data, incx)); 
  case Complex_t: //c
    return rb_float_new(cblas_scasum(n , dx->data, incx)); 
  case Double_Complex_t: //z
    return rb_float_new(cblas_dzasum(n , dx->data, incx)); 
  default:
    sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, error_msg);
    return Qnil; //Never reaches here.
  }
}
