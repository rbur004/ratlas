/*
      REAL or Double FUNCTION XNRM2(N, X, INCX)
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
*  XNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     SNRM2 := sqrt( x'*x ). or DZNRM2 := sqrt( conjg( x' )*x ) for complex numbers

*   RUBY
*     x.nrm2(n = x.ncols, incx = 1)
*/

#include "rb_blas.h"

VALUE rb_blas_xnrm2(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx;
  int incx;
  int incy;
  int n;
  //char error_msg[64];
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
  { //sprintf(error_msg, "Self is not a Vector");
    rb_raise(rb_eRuntimeError, "Self is not a Vector");
  }
  
  switch(dx->data_type)
  {
  case Single_t: //s
    return rb_float_new(cblas_snrm2(n , (float *)dx->data, incx)); 
  case Double_t: //d
    return rb_float_new(cblas_dnrm2(n , (double *)dx->data, incx)); 
  case Complex_t: //c
    return rb_float_new(cblas_scnrm2(n , dx->data, incx)); 
  case Double_Complex_t: //z
    return rb_float_new(cblas_dznrm2(n , dx->data, incx)); 
  default:
    //sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, "Invalid data_type (%d) in Matrix", dx->data_type);
    return Qnil; //Never reaches here.
  }
}
