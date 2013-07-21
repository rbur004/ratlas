/*
      INTEGER FUNCTION IXAMAX(N,X,INCX)
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
*     finds the index of element having max. absolute value.

*   Ruby
*     x.max(n = x.ncols, incx = 1)
*/


#include "rb_blas.h"

VALUE rb_blas_ixamax(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx;
  int incx;
  int incy;
  int n;
  //char error_msg[64];
  VALUE n_value,  incx_value;
  
  rb_scan_args(argc, argv, "02",  &incx_value,  &n_value);
  
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
  { //sprintf(error_msg, );
    rb_raise(rb_eRuntimeError, "Self is not a Vector");
  }
  
  switch(dx->data_type)
  {
  case Single_t: //s
    return INT2FIX(cblas_isamax(n , (float *)dx->data, incx)); 
  case Double_t: //d
    return INT2FIX(cblas_idamax(n , (double *)dx->data, incx)); 
  case Complex_t: //c
    return INT2FIX(cblas_icamax(n , dx->data, incx)); 
  case Double_Complex_t: //z
    return INT2FIX(cblas_izamax(n , dx->data, incx)); 
  default:
    //sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, "Invalid data_type (%d) in Matrix", dx->data_type);
    return Qnil; //Never reaches here.
  }
}

