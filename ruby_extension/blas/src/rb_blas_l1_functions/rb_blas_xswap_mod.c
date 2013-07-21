
/*
      SUBROUTINE XSWAP(N,X,INCX,Y,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      MATRIX X(*),Y(*)
*     ..
*
*  Purpose
*  =======
*
*     interchanges two vectors.
*     uses unrolled loops for increments equal to 1.

*   Ruby
*     x.swap!(y, n = x.ncols, incx = 1, incy = 1)
*/

#include "rb_blas.h"

VALUE rb_blas_xswap_mod(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx, *dy;
  int incx;
  int incy;
  int n;
  //char error_msg[64];
  VALUE vector_y,  n_value,  incx_value,  incy_value;
  
  rb_scan_args(argc, argv, "13", &vector_y, &incx_value, &incy_value, &n_value);
  
  Data_Get_Struct(self, Matrix, dx);
  Data_Get_Struct(vector_y, Matrix, dy);

  if(incx_value == Qnil)
    incx = 1;
  else
    incx = NUM2INT(incx_value);
  
  if(incy_value == Qnil)
    incy = 1;
  else
    incy = NUM2INT(incy_value);

  if(n_value == Qnil)
    n = dx->nrows;
  else
    n = NUM2INT(n_value);

  if(dx == NULL || dx->ncols != 1)
  { //sprintf(error_msg, "Self is not a Vector");
    rb_raise(rb_eRuntimeError, "Self is not a Vector");
  }
  
  if(dy == NULL || dy->ncols != 1)
  { //sprintf(error_msg, "Argument dy is not a Vector");
    rb_raise(rb_eRuntimeError, "Argument dy is not a Vector");
  }
  
  if(dy->data_type != dx->data_type)
  { //sprintf(error_msg, "Vectors are different data_types");
    rb_raise(rb_eRuntimeError, "Vectors are different data_types");
  }
    
  switch(dx->data_type)
  {
  case Single_t: //s
    cblas_sswap(n , (float *)dx->data, incx, (float *)dy->data, incy ); 
    break;
  case Double_t: //d
    cblas_dswap(n , (double *)dx->data, incx, (double *)dy->data, incy ); 
    break;
  case Complex_t: //c
    cblas_cswap(n , dx->data, incx, dy->data, incy ); 
    break;
  case Double_Complex_t: //z
    cblas_zswap(n , dx->data, incx, dy->data, incy ); 
    break;
  default:
    //sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, "Invalid data_type (%d) in Matrix", dx->data_type);
    break; //Never reaches here.
  }

  return self;
}
