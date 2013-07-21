
/*DOUBLE PRECISION FUNCTION XDOT(N,X,INCX,Y,INCY)
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      MATRIX X(*), Y(*)
*     ..
*
*  Purpose
*  =======
*
*     forms the dot product of two vectors.
*     uses unrolled loops for increments equal to one.

*   Ruby
      x.xdot(y, n = x.ncols, incx = 1, incy = 1)
If complex matrix, calls dotu version.
*
*/

#include "rb_blas.h"

VALUE rb_blas_xdot(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx, *dy;
  int incx;
  int incy;
  int n;
  //char error_msg[64];
  VALUE vector_y,  n_value,  incx_value,  incy_value;
  float result_c[2];
  double result_z[2];
  VALUE complex_result[2];
  
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
  { //sprintf(error_msg, );
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
    return rb_float_new(cblas_sdot(n , (float *)dx->data, incx, (float *)dy->data, incy )); 
  case Double_t: //d
    return rb_float_new(cblas_ddot(n , (double *)dx->data, incx, (double *)dy->data, incy )); 
  case Complex_t: //c
    cblas_cdotu_sub(n , dx->data, incx, dy->data, incy, result_c ); 
    complex_result[0] = rb_float_new(result_c[0]);
    complex_result[1] = rb_float_new(result_c[1]);
    return rb_class_new_instance(2, complex_result, rb_intern("Complex"));
  case Double_Complex_t: //z
    cblas_zdotu_sub(n , dx->data, incx, dy->data, incy, result_z ); 
    complex_result[0] = rb_float_new(result_z[0]);
    complex_result[1] = rb_float_new(result_z[1]);
    return rb_class_new_instance(2, complex_result, rb_intern("Complex"));
  default:
    //sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, "Invalid data_type (%d) in Matrix", dx->data_type);
    return Qnil; //Never reaches here.
  }
}
