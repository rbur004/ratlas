/*
      SUBROUTINE XAXPY(N, A, X,INCX, Y,INCY)
  and SUBROUTINE XAXPBY(N, A, X,INCX, Y,INCY) //Atlas extention that allows for A*X[] + B*Y[]
*     .. Scalar Arguments ..
      Constant Type_X  A is the scalar multiplier of each  X element, before the addition
      INTEGER INCX,INCY, are the increments used to step to the next  X[] and next Y[] (usually 1)
      Constant INTEGER N is the size of the two vectors (should be the same)
*     ..
*     .. Array Arguments ..
      Const Vector X[*],
      Vector Y[*] Warning: This call Modifies Y, using it for the result.
*     ..
*
*  Purpose
*  =======
*
*     A * X[] + B * Y[]
*     uses unrolled loops for increments equal to one.

*  Ruby interface varies slighty (alright a lot)
*   X is self, Y is the first argument
*  All the remaining arguments are optional,
*   A defaults to 1.0
*   B defaults to 1.0
*   N  defaults to the size of the DX Blas Vector
*   INCX & INCY default to 1
  i.e. 
    x.axpy(y, a = 1.0, b = 1.0, incx = 1, incy = 1, n = x.ncols)
  There is also an alias for methods + and - which assume all the default values (accept - assumes a = -1.0)
  The blas fortran routine called will be dictated by the data type in the x Blas matrix.
*  
*/

#include "rb_blas.h"
VALUE rb_blas_xaxpy_mod(int argc, VALUE *argv, VALUE self)
{
	return rb_blas_xaxpy_mod_actual(argc, argv, self, 1);
}

//need to modify the code to use neg * db_c so we get a subtract function.
//Currently, neg is ignored.
VALUE rb_blas_xaxpy_mod_actual( int argc, VALUE *argv, VALUE self, int neg)
{
  Matrix *dx, *dy;
  int incx;
  int incy;
  int n;
  float da_c[2], db_c[2];
  double da_z[2], db_z[2];
  //char error_msg[64];
  VALUE vector_y,  da_value,  db_value,  n_value,  incx_value,  incy_value;
  
  rb_scan_args(argc, argv, "15", &vector_y, &da_value,  &db_value, &incx_value, &incy_value, &n_value);
  
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
  { //sprintf(error_msg, );
    rb_raise(rb_eRuntimeError, "Argument dy is not a Vector");
  }
    
  if(dy->data_type != dx->data_type)
  { //sprintf(error_msg, );
    rb_raise(rb_eRuntimeError, "Vectors are different data_types");
  }
    
  switch(dx->data_type)
  {
  case Single_t: //s
    if(da_value == Qnil)
      da_c[0] = (float) 1.0;
    else
      da_c[0] = (float) NUM2DBL(da_value);
    if(db_value == Qnil)
      cblas_saxpy(n , da_c[0], (float *)dx->data, incx, (float *)dy->data, incy ); 
    else
    {
      db_c[0] = (float) NUM2DBL(db_value);
      catlas_saxpby(n , da_c[0], (float *)dx->data, incx, db_c[0], (float *)dy->data, incy ); 
    }
    break;
  case Double_t: //d
    if(da_value == Qnil)
      da_z[0] = (double) 1.0;
    else
      da_z[0] = NUM2DBL(da_value);
    if(db_value == Qnil)
      cblas_daxpy(n , da_z[0], (double *)dx->data, incx, (double *)dy->data, incy ); 
    else
    {
      db_z[0] = NUM2DBL(db_value);
      catlas_daxpby(n , da_z[0], (double *)dx->data, incx, db_z[0], (double *)dy->data, incy ); 
    }
    break;
  case Complex_t: //c
    if(da_value == Qnil)
    {
      da_c[0] = (float) 1.0;
      da_c[1] = (float) 0.0;
    }
    else
    {
      da_c[0] = (float) NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, da_value) );
      da_c[1] = (float) NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, da_value ) );
    }
    if(db_value == Qnil)
      cblas_caxpy(n , da_c, dx->data, incx, dy->data, incy ); 
    else
    {
      db_c[0] = (float) NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, db_value) );
      db_c[1] = (float) NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, db_value ) );
      catlas_caxpby(n , da_c, dx->data, incx, db_c, dy->data, incy ); 
    }
    break;
  case Double_Complex_t: //z
    if(da_value == Qnil)
    {
      da_z[0] = (double) 1.0;
      da_z[1] = (double) 0.0;
    }
    else
    {
      da_z[0] = NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, da_value) );
      da_z[1] = NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, da_value ) );
    }
    if(db_value == Qnil)
      cblas_zaxpy(n , da_z, dx->data, incx, dy->data, incy ); 
    else
    {
      db_z[0] = NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, db_value) );
      db_z[1] = NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, db_value ) );
      catlas_zaxpby(n , da_z, dx->data, incx, db_z, dy->data, incy ); 
    }
    break;
  default:
    //sprintf(error_msg, );
    rb_raise(rb_eRuntimeError, "Invalid data_type (%d) in Matrix", dx->data_type);
    break; //Never reaches here.
  }

  return vector_y;
}

