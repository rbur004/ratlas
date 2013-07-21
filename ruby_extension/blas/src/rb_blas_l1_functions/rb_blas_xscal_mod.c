/*
Multiplies each element of a vector by a constant (single-precision).

void cblas_sscal (
   const int N,
   const float alpha,
   float *X,
   const int incX
);

Parameters

n				Number of elements to scale.
alpha		The constant to multiply by.
X				Vector x.
incX	Stride within X. 
	For example, if incX is 7, every 7th element is multiplied by alpha.
	
	xscal(alpha, incx, n)
*/

#include "rb_blas.h"

VALUE rb_blas_xscal_mod(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx;
  int incx;
  int n;
  float da_f;
  double da_d;
  float da_c[2];
  double da_z[2];
  //char error_msg[64];
  VALUE da_value,  n_value,  incx_value;
  
  rb_scan_args(argc, argv, "12", &da_value, &incx_value, &n_value);
  
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
    if(da_value == Qnil)
      da_f = (float) 1.0;
    else
      da_f = (float) NUM2DBL(da_value);
    cblas_sscal(n , da_f, (float *)dx->data, incx ); 
    break;
  case Double_t: //d
    if(da_value == Qnil)
      da_d = (double) 1.0;
    else
      da_d = NUM2DBL(da_value);
    cblas_dscal(n , da_d, (double *)dx->data, incx ); 
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
    cblas_cscal(n , da_c, dx->data, incx ); 
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
    cblas_zscal(n , da_z, dx->data, incx ); 
    break;
  default:
    //sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, "Invalid data_type (%d) in Matrix", dx->data_type);
    break; //Never reaches here.
  }

  return self;
}

