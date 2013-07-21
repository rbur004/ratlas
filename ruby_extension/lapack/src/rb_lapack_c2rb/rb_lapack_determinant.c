#include "rb_lapack.h"

static VALUE rb_getrf_mod_helper(VALUE self)
{
  Matrix *m;
  VALUE copy, ipiv;
  //Creates a single VALUE function we can wrap rb_protect around.
  //Return an LU matrix for calculating the determinant. Ignore the ipiv result.
  Data_Get_Struct(self, Matrix, m);
  copy = rblas_new_instance(m->class_id, m->data, m->nrows, m->ncols, m->data_type, m->matrix_type, m->cblas_order);
  ipiv = rblas_integer_rblas_new_instance(NULL, m->nrows, 1, GE, m->cblas_order);
  rb_lapack_xgetrf_mod(copy, ipiv);
  return copy;
}

static VALUE rb_lapack_getrf_rescue(VALUE rescue_func_args, VALUE error_info)
{
  VALUE str = rb_funcall(error_info, rb_intern("to_s"), 0, NULL);
  if(strncmp(StringValuePtr(str), "Singular Matrix", strlen("Singular Matrix")) == 0)
    return Qnil;
  else
    rb_exc_raise(error_info);
}

VALUE rb_lapack_determinant(VALUE self)
{
  Matrix *m;
  VALUE copy;
  //char error_msg[64];
  int i;
  float result_c[2];
  double result_z[2];
  float r_tmp[2];
  double d_tmp[2];
  
  Data_Get_Struct(self, Matrix, m);
  if(m->nrows != m->ncols)
  { //sprintf(error_msg, "Matrix not square");
    rb_raise(rb_eRuntimeError, "Matrix not square");
  }
  
  if(m->matrix_type != LU)
  {
    copy = rb_rescue( rb_getrf_mod_helper, self, rb_lapack_getrf_rescue, Qnil);
    if(copy == Qnil)
      return INT2FIX(0); //any error, we return a determinant of 0
    Data_Get_Struct(copy, Matrix, m);
  }
  
  switch(m->data_type)
  {
  case Single_t: 
    result_c[0] = (float) 1.0;
    for(i=0; i < m->nrows; i++) 
    { rb_blas_get_member(r_tmp, m, i, i);
      result_c[0] *= r_tmp[0];
    }  
    return rb_float_new(result_c[0]);
  case Double_t:
    result_z[0] = 1.0;
    for(i=0; i < m->nrows; i++) 
    { rb_blas_get_member(d_tmp, m, i, i);
      result_z[0] *= d_tmp[0];
    }  
    return rb_float_new(result_z[0]);

  case Complex_t:
  case Double_Complex_t:
    //sprintf(error_msg, "Complex yet to be implemented");
    rb_raise(rb_eRuntimeError, "Complex yet to be implemented");
    break;

  }
    
}