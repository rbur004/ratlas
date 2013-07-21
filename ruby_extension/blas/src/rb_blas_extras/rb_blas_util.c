#include "rb_blas.h"

inline void rb_blas_get_member(void *dest, Matrix *m, int row, int col)
{ //If we always store the array as column major,to be fortran compatible, we could avoid the test.
  if(m->cblas_order == CblasRowMajor)
  //!!!!!!!!!!!!  WARNING: NEED TO Override THE FOLLOWING FOR BANDED and PACKED MATRIX_TYPES    
    memcpy((char *)dest, (char *)m->data + (m->data_size * (row * m->ncols + col)), m->data_size );
  else //CblasColMajor
    memcpy((char *)dest, (char *)m->data + (m->data_size * (col * m->nrows + row)), m->data_size );
}

inline void rb_blas_set_member(Matrix *m, void *src , int row, int col)
{ //If we always store the array as column major,to be fortran compatible, we could avoid the test.
  //!!!!!!!!!!!!  WARNING: NEED TO Override THE FOLLOWING FOR BANDED and PACKED MATRIX_TYPES
  if(m->cblas_order == CblasRowMajor)
    memcpy((char *)m->data + (m->data_size * (row * m->ncols + col)), (char *)src, m->data_size );
  else //CblasColMajor
    memcpy((char *)m->data + (m->data_size * (col * m->nrows + row)), (char *)src, m->data_size );
}

VALUE rb_blas_member_to_value(Matrix *matrix, int row, int col)
{
  double d_value, d_value_i[2];
  float f_value, f_value_i[2];
  int i_value;
  VALUE complex_result[2];
  
  switch(matrix->data_type)
  {
    case  Double_t: rb_blas_get_member((void *)&d_value, matrix,  row,  col); return rb_float_new(d_value);
    case  Single_t: rb_blas_get_member((void *)&f_value, matrix,  row,  col); return rb_float_new((double)f_value);
    case  Integer_t: rb_blas_get_member((void *)&i_value, matrix,  row,  col); return INT2FIX(i_value);
    case  Complex_t:
      rb_blas_get_member((void *)&f_value_i, matrix,  row,  col); 
      complex_result[0] = rb_float_new((double)f_value_i[0]);
      complex_result[1] = rb_float_new((double)f_value_i[1]);
      return rb_class_new_instance(2, complex_result, rb_const_get(rb_mMath, rb_intern("Complex")));
    case  Double_Complex_t: 
      rb_blas_get_member((void *)&d_value_i, matrix,  row,  col);
      complex_result[0] = rb_float_new(d_value_i[0]);
      complex_result[1] = rb_float_new(d_value_i[1]);
      return rb_class_new_instance(2, complex_result, rb_const_get(rb_mMath, rb_intern("Complex")));
    default:  rb_raise(rb_eTypeError, "Unknown data type"); //need to make exception more sensible.
  }
}

void rb_blas_value_to_member(Matrix *matrix, int row, int col, VALUE v)
{
  double d_value, d_value_i[2];
  float f_value, f_value_i[2];
  int i_value;
  
  switch(matrix->data_type)
  {
    case  Double_t: 
      d_value = NUM2DBL(v);
      rb_blas_set_member( matrix,  (void *)&d_value, row,  col); 
      break;
    case  Single_t: 
      f_value = (float) NUM2DBL(v);
      rb_blas_set_member( matrix,  (void *)&f_value, row,  col); 
      break;
    case  Integer_t:
      i_value = NUM2INT(v);
      rb_blas_set_member( matrix,  (void *)&i_value, row,  col); 
      break;
    case  Complex_t:
      f_value_i[0] = (float) NUM2DBL(rb_funcall( v,  rb_intern("real"),  0 ) );
      f_value_i[1] = (float) NUM2DBL(rb_funcall( v,  rb_intern("image"),  0 ) );
      rb_blas_set_member( matrix,  (void *)&f_value_i, row,  col); 
      break;
    case  Double_Complex_t: 
      d_value_i[0] = NUM2DBL(rb_funcall( v,  rb_intern("real"),  0 ));
      d_value_i[1] = NUM2DBL(rb_funcall( v,  rb_intern("image"),  0 ));
      rb_blas_set_member( matrix,  (void *)&d_value_i, row,  col); 
      break;
    default:  rb_raise(rb_eTypeError, "Unknown data type"); //need to make exception more sensible.
  }
}

char * rb_blas_member_to_s(char *buffer, Matrix *matrix, int row, int col)
{
  double d_value, d_value_i[2];
  float f_value, f_value_i[2];
  int i_value;
  
  switch(matrix->data_type)
  {
    case  Double_t: 
      rb_blas_get_member((void *)&d_value, matrix,  row,  col); 
      sprintf(buffer, "%E", d_value);
      break;
    case  Single_t: 
        rb_blas_get_member((void *)&f_value, matrix,  row,  col);
        sprintf(buffer, "%E", f_value);
        break;
    case  Integer_t: 
      rb_blas_get_member((void *)&i_value, matrix,  row,  col);
      sprintf(buffer, "%d", i_value);
      break;
    case  Complex_t:
      rb_blas_get_member((void *)&f_value_i, matrix,  row,  col); 
      if(f_value_i[1] < 0)
        sprintf(buffer, "%E - i%E", f_value_i[0],-f_value_i[1]);
      else
        sprintf(buffer, "%E + i%E", f_value_i[0],f_value_i[1]);
      break;
    case  Double_Complex_t: 
      rb_blas_get_member((void *)&d_value_i, matrix,  row,  col);
      if(d_value_i[1] < 0)
        sprintf(buffer, "%E - i%E", d_value_i[0],-d_value_i[1]);
      else
        sprintf(buffer, "%E + i%E", d_value_i[0],d_value_i[1]);
      break;
    default:  rb_raise(rb_eTypeError, "Unknown data type"); //need to make exception more sensible.
  }
}


#ifdef xxxx
/*
When parsing the incoming arguments from ruby calls to the library, 
we have a lot more information available to us than a fortran function would have.
The Blas object knows the Matrix type, dimensions and other properties, so these need not be passed in.

Many of the Blas routines have more than one result, depending on the UPLO, TRANS, DIAG and SIDE flags.
We need these passed in to us, so we know which variant of the function to call. We could have done this
with multiple methods, but thought it simpler to retain the BLAS convention of setting option flags.

The JumpTable structure tells us what flags we should expect from the user. These are passed from ruby as 
class constants defined in rb_blas.c. e.g. TRANS is Blass::Trans

We need a function to retrieve the arguments, and combine these with the Matrix and call the appropriate function.
There is no nice way of calling a function with variable arguments (that is determining at runtime what arguments
should be in the call, as apposed to using varargs to determine what arguments the function was called with.). We 
therefore have used a less than ideal method and created a function call for each of the possible blas library functions'
argument list. This isn't that many.

*/

void rb_blas_lib_caller(VALUE self, int argc, VALUE *argv, JumpTable *jt)
{
  Matrix *a, *x, *y;
  int incx, incy;
  float  alpha_f[2], temp_f[2];  //need space for complex, so we use and array[2]
  double alpha_d[2], temp_d[2];
  int uplo, trans, diag, side; //possible flags passed in
 
  int n;
  //char error_msg[64];
  VALUE vector_x,  n_value,  incx_value;
  VALUE uplo_value, trans_value, diag_value;

  Data_Get_Struct(self, Matrix, a);


  Data_Get_Struct(vector_x, Matrix, x);

  if(incx_value == Qnil)
    incx = 1;
  else
    incx = NUM2INT(incx_value);

  if(n_value == Qnil)
    n = a->ncols;
  else
    n = NUM2INT(n_value);

  switch(jt->flags)
  {
    case F1: //ruby call has args ( alpha, Vector x, incx, beta, Vector y, incy, TRANS=N)
        rb_scan_args(argc, argv, "15", &vector_x, &n_value,  &incx_value, &trans_value);
      
    case F2: //ruby call has args ( alpha, Vector x, incx, beta, Vector y, incy, TRANS=N)
    case F3: //ruby call has args ( alpha, Vector x, incx, beta, Vector y, incy, TRANS=N)
    case F4:
    case F5:
    case F6:
    case F7:
    case F8:
    default:
      //sprintf(error_msg, "x is not a Vector");
      rb_raise(rb_eRuntimeError, "x is not a Vector");
      break; //never reach here
  }
}

#endif xxxx
