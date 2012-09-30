#include "rb_blas.h"

#ifdef xxxx
#include "rb_blas_util.h"
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

void blas_lib_caller(VALUE self, int argc, VALUE *argv, JumpTable *jt)
{
  Matrix *a, *x, *y;
  int incx, incy;
  float  alpha_f[2], temp_f[2];  //need space for complex, so we use and array[2]
  double alpha_d[2], temp_d[2];
  int uplo, trans, diag, side; //possible flags passed in
 
  int n;
  char error_msg[64];
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
      sprintf(error_msg, "x is not a Vector");
      rb_raise(rb_eRuntimeError, error_msg);
      break; //never reach here
  }
}

#endif xxxx
