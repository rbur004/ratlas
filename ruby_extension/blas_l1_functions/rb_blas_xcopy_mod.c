
/*
      SUBROUTINE XCOPY(N, X,INCX, Y,INCY)
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
*     copies a vector, x, to a vector, y.
*     uses unrolled loops for increments equal to one.

*  Ruby interface 
*  Y is self
*  X is the argument
*  N, INCX, and INCY are optional, and default to nelements in y, 1, 1
*  i.e.
*   y.xcopy(x, n = y.ncols, incx = 1, incy = 1)
*  Modifies y
*/

static VALUE rb_blas_xcopy_mod(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx, *dy;
  int incx;
  int incy;
  int n;
  char error_msg[64];
  VALUE vector_y, n_value,  incx_value,  incy_value;
  
  rb_scan_args(argc, argv, "13", &vector_y,  &incx_value, &incy_value, &n_value);
  
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
  { sprintf(error_msg, "Self is not a Vector");
    rb_raise(rb_eRuntimeError, error_msg);
  }
  
  if(dy == NULL || dy->ncols != 1)
  { sprintf(error_msg, "Argument dy is not a Vector");
    rb_raise(rb_eRuntimeError, error_msg);
  }
  
  if(dy->data_type != dx->data_type)
  { sprintf(error_msg, "Vectors are different data_types");
    rb_raise(rb_eRuntimeError, error_msg);
  }
    
  switch(dx->data_type)
  {
  case Single_t: //s
    cblas_scopy(n , (float *)dy->data, incy, (float *)dx->data, incx ); 
    break;
  case Double_t: //d
    cblas_dcopy(n , (double *)dy->data, incy, (double *)dx->data, incx ); 
    break;
  case Complex_t: //c
    cblas_ccopy(n , dy->data, incy, dx->data, incx ); 
    break;
  case Double_Complex_t: //z
    cblas_zcopy(n , dy->data, incy, dx->data, incx ); 
    break;
  case Integer_t: // Not at valid type for blas routines
  default:
    sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, error_msg);
    break; //Never reaches here.
  }

  return self;
}
