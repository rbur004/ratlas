/*
FORMAT
  {S,D,C,Z}ROT (x, y, rotg, incx, incy, n) 

Arguments


x                    real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array X of length at least
                      (1+(n-1)*|incx|), containing the elements of the vector
                      x.
                      On exit, if n<=0 or if c is 1.0 and s is 0.0, x is
                      unchanged.  Otherwise, x is overwritten; X contains the
                      rotated vector x.

 y                    real*4 | real*8 | complex*8 | complex*16
                      On entry, a one-dimensional array Y of length at least
                      (1+(n-1)*|incy|).  Y contains the n elements of the
                      vector y.
                      On exit, if n<=0 or if c is 1.0 and s is 0.0, y is
                      unchanged.  Otherwise, y is overwritten; Y contains the
                      rotated vector y.
                      
rotg                  class created by [SD]rotg.new

Optional Arguments

incx                 integer*4 can be nil and will default to 1.
                    On entry, the increment for the array X.
                    If incx >= 0, vector x is      stored forward in the array,
                    so that x(i) is stored in location X(1+(i-1)*incx).
                    If incx < 0, vector x is stored backward in the array,
                    so that x(i) is stored in location X(1+(n-i)*|incx|).
                    On exit, incx is unchanged.

incy                 integer*4 can be nil and will default to 1.
                    On entry, the increment for the array Y.
                    If incy >= 0, vector y is      stored forward in the array,
                    so that y(i) is stored in location Y(1+(i-1)*incy).
                    If incy < 0, vector y is stored backward in the array,
                    so that y(i) is stored in location Y(1+(n-i)*|incy|).
                    On exit, incy is unchanged.

n                    integer*4 can be nil and will default to size of x.
                    On entry, the number of elements in the vectors x and
                    y.
                    On exit, n is unchanged.

Description
  SROT and DROT apply a real Givens plane rotation to each element in the
  pair of real vectors, x and y. CSROT and ZDROT apply a real Givens plane
  rotation to elements in the complex vectors, x and y.  CROT and ZROT apply
  a complex Givens plane rotation to each element in the pair of complex vec-
  tors x and y.

  The cosine and sine of the angle of rotation are c and s, respectively, and
  are provided by the BLAS Level 1 _ROTG subroutines.

  The Givens plane rotation for SROT, DROT, CSROT, and ZDROT follows: x(i) =
  c*x(i) + s*y(i) y(i) = -s*x(i) + c*y(i)

  The elements of the rotated vector x are x(i)  = cx(i) + sy(i).
  The elements of the rotated vector y are y(i)  =  -sx(i) + cy(i).

  The Givens plane rotation for CROT and ZROT follows: x(i) = c*x(i) + s*y(i)
  y(i) = -conjugate(s)*x(i) + c*y(i)

  The elements of the rotated vector x are x(i)  = cx(i) + sy(i).
  The elements of the rotated vector y are y(i)  =  -conjugate(s)x(i) +
  cy(i).

  If n<=0 or if      c = 1.0 and s = 0.0, x and y are unchanged.  If any element
  of x shares a memory location with an element of y, the results are
  unpredictable.

  These subroutines can be used to introduce zeros selectively into a matrix.
*/
static VALUE rb_blas_xrotm_mod(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx, *dy;
  int incx;
  int incy;
  int n;
  sRotmg *srotmg;
  dRotmg *drotmg;
  char error_msg[64];
  VALUE vector_y, rotmg,  n_value,  incx_value,  incy_value;

  rb_scan_args(argc, argv, "23", &vector_y, &rotmg,  &incx_value, &incy_value, &n_value);

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
  
  if(rotmg == Qnil)
  { sprintf(error_msg, "[SD]rotg argument is nil?");
    rb_raise(rb_eRuntimeError, error_msg);
  }

  switch(dx->data_type)
  {
  case Single_t: //s
    Data_Get_Struct(rotmg, sRotmg, srotmg);
    cblas_srotm(n , (float *)dx->data, incx, (float *)dy->data, incy, srotmg->param ); 
    break;
  case Double_t: //d
    Data_Get_Struct(rotmg, dRotmg, drotmg);
    cblas_drotm(n , (double *)dx->data, incx, (double *)dy->data, incy, drotmg->param ); 
    break;
  default:
    sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, error_msg);
    break; //Never reaches here.
  }

  return self;
}

