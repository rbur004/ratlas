
/*SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S)
*     .. Scalar Arguments ..
      REAL C,S
      INTEGER INCX,INCY,N
*     ..
*     .. Vector Arguments ..
      REAL SX(*),SY(*)
*     ..
*
*  Purpose
*  =======
*
*     applies a plane rotation.
*     i.e. SX <- C*SX + S*SY
*          SY <- C*SY - S*SX
*/
static VALUE rb_blas_xrot_mod(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx, *dy;
  int incx;
  int incy;
  int n;
  float c_f[2], s_f[2];
  double c_d[2], s_d[2];
  char error_msg[64];
  VALUE vector_y,  c_value, s_value,  n_value,  incx_value,  incy_value;

  rb_scan_args(argc, argv, "15", &vector_y, &c_value,  &s_value,  &incx_value, &incy_value, &n_value);

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
    if(c_value == Qnil)
      c_f[0] = 1.0;
    else
      c_f[0] = NUM2DBL(c_value);
    if(s_value == Qnil)
      s_f[0] = 1.0;
    else
      s_f[0] = NUM2DBL(s_value);
    cblas_srot(n , (float *)dx->data, incx, (float *)dy->data, incy, c_f[0], s_f[0]); 
    break;
  case Double_t: //d
    if(c_value == Qnil)
      c_d[0] = 1.0;
    else
      c_d[0] = NUM2DBL(c_value);
    if(s_value == Qnil)
      s_d[0] = 1.0;
    else
      s_d[0] = NUM2DBL(s_value);
    cblas_drot(n , (double *)dx->data, incx, (double *)dy->data, incy, c_d[0], s_d[0]); 
    break;
  default:
    sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, error_msg);
    break; //Never reaches here.
  }

  return self;
}

