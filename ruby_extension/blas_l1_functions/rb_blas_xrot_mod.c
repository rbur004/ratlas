
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
  sRotg *srotg;
  dRotg *drotg;
  char error_msg[64];
  VALUE vector_y, rotg,  n_value,  incx_value,  incy_value;

  rb_scan_args(argc, argv, "23", &vector_y, &rotg,  &incx_value, &incy_value, &n_value);

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
  
  if(rotg == Qnil)
  { sprintf(error_msg, "[SD]rotg argument is nil?");
    rb_raise(rb_eRuntimeError, error_msg);
  }

  switch(dx->data_type)
  {
  case Single_t: //s
    Data_Get_Struct(rotg, sRotg, srotg);
    cblas_srot(n , (float *)dx->data, incx, (float *)dy->data, incy, srotg->c, srotg->s); 
    break;
  case Double_t: //d
    Data_Get_Struct(rotg, dRotg, drotg);
    cblas_drot(n , (double *)dx->data, incx, (double *)dy->data, incy, drotg->c, drotg->s); 
    break;
  default:
    sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, error_msg);
    break; //Never reaches here.
  }

  return self;
}

