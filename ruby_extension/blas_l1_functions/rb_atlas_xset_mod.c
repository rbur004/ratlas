
static VALUE rb_atlas_xset_mod(int argc, VALUE *argv, VALUE self)
{
  Matrix *dx;
  int incx;
  int n;
  float a_f;
  double a_d;
  float a_c[2];
  double a_z[2];
  char error_msg[64];
  VALUE a_value,  n_value,  incx_value;
  
  rb_scan_args(argc, argv, "12", &a_value, &incx_value, &n_value);
  
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
  { sprintf(error_msg, "Self is not a Vector");
    rb_raise(rb_eRuntimeError, error_msg);
  }
    
  switch(dx->data_type)
  {
  case Single_t: //s
    if(a_value == Qnil)
      a_f = 1.0;
    else
      a_f = NUM2DBL(a_value);
    catlas_sset(n , a_f, (float *)dx->data, incx ); 
    break;
  case Double_t: //d
    if(a_value == Qnil)
      a_d = 1.0;
    else
      a_d = NUM2DBL(a_value);
    catlas_dset(n , a_d, (double *)dx->data, incx); 
    break;
  case Complex_t: //c
    if(a_value == Qnil)
    {
      a_c[0] = 1.0;
      a_c[1] = 0.0;
    }
    else
    {
      a_c[0] = NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, a_value) );
      a_c[1] = NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, a_value ) );
    }
    catlas_cset(n , a_c, dx->data, incx ); 
    break;
  case Double_Complex_t: //z
    if(a_value == Qnil)
    {
      a_z[0] = 1.0;
      a_z[1] = 0.0;
    }
    else
    {
      a_z[0] = NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, a_value) );
      a_z[1] = NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, a_value ) );
    }
    catlas_zset(n , a_z, dx->data, incx); 
    break;
  default:
    sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, error_msg);
    break; //Never reaches here.
  }

  return self;
}

