
/*Need to decide what this is. A dot or cross product?*/

static VALUE rb_blas_times(VALUE self, VALUE multiplier)
{
  char argc = 1;
  VALUE argv[] = { multiplier };
  Matrix *dx, *dy;
  
  switch(TYPE(multiplier))
  {
  case T_FIXNUM:
  case T_BIGNUM:
  case T_FLOAT:
    Data_Get_Struct(self, Matrix, dx);
    if(dx->ncols == 1)
      return rb_blas_xscal(argc, argv, self); //Vector * Scalar
    else
      rb_raise(rb_eRuntimeError, "Matrix * Scalar not implemented");
  default:
    Data_Get_Struct(self, Matrix, dx);
    Data_Get_Struct(multiplier, Matrix, dy);
    if(dx->ncols == 1)
    {
      if(dy->ncols == 1)
        return rb_blas_xdot(argc, argv, self ); //Vector * Vector (dot product)
      else
        rb_raise(rb_eRuntimeError, "Vector * Matrix not implemented");
    }
    else
    {      
      rb_raise(rb_eRuntimeError, "Matrix * Vector and Matrix * Matrix not implemented");
#ifdef gotaroundtoit      
      if(dy->ncols == 1)
        return rb_blas_xtrmv(argc,argv, self); //Matrix * Vector
      else
        return rb_blas_xtrmm(argc,argv, self); //Matrix * Matrix
#endif
    }  
  }
  return Qnil; //Never reaches here.
}
