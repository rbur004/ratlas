static VALUE rb_blas_xaxpy_sub(VALUE self, VALUE vector_y)
{
  char argc = 2;
  VALUE argv[] = { self, rb_float_new(-1.0) };
  
  return rb_blas_xaxpy(argc, argv, vector_y);
}
