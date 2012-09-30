static VALUE rb_blas_xaxpy(int argc, VALUE *argv, VALUE self)
{Matrix *dy;
    
  Data_Get_Struct(argv[0], Matrix, dy);
  argv[0] = rblas_new_instance(dy->class_id, dy->data, dy->nrows, dy->ncols, dy->data_type, dy->matrix_type, dy->cblas_order);
  return rb_blas_xaxpy_mod(argc, argv, self);
}
