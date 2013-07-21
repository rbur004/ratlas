
#include "rb_blas.h"

VALUE rb_blas_xscal(int argc, VALUE *argv, VALUE self)
{
  VALUE vector_y_copy;  
  Matrix *dy;
    
  Data_Get_Struct(self, Matrix, dy);
  vector_y_copy = rblas_new_instance(dy->class_id, dy->data, dy->nrows, dy->ncols, dy->data_type, dy->matrix_type, dy->cblas_order);
  return rb_blas_xscal_mod(argc, argv, vector_y_copy);
}

