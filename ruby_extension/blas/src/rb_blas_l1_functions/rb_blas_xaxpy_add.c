#include "rb_blas.h"

//Quick hack for ruby +
VALUE rb_blas_xaxpy_add(VALUE self, VALUE vector_y)
{ char argc = 1;
  VALUE argv[] = { vector_y };
  
  return rb_blas_xaxpy(argc, argv, self);
}
