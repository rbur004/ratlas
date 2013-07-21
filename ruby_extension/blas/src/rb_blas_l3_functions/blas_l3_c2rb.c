#include "rb_blas.h"
#include "rb_lapack.h"

void Sub_Init_blas_l3(VALUE myClass)
{
  rb_define_method(myClass, "gemm", rb_blas_xgemm, -1);
  rb_define_method(myClass, "symm", rb_blas_xsymm, -1);
  rb_define_method(myClass, "hemm!", rb_blas_xhemm, -1);
  rb_define_method(myClass, "syrk", rb_blas_xsyrk, -1);
  rb_define_method(myClass, "herk", rb_blas_xherk, -1);
  rb_define_method(myClass, "syr2k", rb_blas_xsyr2k, -1);
  rb_define_method(myClass, "her2k", rb_blas_xher2k, -1);
  rb_define_method(myClass, "trmm", rb_blas_xtrmm, -1);
  rb_define_method(myClass, "trsm", rb_blas_xtrsm, -1);  
}
