#include "rb_blas.h"
#include "rb_lapack.h"


void Sub_Init_blas_l2(VALUE myClass)
{
//"gemv", "gbmv", "hemv","hbmv", "symv", "sbmv", "spmv", "trmv", "tbmv", "tpmv"
  rb_define_method(myClass, "mv", rb_blas_mv, -1);

//"trsv", "tbsv", "tpsv"
  rb_define_method(myClass, "sv", rb_blas_sv, -1);

// "ger", "geru", "her", "hpr", "syr", "spr" 
  rb_define_method(myClass, "r", rb_blas_r, -1);

//"gerc"
  rb_define_method(myClass, "rc", rb_blas_rc, -1);
  
// "her2", "hpr2", "syr2", "spr2"
  rb_define_method(myClass, "r2", rb_blas_r2, -1);
}
