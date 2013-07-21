#include "rb_blas.h"
#include "rb_lapack.h"

void Sub_Init_blas_l1(VALUE myClass)
{
  //rotg and rotmg are now classes. Use [S,D]rotg.new and [S,D]rotmg.new
  //These then get passed to rot and rotm as the second parameter.
  rb_define_method(myClass, "rot!", rb_blas_xrot_mod, -1);
  rb_define_method(myClass, "rotm", rb_blas_xrotm_mod, -1);
  rb_define_method(myClass, "axpy", rb_blas_xaxpy, -1);
  rb_define_method(myClass, "axpy!", rb_blas_xaxpy_mod, -1);
  rb_define_method(myClass, "+", rb_blas_xaxpy_add, 1);
  rb_define_method(myClass, "-", rb_blas_xaxpy_sub, 1);
  rb_define_method(myClass, "copy!", rb_blas_xcopy_mod, -1);
  rb_define_method(myClass, "scal", rb_blas_xscal, -1);
  rb_define_method(myClass, "scal!", rb_blas_xscal_mod, -1);
  rb_define_method(myClass, "dot", rb_blas_xdot, -1);  
//  rb_define_method(myClass, "dotu", rb_blas_xdotu, -1);  //included in dot as the default complex case.
  rb_define_method(myClass, "dotc", rb_blas_xdotc, -1);  //only has meaning for complex and double_complex matrix
//  rb_define_method(myClass, "dotsds", rb_blas_xdotsds, -1);  
  rb_define_method(myClass, "amax", rb_blas_ixamax, -1);  
//  rb_define_method(myClass, "max", rb_blas_ixmax, -1);   //Check this exists in our blas
//  rb_define_method(myClass, "amin", rb_blas_ixamin, -1);  //Check this exists
//  rb_define_method(myClass, "min", rb_blas_ixmin, -1);  //Check this exists
  rb_define_method(myClass, "nrm2", rb_blas_xnrm2, -1);  
  rb_define_alias(myClass, "length", "nrm2");  
  rb_define_method(myClass, "asum", rb_blas_xasum, -1);  
  rb_define_method(myClass, "swap!", rb_blas_xswap_mod, -1);  
  rb_define_method(myClass, "set!", rb_atlas_xset_mod, -1);  
 
}
