#include "ruby.h"
#include "rb_lapack.h"

//Define the ruby names for the C library functions.
void Sub_Init_lapack(VALUE myClass)
{
  rb_define_method(myClass, "xgesv!", rb_lapack_xgesv_mod, 2);
  rb_define_method(myClass, "xgesv", rb_lapack_xgesv, 2);
  
  rb_define_method(myClass, "xgetrf!", rb_lapack_xgetrf_mod, 1);
 
  rb_define_method(myClass, "xgetrs!", rb_lapack_xgetrs_mod, 3);
  rb_define_method(myClass, "xgetrs", rb_lapack_xgetrs, 3);
  
  rb_define_method(myClass, "determinant", rb_lapack_determinant,0);
  //not in the lapack library
//  rb_define_method(myClass, "determinant!", rb_lapack_determinant_mod, 0);
//  rb_define_method(myClass, "determinant", rb_lapack_determinant, 0);
}

