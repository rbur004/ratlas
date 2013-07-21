#include <ruby.h>
#include "rb_blas.h"
#include "rb_lapack.h"
#include "rb_rotg.h"


static VALUE myModule;

void Init_ratlas() 
{
  VALUE blass_class, lapack_class, rotg_class;
  // Create a Ruby module.
  myModule = rb_define_module("RAtlas");
  
  blass_class = Init_blas(myModule);
  
  lapack_class = Init_lapack(myModule, blass_class); 
  
  rotg_class= Init_Rotg(myModule);
}
