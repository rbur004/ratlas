#include <ruby.h>
#include "rb_blas.h"
#include "rb_lapack.h"
#include "rotg/rb_rotg.h"


static VALUE myModule;

void Init_ratlas() 
{
  VALUE blass_class, lapack_class, rotg_class;
  // Create a Ruby module.
  myModule = rb_define_module("RAtlas");
  
  blass_class = Init_blas(myModule);
  Init_IntegerBlas(myModule, blass_class);
  Init_SingleBlas(myModule, blass_class);
  Init_DoubleBlas(myModule, blass_class);
  Init_ComplexBlas(myModule, blass_class);
  Init_DoubleComplexBlas(myModule, blass_class);
  
  lapack_class = Init_lapack(myModule, blass_class); 
  Init_SingleLapack(myModule, lapack_class);
  Init_DoubleLapack(myModule, lapack_class);
  
  rotg_class= Init_Rotg(myModule);
  Init_srotg( myModule, rotg_class ) ;
  Init_drotg( myModule, rotg_class ) ;
  Init_crotg( myModule, rotg_class ) ;
  Init_zrotg( myModule, rotg_class ) ;
  Init_srotgm( myModule, rotg_class );
  Init_drotgm( myModule, rotg_class ) ;
}
