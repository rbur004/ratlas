// This is called when the Ruby interpreter loads this C extension.
// The part after "Init_" is the name of the C extension specified
// in extconf.rb, not the name of the C source file.
#include "rb_rotg.h"

static VALUE myClass;


static void rblas_Rotg_mark(Rotg *rotg)
{
  //Empty
}

static void rblas_Rotg_free(Rotg *rotg)
{
  if(rotg)
  {
    free(rotg);
  }
}

static VALUE rblas_Rotg_alloc(VALUE klass)
{
  Rotg *rotg = calloc(1, sizeof(Rotg));
  return Data_Wrap_Struct(klass, rblas_Rotg_mark, rblas_Rotg_free, rotg); //wraps a c struct in VALUE
  //recover with Data_Get_Struct()
}

static VALUE rblas_Rotg_initialize(int argc, VALUE *argv, VALUE obj)
{
  Rotg *rotg;  
    
  Data_Get_Struct(obj, Rotg, rotg);
  rotg->class_id = myClass; //Remember how we came into this world.
  
  return obj;
}

VALUE Init_Rotg(VALUE myModule) 
{
  // Create a Ruby class in this module.
  // rb_cObject is defined in ruby.h
  
  myClass = rb_define_class_under(myModule, "Rotg", rb_cObject);
  rb_define_alloc_func(myClass, rblas_Rotg_alloc); 
  rb_define_method(myClass, "initialize", rblas_Rotg_initialize, -1);

  Init_srotg( myModule, myClass ) ;
  Init_drotg( myModule, myClass ) ;
  Init_crotg( myModule, myClass ) ;
  Init_zrotg( myModule, myClass ) ;
  Init_srotmg( myModule, myClass );
  Init_drotmg( myModule, myClass ) ;

  return myClass;
  
}