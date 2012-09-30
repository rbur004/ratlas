#include <ruby.h>
#include "rb_blas.h"
#include "rb_lapack.h"
#include "lapack_c2rb.h"


static VALUE myClass;

static VALUE lapack_initialize(int argc, VALUE *argv, VALUE obj)
{
  Matrix *matrix;
  
  rb_call_super(argc,argv);
  
  Data_Get_Struct(obj, Matrix, matrix);
  matrix->class_id = myClass;
  
  return obj;
}

static VALUE lapack_to_s(VALUE obj)
{
  return rb_str_new2("Not Implemented");
}

VALUE Init_lapack(VALUE myModule, VALUE parent_class)
{
  myClass = rb_define_class_under(myModule, "Lapack", parent_class);
  rb_define_method(myClass, "initialize", lapack_initialize, -1);
  //rb_define_method(myClass, "to_s", lapack_to_s, 0);
  
  Sub_Init_lapack(myClass);
  
  return myClass;
}


