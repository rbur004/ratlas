//SingleLapack class. Inherits from Lapack
#include "rb_lapack.h"

static VALUE myClass;


static VALUE single_row_major_matrix_initialize(int argc, VALUE *argv, VALUE obj)
{
  return new_row_major_matrix_helper(GE, Single_t, argc, argv, obj);
}

static VALUE single_column_major_matrix_initialize(int argc, VALUE *argv, VALUE obj)
{
  return new_column_major_matrix_helper(GE, Single_t, argc, argv, obj);
}

static VALUE single_lapack_initialize(int argc, VALUE *argv, VALUE obj)
{
  Matrix *matrix;
  int data_type;
    
  if(argc == 5) 
  {  
    if((data_type = NUM2INT(argv[2])) != Single_t)
      rb_raise (rb_eRuntimeError, "Data Type should be Single_t, but was %s\n", DATA_TYPE_Text[data_type]);
    rb_call_super(argc, argv);
  }
  else
    new_subclass_init_helper(GE, Single_t, argc, argv, obj);
  
  Data_Get_Struct(obj, Matrix, matrix);
  matrix->class_id = myClass;
  
  return obj;
}

void Init_SingleLapack(VALUE myModule, VALUE parent_class)
{
  myClass = rb_define_class_under(myModule, "SingleLapack", parent_class);
  rb_define_method(myClass, "initialize", single_lapack_initialize, -1);

   // Add a class method to the Ruby class.
   // Create a Matrix, given an array of row arrays, or an Array of elements for the Vector
   rb_define_module_function(myClass, "rows", single_row_major_matrix_initialize, -1);
   rb_define_module_function(myClass, "elements", single_row_major_matrix_initialize, -1);
   // create a matrix given an array of column arrays
   // Will also create a Vector if given an array of elements, rather than of arrays.
   rb_define_module_function(myClass, "columns", single_column_major_matrix_initialize, -1);

   rb_define_module_function(myClass, "[]", single_row_major_matrix_initialize, -1);
   rb_define_alias(myClass, "[]",  "get");
   rb_define_alias(myClass, "[]=",  "set");
   
}