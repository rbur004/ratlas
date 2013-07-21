#include "rb_blas.h"

static VALUE myClass;

VALUE rblas_integer_rblas_new_instance(void *data, int nrows, int ncols, int matrix_type, int cblas_order)
{
  return rblas_new_instance(myClass, data, nrows, ncols, Integer_t, matrix_type, cblas_order);
}

static VALUE rblas_integer_row_major_matrix_intialize(int argc, VALUE *argv, VALUE obj)
{
  return new_row_major_matrix_helper(GE, Integer_t, argc, argv, obj);
}

static VALUE rblas_integer_column_major_matrix_intialize(int argc, VALUE *argv, VALUE obj)
{
  return new_column_major_matrix_helper(GE, Integer_t, argc, argv, obj);
}

static VALUE integer_blas_initialize(int argc, VALUE *argv, VALUE obj)
{
  Matrix *matrix;
  int data_type;
  
  if(argc == 5) 
  {  
    if((data_type = NUM2INT(argv[2])) != Integer_t)
      rb_raise (rb_eRuntimeError, "Data Type should be Integer_t, but was %s\n", DATA_TYPE_Text[data_type]);
    rb_call_super(argc, argv);
  }
  else  
    new_subclass_init_helper(GE, Integer_t, argc, argv, obj);
  
  Data_Get_Struct(obj, Matrix, matrix);
  matrix->class_id = myClass;
  
  return obj;
}

void Init_IntegerBlas(VALUE myModule, VALUE parent_class)
{
  myClass = rb_define_class_under(myModule, "IntegerBlas", parent_class);
  rb_define_method(myClass, "initialize", integer_blas_initialize, -1);

  // Add a class method to the Ruby class.
  // Create a Matrix, given an array of row arrays, or an Array of elements for the Vector
  rb_define_module_function(myClass, "rows", rblas_integer_row_major_matrix_intialize, -1);
  rb_define_module_function(myClass, "elements", rblas_integer_row_major_matrix_intialize, -1);
  // create a matrix given an array of column arrays
  // Will also create a Vector if given an array of elements, rather than of arrays.
  rb_define_module_function(myClass, "columns", rblas_integer_column_major_matrix_intialize, -1);
  
  rb_define_module_function(myClass, "[]", rblas_integer_row_major_matrix_intialize, -1);
  rb_define_alias(myClass, "[]",  "get");
  rb_define_alias(myClass, "[]=",  "set");

}
