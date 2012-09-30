#include <ruby.h>
#define RB_BLAS_C
#include "rb_blas.h"
#include "blas_l1_c2rb.h"

static VALUE myClass;

static void rblas_mark(Matrix *matrix)
{
  //Empty
}

static void rblas_free(Matrix *matrix)
{
  if(matrix)
  {
    if(matrix->unaligned) 
      free(matrix->unaligned);
    free(matrix);
  }
}

static VALUE rblas_alloc(VALUE klass)
{
  Matrix *matrix = calloc(1, sizeof(Matrix));
  return Data_Wrap_Struct(klass, rblas_mark, rblas_free, matrix); //wraps a c struct in VALUE
  //recover with Data_Get_Struct()
}

static VALUE rblas_init_copy(VALUE copy, VALUE orig) 
{
Matrix *orig_matrix;
Matrix *copy_matrix;
void  *ptr;
  if (copy == orig)
    return copy;
  // we can initialize the copy from other CDPlayers
  // or their subclasses only
  if (  TYPE(orig) != T_DATA 
  || RDATA(orig)->dfree != (RUBY_DATA_FUNC)rblas_free) 
  {
    rb_raise(rb_eTypeError, "wrong argument type");
  }
  //   copy all the fields from the original
  // object's CDJukebox structure to the
  // new object
  Data_Get_Struct(orig, Matrix, orig_matrix);
  Data_Get_Struct(copy, Matrix, copy_matrix);
    
  if(copy_matrix->unaligned != NULL)
    free(copy_matrix->unaligned);
  
  if(orig_matrix->ncols && orig_matrix->nrows)
  {
    copy_matrix->unaligned = ptr = calloc((size_t)(orig_matrix->nrows * orig_matrix->ncols + 16), orig_matrix->data_size );
    copy_matrix->data = ( (size_t)ptr%16 ? (void*) ((size_t)ptr+16-((size_t)ptr%16)) : ptr );    
    MEMCPY(copy_matrix->data, orig_matrix->data, orig_matrix->data_size, orig_matrix->nrows*orig_matrix->ncols );
  }
  else
    copy_matrix->unaligned = copy_matrix->data = NULL;
  copy_matrix->nrows = orig_matrix->nrows;
  copy_matrix->ncols = orig_matrix->ncols;
  copy_matrix->data_type = orig_matrix->data_type;
  copy_matrix->data_size = orig_matrix->data_size;
  copy_matrix->cblas_order = orig_matrix->cblas_order;
  
  return copy;
}

static VALUE rblas_initialize(int argc, VALUE *argv, VALUE obj)
{
  Matrix *matrix;
  VALUE nrows, ncols, matrix_type, data_type, cblas_order;
  void *ptr;
  
  
  rb_scan_args(argc, argv, "14", &nrows, &ncols, &data_type, &matrix_type, &cblas_order);
  
  Data_Get_Struct(obj, Matrix, matrix);

  if(TYPE(nrows) == T_FIXNUM)
  {
    matrix->class_id = myClass; //Remember how we came into this world.
    matrix->nrows =  NUM2INT(nrows);
    
    if(ncols == Qnil) //store Vectors as a 1 column motrix with nrows.
    { 
      matrix->ncols = 1;
    }
    else
      matrix->ncols =  NUM2INT(ncols);
  
    if(data_type == Qnil) //default to type double
      matrix->data_type = Double_t;
    else
      matrix->data_type = NUM2INT(data_type);
  

    if(matrix_type == Qnil) //default to type double
      matrix->matrix_type = GE;
    else
      matrix->matrix_type = NUM2INT(matrix_type);

    if(cblas_order == Qnil) //default to Col Major array storage (i.e. Fortran based array, not C)
      matrix->cblas_order = CblasColMajor;
    else
      matrix->cblas_order = NUM2INT(cblas_order);

    //!!!!!!!!!!!!  WARNING: NEED TO CHANGE THE FOLLOWING FOR BANDED and PACKED MATRIX_TYPES
    switch(matrix->data_type)
    {
      case  Double_t: matrix->data_size = sizeof(double);  break;
      case  Single_t: matrix->data_size = sizeof(float);  break;
      case  Integer_t: matrix->data_size = sizeof(int);  break;
      case  Complex_t: matrix->data_size = sizeof(float)*2;  break;
      case  Double_Complex_t: matrix->data_size = sizeof(double)*2;  break;
      default:  rb_raise(rb_eTypeError, "Unknown data type"); //need to make exception more sensible.
    }
    if(matrix->nrows && matrix->ncols)
    {
      /*intel recommends this hack to ensure 16bit alignment for doubles. It shouldn't hurt for other architectures*/
      matrix->unaligned = ptr = calloc((size_t)(matrix->nrows * matrix->ncols + 16), matrix->data_size);
      matrix->data = ( (size_t)ptr%16 ? (void*) ((size_t)ptr+16-((size_t)ptr%16)) : ptr ); 
    }
    else
      matrix->unaligned = matrix->data = NULL; //Allow a null Vector/Matrix
  }
  else //Make wild assumption that we are attempting to create a copy of another Blas object.
  {
    matrix->class_id = myClass; //Remember how we came into this world.
    matrix->data = NULL;
    return rblas_init_copy(obj, nrows);
  }
  return obj;
}

//Create a Matrix using the storage type of matrix_type, data type of data_type and
//the argv values are the array of values to store.
VALUE new_subclass_init_helper(int matrix_type, int data_type, int argc, VALUE *argv, VALUE obj)
{
  char new_argc = 5;
  VALUE *new_argv = (VALUE *) malloc(sizeof(VALUE)*new_argc);
  VALUE return_value;

  memcpy(new_argv, argv, sizeof(VALUE)*argc);
  if(argc == 1)
    new_argv[1] = Qnil;
  new_argv[2] = INT2FIX(data_type);
  new_argv[3] = INT2FIX(matrix_type);
  new_argv[4] = INT2FIX(CblasColMajor);
  //Assume default of row major.
  return_value = rb_call_super(new_argc, new_argv);
  free(new_argv);
  return return_value;
}

inline void get_member(void *dest, Matrix *m, int row, int col)
{ //If we always store the array as column major,to be fortran compatible, we could avoid the test.
  if(m->cblas_order == CblasRowMajor)
  //!!!!!!!!!!!!  WARNING: NEED TO Override THE FOLLOWING FOR BANDED and PACKED MATRIX_TYPES    
    memcpy(dest, m->data + (m->data_size * (row * m->ncols + col)), m->data_size );
  else //CblasColMajor
    memcpy(dest, m->data + (m->data_size * (col * m->nrows + row)), m->data_size );
}

inline void set_member(Matrix *m, void *src , int row, int col)
{ //If we always store the array as column major,to be fortran compatible, we could avoid the test.
  //!!!!!!!!!!!!  WARNING: NEED TO Override THE FOLLOWING FOR BANDED and PACKED MATRIX_TYPES
  if(m->cblas_order == CblasRowMajor)
    memcpy(m->data + (m->data_size * (row * m->ncols + col)), src, m->data_size );
  else //CblasColMajor
    memcpy(m->data + (m->data_size * (col * m->nrows + row)), src, m->data_size );
}

static VALUE member_to_value(Matrix *matrix, int row, int col)
{
  double d_value, d_value_i[2];
  float f_value, f_value_i[2];
  int i_value;
  VALUE complex_result[2];
  
  switch(matrix->data_type)
  {
    case  Double_t: get_member((void *)&d_value, matrix,  row,  col); return rb_float_new(d_value);
    case  Single_t: get_member((void *)&f_value, matrix,  row,  col); return rb_float_new((double)f_value);
    case  Integer_t: get_member((void *)&i_value, matrix,  row,  col); return INT2FIX(i_value);
    case  Complex_t:
      get_member((void *)&f_value_i, matrix,  row,  col); 
      complex_result[0] = rb_float_new((double)f_value_i[0]);
      complex_result[1] = rb_float_new((double)f_value_i[1]);
      return rb_class_new_instance(2, complex_result, rb_const_get(rb_mMath, rb_intern("Complex")));
    case  Double_Complex_t: 
      get_member((void *)&d_value_i, matrix,  row,  col);
      complex_result[0] = rb_float_new(d_value_i[0]);
      complex_result[1] = rb_float_new(d_value_i[1]);
      return rb_class_new_instance(2, complex_result, rb_const_get(rb_mMath, rb_intern("Complex")));
    default:  rb_raise(rb_eTypeError, "Unknown data type"); //need to make exception more sensible.
  }
}

static void value_to_member(Matrix *matrix, int row, int col, VALUE v)
{
  double d_value, d_value_i[2];
  float f_value, f_value_i[2];
  int i_value;
  
  switch(matrix->data_type)
  {
    case  Double_t: 
      d_value = NUM2DBL(v);
      set_member( matrix,  (void *)&d_value, row,  col); 
      break;
    case  Single_t: 
      f_value = NUM2DBL(v);
      set_member( matrix,  (void *)&f_value, row,  col); 
      break;
    case  Integer_t:
      i_value = NUM2INT(v);
      set_member( matrix,  (void *)&i_value, row,  col); 
      break;
    case  Complex_t:
      f_value_i[0] = NUM2DBL(rb_funcall( v,  rb_intern("real"),  0 ) );
      f_value_i[1] = NUM2DBL(rb_funcall( v,  rb_intern("image"),  0 ) );
      set_member( matrix,  (void *)&f_value_i, row,  col); 
      break;
    case  Double_Complex_t: 
      d_value_i[0] = NUM2DBL(rb_funcall( v,  rb_intern("real"),  0 ));
      d_value_i[1] = NUM2DBL(rb_funcall( v,  rb_intern("image"),  0 ));
      set_member( matrix,  (void *)&d_value_i, row,  col); 
      break;
    default:  rb_raise(rb_eTypeError, "Unknown data type"); //need to make exception more sensible.
  }
}

static char * member_to_s(char *buffer, Matrix *matrix, int row, int col)
{
  double d_value, d_value_i[2];
  float f_value, f_value_i[2];
  int i_value;
  
  switch(matrix->data_type)
  {
    case  Double_t: 
      get_member((void *)&d_value, matrix,  row,  col); 
      sprintf(buffer, "%E", d_value);
      break;
    case  Single_t: 
        get_member((void *)&f_value, matrix,  row,  col);
        sprintf(buffer, "%E", f_value);
        break;
    case  Integer_t: 
      get_member((void *)&i_value, matrix,  row,  col);
      sprintf(buffer, "%d", i_value);
      break;
    case  Complex_t:
      get_member((void *)&f_value_i, matrix,  row,  col); 
      if(f_value_i[1] < 0)
        sprintf(buffer, "%E - i%E", f_value_i[0],-f_value_i[1]);
      else
        sprintf(buffer, "%E + i%E", f_value_i[0],f_value_i[1]);
      break;
    case  Double_Complex_t: 
      get_member((void *)&d_value_i, matrix,  row,  col);
      if(d_value_i[1] < 0)
        sprintf(buffer, "%E - i%E", d_value_i[0],-d_value_i[1]);
      else
        sprintf(buffer, "%E + i%E", d_value_i[0],d_value_i[1]);
      break;
    default:  rb_raise(rb_eTypeError, "Unknown data type"); //need to make exception more sensible.
  }
}


//Creates a Blas class,
//if c_array is not NULL, then it sets the data field of the Blas class to a copy of c_array.
VALUE rblas_new_instance( VALUE self, void *c_array, int nrows, int ncols,  int data_type, int matrix_type, int cblas_order)
{
  char argc = 5;
  VALUE argv[] = { INT2FIX(nrows), INT2FIX(ncols),  INT2FIX(data_type), INT2FIX(matrix_type), INT2FIX(cblas_order)  };
  VALUE obj;
  
  obj = rb_class_new_instance( argc, argv, self); //new instance of this class.
  
  if( c_array != NULL )
  {
    Matrix *m;
    Data_Get_Struct(obj, Matrix, m);
    memcpy(m->data, c_array, m->data_size*nrows*ncols );
  }
  
  return obj;
}


static VALUE rblas_row_major_matrix_intialize(int argc, VALUE *argv, VALUE obj)
{ //Note that this still stores the array in the row or column major form defined in initialize()
  VALUE new_matrix = Qnil;
  Matrix *m;
  int i,r,c;
  int ncols = 0;
  VALUE* row_dataPtr ; //we know the bounds, so work with the row data values.
  VALUE row_v;
  double value;
  int data_type;
  int matrix_type;
  int cblas_order;
  
  if(argc < 2)
    rb_raise (rb_eRuntimeError, "too few arguments (got %d. expected 2 or more)", argc);
    
  matrix_type = NUM2INT(argv[0]); //first argument is the matrix_type.
  data_type = NUM2INT(argv[1]); //secand argument is the data_type.
  cblas_order = CblasColMajor; //We want the internal storage to match the fortran default. Ruby need not know:)

  if(argc == 2 || (argc == 3 && argv[2] == Qnil) )
  {
    new_matrix = rblas_new_instance( obj, NULL, 0, 0, data_type, matrix_type, cblas_order);   
  }
  else
  {
    for(i = 2; i < argc; i++)
    {   
      if (TYPE(argv[i]) == T_ARRAY) //then this is Matrix.
      { 
        row_dataPtr = RARRAY(argv[i])->ptr; //we know the bounds, so work with the col data values.
        if(ncols == 0)
        {
          ncols = RARRAY(argv[i])->len; 
          new_matrix = rblas_new_instance( obj, NULL, argc - 2, ncols, data_type, matrix_type, cblas_order);
          Data_Get_Struct(new_matrix, Matrix, m);
        }
        else if(ncols != RARRAY(argv[i])->len)
          rb_raise(rb_eRuntimeError, "Number of Columns in Matrix varied"); //need to make exception more sensible.
  
        for(c = 0; c < ncols; c++)
        {
          value_to_member(m, i-2, c, row_dataPtr[c]);
        }
      }
      else
      {   //Then we have a Vector
        if(ncols == 0)
        {
         ncols = 1;
         new_matrix = rblas_new_instance( obj, NULL, argc - 2 , ncols, data_type, matrix_type, CblasColMajor); //store Vector as a Fortran Column
         Data_Get_Struct(new_matrix, Matrix, m);
        }
        else if(ncols != 1)
          rb_raise(rb_eRuntimeError, "Number of columns in Vector varied"); //need to make exception more sensible.

        value_to_member(m, i-2, 0, argv[i]);
      }
    }
  }
  
  return new_matrix;
}

static VALUE rblas_column_major_matrix_intialize(int argc, VALUE *argv, VALUE obj)
{ //Note that this still stores the array in the row or column major form defined in initialize()
  VALUE new_matrix = Qnil;
  Matrix *m;
  int i,r,c;
  int nrows = 0;
  VALUE* col_dataPtr ; //we know the bounds, so work with the row data values.
  VALUE row_v;
  double value;
  int data_type;
  int matrix_type;
  int cblas_order;
  
  if(argc < 2)
    rb_raise (rb_eRuntimeError, "too few arguments (expected at least 2)");
    
  matrix_type = NUM2INT(argv[0]);   //first argument is the matrix_type.
  data_type = NUM2INT(argv[1]); //second argument is the data_type.
  cblas_order = CblasColMajor; //We want the internal storage to match the fortran default. Ruby need not know:)
  
  if(argc == 2|| (argc == 3 && argv[2] == Qnil))
  {
    new_matrix = rblas_new_instance( obj, NULL, 0, 0, data_type, matrix_type, cblas_order);   
  }
  else
  {  for(i = 2; i < argc; i++)
    {   
      if (TYPE(argv[i]) == T_ARRAY) //then this is Matrix.
      { 
        col_dataPtr = RARRAY(argv[i])->ptr; //we know the bounds, so work with the col data values.
        if(nrows == 0)
        {
          nrows = RARRAY(argv[i])->len; 
          new_matrix = rblas_new_instance( obj, NULL, nrows, argc - 2, data_type, matrix_type, CblasColMajor);
          Data_Get_Struct(new_matrix, Matrix, m);
        }
        else if(nrows != RARRAY(argv[i])->len)
          rb_raise(rb_eRuntimeError, "Number of rows in Matrix varied"); //need to make exception more sensible.
  
        for(r = 0; r < nrows; r++)
        {
          value_to_member(m, r, i-2,  col_dataPtr[r]);
        }
      }
      else
      {   //Then we have a Vector, which we store as a row
        if(nrows == 0)
        {
         //ncols = 1;
         new_matrix = rblas_new_instance( obj, NULL, argc - 2, 1, data_type, matrix_type, CblasColMajor); //store Vector as a row, not a column
         Data_Get_Struct(new_matrix, Matrix, m);
        }
        else if(nrows != 1)
          rb_raise(rb_eRuntimeError, "Number of rows in Vector varied"); //need to make exception more sensible.

        value_to_member(m,  i-2, 0,  argv[i]);
      }
    }
  }  
  return new_matrix;
}


VALUE new_row_major_matrix_helper(int matrix_type, int data_type, int argc, VALUE *argv, VALUE obj)
{
  int new_argc = argc + 2;
  VALUE *new_argv = (VALUE *) malloc(sizeof(VALUE)*new_argc);
  
  memcpy(&new_argv[2], argv, sizeof(VALUE)*argc);
  new_argv[0] = INT2FIX(matrix_type);
  new_argv[1] = INT2FIX(data_type);
  //Assume default of row major.
  return rblas_row_major_matrix_intialize(new_argc, new_argv, obj);
}

VALUE new_column_major_matrix_helper(int matrix_type, int data_type, int argc, VALUE *argv, VALUE obj)
{
  int new_argc = argc + 2;
  VALUE *new_argv = (VALUE *) malloc(sizeof(VALUE)*new_argc);

  memcpy(&new_argv[2], argv, sizeof(VALUE)*argc);
  new_argv[0] = INT2FIX(matrix_type);
  new_argv[1] = INT2FIX(data_type);
  //Assume default of row major.
  return rblas_column_major_matrix_intialize(new_argc, new_argv, obj);
}


static VALUE rblas_rows(VALUE obj) 
{
  Matrix *m;
 
  Data_Get_Struct(obj, Matrix, m);
  return INT2FIX(m->nrows);
}

static VALUE rblas_columns(VALUE obj) 
{
  Matrix *m;
 
  Data_Get_Struct(obj, Matrix, m);
  return INT2FIX(m->ncols);
}

static VALUE rblas_data_size(VALUE obj)
{
  Matrix *m;
  
  Data_Get_Struct(obj, Matrix, m);
  return INT2FIX(m->data_size);  
}

static VALUE rblas_data_type(VALUE obj)
{
  Matrix *m;
  
  Data_Get_Struct(obj, Matrix, m);
  return INT2FIX(m->data_type);  
}


static VALUE rblas_get_by_index(int argc, VALUE *argv, VALUE obj)
{
  VALUE row,  col;
  Matrix *m;
  int the_row;
  int the_col;

  rb_scan_args(argc, argv, "11",   &row, &col);
  Data_Get_Struct(obj, Matrix, m);
  
  the_row = NUM2INT(row);
  if(col == Qnil)
    the_col = 0;
  else
    the_col = NUM2INT(col);
  
  if( the_row < 0 || the_row >= m->nrows)
    rb_raise(rb_eRuntimeError, "Row value out of bounds"); //need to make exception more sensible.
  if( the_col < 0 || the_col >= m->ncols)
    rb_raise(rb_eRuntimeError, "Col value out of bounds"); //need to make exception more sensible.

  return member_to_value(m, the_row, the_col);
}

static VALUE rblas_set_by_index(int argc, VALUE *argv, VALUE obj) 
{
  VALUE row,  col,  value;
  Matrix *m;
  int the_row;
  int the_col;

  rb_scan_args(argc, argv, "21",   &row, &col, &value);
  Data_Get_Struct(obj, Matrix, m);
  
  the_row = NUM2INT(row);
  if(value == Qnil)
  {
    value = col;
    the_col = 0;
  }
  else
    the_col = NUM2INT(col);
  
  if( the_row < 0 || the_row >= m->nrows)
    rb_raise(rb_eRuntimeError, "Row value out of bounds"); //need to make exception more sensible.
  if( the_col < 0 || the_col >= m->ncols)
    rb_raise(rb_eRuntimeError, "Col value out of bounds"); //need to make exception more sensible.

  value_to_member(m, the_row, the_col, value);
  return value;
}

static VALUE rblas_each_by_row(VALUE obj)
{
  Matrix *m;
  int r,c;
 
  Data_Get_Struct(obj, Matrix, m);
  for(r = 0; r < m->nrows; r++)
    for(c = 0; c < m->ncols; c++)
    { //yield the double at [r,c], and give the row number, and column number.
      rb_yield(  member_to_value(m, r, c) );
    } 
  return Qnil; 
}

//Yeilding each row as an RBlas one row Matrix, and the row index
//Unless, there is only one row, and then it returns the values.

static VALUE rblas_each_row(VALUE obj) 
{
  Matrix *m;
  int c,r;
  void *row_data;
 
  Data_Get_Struct(obj, Matrix, m);
  

  if(m->ncols >= 1)
  {
    if((row_data = calloc(m->ncols, m->data_size)) == NULL)
      rb_raise(rb_eRuntimeError, "Can't allocate memory"); //need to make exception more sensible.
        
    for(r = 0; r < m->nrows; r++)
    {
      for(c = 0; c < m->ncols; c++)
      {
        get_member(row_data + (c * m->data_size), m, r, c);
      }
      rb_yield( rb_ary_new3( 2, rblas_new_instance( m->class_id, row_data,  m->ncols, 1, m->data_type, m->matrix_type, m->cblas_order), INT2FIX(r) ) );
    } 
  
    free(row_data);
  }
  return Qnil;  
}



static VALUE rblas_each_by_row_with_index(VALUE obj)
{
  Matrix *m;
  int r,c;
 
  Data_Get_Struct(obj, Matrix, m);
  for(r = 0; r < m->nrows; r++)
    for(c = 0; c < m->ncols; c++)
    { //yield the double at [r,c], and give the row number, and column number.
      rb_yield( rb_ary_new3( 3, member_to_value(m, r, c), INT2FIX(r), INT2FIX(c) ) );
    } 
  return Qnil; 
}


static VALUE rblas_each_by_col(VALUE obj)
{
  Matrix *m;
  int r,c;
 
  Data_Get_Struct(obj, Matrix, m);
  for(c = 0; c < m->ncols; c++)
    for(r = 0; r < m->nrows; r++)
    { //yield the double at [r,c]
      rb_yield( member_to_value(m, r, c)  );
    } 
  return Qnil; 
}

//Yeilding each column as an RBlas one row Matrix, and the column index
//Unless, there is only one column, and then it returns the values.
static VALUE rblas_each_col(VALUE obj)
{
  Matrix *m;
  int c,r;
  void *column_data;
 
  Data_Get_Struct(obj, Matrix, m);
  if(m->ncols >= 1)
  {
    if((column_data = calloc(m->nrows, m->data_size)) == NULL)
      rb_raise(rb_eRuntimeError, "Can't allocate memory"); //need to make exception more sensible.
        
    for(c = 0; c < m->ncols; c++)
    {
      for(r = 0; r < m->nrows; r++)
      {
        get_member(column_data + (r * m->data_size), m, r, c);
      }
      rb_yield( rb_assoc_new(rblas_new_instance( m->class_id, column_data, m->nrows, 1,  m->data_type, m->matrix_type, m->cblas_order), INT2FIX(c) ) );
    } 
  
    free(column_data);
  }
  return Qnil;  
}

static VALUE rblas_each_by_col_with_index(VALUE obj)
{
  Matrix *m;
  int r,c;
 
  Data_Get_Struct(obj, Matrix, m);
  for(c = 0; c < m->ncols; c++)
    for(r = 0; r < m->nrows; r++)
    { //yield the double at [r,c], and give the row number, and column number.
      rb_yield( rb_ary_new3( 3, member_to_value(m, r, c), INT2FIX(r), INT2FIX(c) ) );
    } 
  return Qnil; 
}

static VALUE rblas_abs(VALUE obj)
{
  //Definitely a hack. We want to be able to test the results with the same syntax as scalars
  //i.e. if (error = (result - expected)).abs > bound ...
  Matrix *m, *abs_m;
  int r,c;
  VALUE abs_matrix;

  Data_Get_Struct(obj, Matrix, m);
  //argv[0] = rblas_new_instance(dy->class_id, dy->data, dy->nrows, dy->ncols, dy->data_type);
  abs_matrix = rblas_new_instance( m->class_id, NULL, m->nrows, m->ncols, m->data_type, m->matrix_type, m->cblas_order);
  Data_Get_Struct(abs_matrix, Matrix, abs_m);

  for(r = 0; r < m->nrows; r++)
  {
    for(c = 0; c < m->ncols; c++)
    { //use ruby to do this, as it will happily handle the complex data types.
      value_to_member(abs_m, r, c, rb_funcall(member_to_value(m, r, c),  rb_intern("abs"), 0) );
    } 
  }

  return abs_matrix;
}

static VALUE rblas_greater(VALUE obj, VALUE bound)
{
  //Definitely a hack. We want to be able to test the results with the same syntax as scalars
  //i.e. if (error = (result - expected)).abs > bound ...
  Matrix *m;
  int r,c;

  Data_Get_Struct(obj, Matrix, m);

  for(r = 0; r < m->nrows; r++)
  {
    for(c = 0; c < m->ncols; c++)
    { //use ruby to do this, as it will happily handle the complex data types.
      if (TYPE(rb_funcall(member_to_value(m, r, c),rb_intern(">"), 1, bound)) == TYPE(Qfalse)) return Qfalse;
    } 
  }
  return Qtrue;
}

static VALUE rblas_to_a(VALUE obj)
{
  Matrix *m;
  int r,c;
  VALUE matrix_array, row_array;
  
  Data_Get_Struct(obj, Matrix, m);
  
  matrix_array = rb_ary_new2(m->nrows); //we will add an array[ncols] per row.
  
  if(m->ncols == 1) //a vector, so just store it in the matrix_array
  {
    for(r = 0; r < m->nrows; r++)
      rb_ary_store(matrix_array, r, member_to_value(m, r, 0) );
  }
  else
  {
    for(r = 0; r < m->nrows; r++)
    {
      row_array = rb_ary_new2(m->ncols); //create an array to hold this row.
      rb_ary_store(matrix_array, r, row_array); //add the row to the matrix array.

      for(c = 0; c < m->ncols; c++)
      { 
        rb_ary_store(row_array, c, member_to_value(m, r, c) ); //add in each column element of this row
      } 
    }
  }
  
  return matrix_array;
}

static VALUE rblas_to_s(VALUE obj)
{
  Matrix *m;
  int r,c;
  char buffer[16];
  VALUE s;
 
  Data_Get_Struct(obj, Matrix, m);
  
  if(m->ncols == 1) //a vector, so just store it in the matrix_array
  {
    s = rb_str_new2("[ ");
    for(r = 0; r < m->nrows; r++)
    {
      member_to_s(buffer, m, r, 0);
      rb_str_cat(s, buffer, strlen(buffer) );
      if(r != m->nrows-1)
        rb_str_cat(s, ", ", 2 );
    }
  }
  else
  {
    s = rb_str_new2("[\n");
    for(r = 0; r < m->nrows; r++)
    {
      rb_str_cat(s, "  [ ", 4 );
      for(c = 0; c < m->ncols; c++)
      { 
        member_to_s(buffer, m, r, c);
        rb_str_cat(s, buffer, strlen(buffer) );
        if(c != m->ncols-1)
          rb_str_cat(s, ", ", 2 );
      } 
      rb_str_cat(s, " ]\n", 3 ); 
    }
  }
  rb_str_cat(s, "]", 1 );
  return s;
}


// This is called when the Ruby interpreter loads this C extension.
// The part after "Init_" is the name of the C extension specified
// in extconf.rb, not the name of the C source file.

VALUE Init_blas(VALUE myModule) 
{
  (void) rb_require( "complex" );
  // Create a Ruby class in this module.
  // rb_cObject is defined in ruby.h
  
  myClass = rb_define_class_under(myModule, "Blas", rb_cObject);
  rb_define_alloc_func(myClass, rblas_alloc); 
  // The first and second arguments are the dimensions of the Matrix, nrows , ncolumns.
  // the third arg is optional and is the data type i.e. Double_t(default), Single_t, Integer_t, Double_Complex_t, Complex_t
  // The last arg is optional and is the matrix_type. i.e. GB, GE (default), HB, HE, HP, SB, SP, SY, TB, TP, TR, US 
  rb_define_method(myClass, "initialize", rblas_initialize, -1);
  rb_define_method(myClass, "initialize_copy", rblas_init_copy, 1);

  // Add a class method to the Ruby class.
  // Create a Matrix, given an array of row arrays, or an Array of elements for the Vector
  // The first arg is the data type to us in C. Double_t, Single_t, Integer_t, Double_Complex_t, Complex_t
  // The first arg is the matrix_type. i.e. GB, GE, HB, HE, HP, SB, SP, SY, TB, TP, TR, US 
  // the second is the data type i.e. Double_t, Single_t, Integer_t, Double_Complex_t, Complex_t
  rb_define_module_function(myClass, "rows", rblas_row_major_matrix_intialize, -1); 
  rb_define_module_function(myClass, "elements", rblas_row_major_matrix_intialize, -1);
  // create a matrix given an array of column arrays
  // Will also create a Vector if given an array of elements, rather than of arrays.
  // The first arg is the matrix_type. i.e. GB, GE, HB, HE, HP, SB, SP, SY, TB, TP, TR, US 
  // the second is the data type i.e. Double_t, Single_t, Integer_t, Double_Complex_t, Complex_t
  rb_define_module_function(myClass, "columns", rblas_column_major_matrix_intialize, -1);


  // Add an instance method to the Ruby class.

  rb_define_method(myClass, "rows", rblas_rows, 0);
  rb_define_alias(myClass, "nrows",  "rows");
  rb_define_method(myClass, "columns", rblas_columns, 0);
  rb_define_alias(myClass, "ncols",  "columns");

  rb_define_method(myClass, "data_size", rblas_data_size, 0);
  rb_define_method(myClass, "data_type", rblas_data_type, 0);
  
  rb_define_method(myClass, "[]", rblas_get_by_index, -1);
  rb_define_method(myClass, "get", rblas_get_by_index, -1);
  rb_define_method(myClass, "[]=", rblas_set_by_index, -1);
  rb_define_method(myClass, "set", rblas_set_by_index, -1);

  rb_define_method(myClass, "each_row", rblas_each_row, 0);
  rb_define_method(myClass, "each_col", rblas_each_col, 0);
  rb_define_method(myClass, "each_by_row", rblas_each_by_row, 0);
  rb_define_alias( myClass, "each",  "each_by_row");
  rb_define_method(myClass, "each_by_row_with_index", rblas_each_by_row_with_index, 0);
  rb_define_method(myClass, "each_by_col", rblas_each_by_col, 0);
  rb_define_method(myClass, "each_by_col_with_index", rblas_each_by_col_with_index, 0);
  
  //Couple of hacks to test the matrix results are in bounds.
  //The intent is to get the same behaviour as testing a scalar is in bounds.
  rb_define_method(myClass, "abs", rblas_abs, 0); //returns abs values of the array elements for testing results are in bounds.
  rb_define_method(myClass, ">", rblas_greater, 1);  //validates that Matricies elements are all > bound testing results are in bounds.
  
  rb_define_method(myClass, "to_a", rblas_to_a, 0); //returns a ruby array from the matrix.
  rb_define_method(myClass, "to_s", rblas_to_s, 0); //returns a ruby string.
  
  //Class constants for use in creating the Matrix
  rb_define_const( myClass, "Double_t", INT2FIX(Double_t) );
  rb_define_const( myClass, "Single_t", INT2FIX(Single_t) );
  rb_define_const( myClass, "Integer_t", INT2FIX(Integer_t) );
  rb_define_const( myClass, "Complex_t", INT2FIX(Complex_t) );
  rb_define_const( myClass, "Double_Complex_t", INT2FIX(Double_Complex_t) );
  rb_define_const( myClass, "D", INT2FIX(Double_t) );
  rb_define_const( myClass, "S", INT2FIX(Single_t) );
  rb_define_const( myClass, "I", INT2FIX(Integer_t) );
  rb_define_const( myClass, "C", INT2FIX(Complex_t) );
  rb_define_const( myClass, "Z", INT2FIX(Double_Complex_t) );
  
  //Class constants from the cblas.h definitions.
  rb_define_const( myClass, "GB", INT2FIX(GB));
  rb_define_const( myClass, "GE", INT2FIX(GE));
  rb_define_const( myClass, "General_Banded", INT2FIX(GB));
  rb_define_const( myClass, "General", INT2FIX(GE));
  
  rb_define_const( myClass, "HB", INT2FIX(HB));
  rb_define_const( myClass, "HE", INT2FIX(HE));
  rb_define_const( myClass, "HP", INT2FIX(HP));
  rb_define_const( myClass, "Hermitian_Banded", INT2FIX(HB));
  rb_define_const( myClass, "Hermitian", INT2FIX(HE));
  rb_define_const( myClass, "Hermitian_Packed", INT2FIX(HP));
  
  rb_define_const( myClass, "SB", INT2FIX(SB));
  rb_define_const( myClass, "SY", INT2FIX(SY));
  rb_define_const( myClass, "SP", INT2FIX(SP));
  rb_define_const( myClass, "Symmetric_banded", INT2FIX(SB));
  rb_define_const( myClass, "Symmetric", INT2FIX(SY));
  rb_define_const( myClass, "Symmetric_Packed", INT2FIX(SP));
  
  rb_define_const( myClass, "TB", INT2FIX(TB));
  rb_define_const( myClass, "TR", INT2FIX(TR));
  rb_define_const( myClass, "TP", INT2FIX(TP));
  rb_define_const( myClass, "Triangular_Banded", INT2FIX(TB));
  rb_define_const( myClass, "Triangular", INT2FIX(TR));
  rb_define_const( myClass, "Triangular_Packed", INT2FIX(TP));

  rb_define_const( myClass, "RowMajor", INT2FIX(CblasRowMajor));
  rb_define_const( myClass, "ColMajor", INT2FIX(CblasColMajor));
  
  rb_define_const( myClass, "NoTrans", INT2FIX(CblasNoTrans));
  rb_define_const( myClass, "Trans", INT2FIX(CblasTrans));
  rb_define_const( myClass, "TransA", INT2FIX(CblasTrans));
  rb_define_const( myClass, "ConjTrans", INT2FIX(CblasConjTrans));
  rb_define_const( myClass, "TransB", INT2FIX(CblasConjTrans));
  rb_define_const( myClass, "AtlasConj", INT2FIX(AtlasConj));
  rb_define_const( myClass, "Upper", INT2FIX(CblasUpper));
  rb_define_const( myClass, "Lower", INT2FIX(CblasLower));
  rb_define_const( myClass, "NonUnit", INT2FIX(CblasNonUnit));
  rb_define_const( myClass, "Unit", INT2FIX(CblasUnit));
  rb_define_const( myClass, "Left", INT2FIX(CblasLeft));
  rb_define_const( myClass, "Right", INT2FIX(CblasRight));
  
  Sub_Init_blas_l1(myClass);
  Sub_Init_blas_l2(myClass);
  Sub_Init_blas_l3(myClass);
  
  return myClass;
}