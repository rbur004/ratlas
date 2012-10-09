/*
FORMAT
  {S,D}rotmg.new(d1, d2, x1, y1, param)

Arguments

  d1                   real*4 | real*8
                      On entry, the first scale factor for the modified
                      Givens transform.
                      On exit, d1 is updated.

  d2                   real*4 | real*8
                      On entry, the second scale factor for the modified
                      Givens transform.
                      On exit, d2 is updated.

  x1                   real*4 | real*8
                      On entry, the first element x(1) of the input vector.
                      On exit, x1 is overwritten with the rotated element.

  y1                   real*4 | real*8
                      On entry, the second element y(1) of the input vector.
                      On exit, y1 is unchanged.

  param                real*4 | real*8
                      On entry, param is unspecified.
                      On exit, param contains an array defining the type of
                      transform matrix H constructed:

                      PARAM(1) specifies the flag characteristic: -1.0, 0.0,
                      1.0, -2.0

                      PARAM(2) specifies H(11) value

                      PARAM(3) specifies H(21) value

                      PARAM(4) specifies H(12) value

                      PARAM(5) specifies H(22) value

Returns {S,D}rotmg object with accessors:

d1                   real*4 | real*8
                    d1 updated.

d2                   real*4 | real*8
                    d2 is updated.

x1                   real*4 | real*8
                    the rotated element.

y1                   real*4 | real*8
                    the second element y(1) of the input vector.

param                real*4 | real*8
                     an array defining the type of
                    transform matrix H constructed:

                    param[1] specifies the flag characteristic: -1.0, 0.0, 1.0, -2.0

                    param[2] specifies H(11) value

                    param[3] specifies H(21) value

                    param[4] specifies H(12) value

                    param[5] specifies H(22) value

Description
  The _ROTMG subroutines construct a modified Givens transform that elim-
  inates the second element of a two-element vector and can be used to intro-
  duce zeros selectively into a matrix.  These routine use the modification
  due to Gentleman of the Givens plane rotations. This modification elim-
  inates the square root from the construction of the plane rotation and
  reduces the operation count when the modified Givens rotation, rather than
  the standard Givens rotations are applied.  In most applications, the
  scale factors d1 and d2 are initially set to 1 and then modified by _ROTMG
  as necessary.

  Given real a and b in factored form:

       a = sqrt(d1) * x(1)
       b = sqrt(d2) * y(1)

  SROTMG and DROTMG construct the modified Givens plane rotation, d1', d2'
  and  H as

                H(11)  H(12)
                H(21)  H(22)

  such that

         |-                 -|        |-    -|          |-   -|    |-   -|
         |sqrt(d1')    0     |        | x(1) |          |  a  |    |  r  |
         |                   | * H *  |      |   =  G * |     | =  |     |
         |  0      sqrt(d2') |        | y(1) |          |  b  |    |  0  |
         |_                 _|        |_    _|          |_   _|    |_   _|

  where G is a 2 by 2 Givens plane rotation matrix which annihilates b, and
  where H is chosen for numerical stability and computational efficiency.

  The routine _ROTM  applies the matrix H, as constructed by _ROTMG, to a
  pair of real vectors, x and y, each with n elements, as follows:

                   |-    -|         |-    -|
                   | x(i) |         | x(i) |
                   |      |   = H * |      |
                   | y(i) |         | y(i) |
                   |_    _|         |_    _|

  These vectors may be either rows or columns of matrices and the indexing of
  the vectors may be either forwards or backwards.

  Depending on the value of IPARAM(1), the matrix H is defined as follows:


  PARAM(1)= -1.0
         H(11)  H(12)
         H(21)  H(22)

  PARAM(1)= 0.0
           1.0  H(12)
         H(21)    1.0

  PARAM(1)= 1.0
         H(11)    1.0
          -1.0  H(22)

  PARAM(1)= -2.0
           1.0    0.0
           0.0    1.0

  The routines _ROTMG and _ROTM perform similar tasks to the routines _ROTG
  and _ROT, which construct and apply the standard Givens plane rotations.
  The modified Givens rotations reduce the operation count of constructing
  and applying the rotations at the cost of increased storage to represent
  the rotations.
*/

#include <ruby.h>
#include "rb_blas.h"


static VALUE myClass;


static void rblas_drotmg_mark(dRotmg *drotmg)
{
  //Empty
}

static void rblas_drotmg_free(dRotmg *drotmg)
{
  if(drotmg)
  {
    free(drotmg);
  }
}

static VALUE rblas_drotmg_alloc(VALUE klass)
{
  dRotmg *drotmg = calloc(1, sizeof(dRotmg));
  return Data_Wrap_Struct(klass, rblas_drotmg_mark, rblas_drotmg_free, drotmg); //wraps a c struct in VALUE
  //recover with Data_Get_Struct()
}

static VALUE rblas_drotmg_initialize(int argc, VALUE *argv, VALUE obj)
{
  dRotmg *drotmg;
  VALUE d1, d2, x1, y1;


  rb_scan_args(argc, argv, "4", &d1, &d2, &x1, &y1);

  Data_Get_Struct(obj, dRotmg, drotmg);
  drotmg->class_id = myClass; //Remember how we came into this world.

  drotmg->d1 = NUM2DBL(d1);
  drotmg->d2 = NUM2DBL(d2);
  drotmg->x1 = NUM2DBL(x1);
  drotmg->y1 = NUM2DBL(y1);
  cblas_drotmg(&drotmg->d1, &drotmg->d2, &drotmg->x1, drotmg->y1, drotmg->param);

  return obj;
}

static VALUE rblas_drotmg_get_d1(VALUE obj)
{
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);
  return rb_float_new((double)drotmg->d1 );
}
static VALUE rblas_drotmg_set_d1(VALUE obj, VALUE d1)
{
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);
  drotmg->d1 = NUM2DBL(d1);
  return d1;
}

static VALUE rblas_drotmg_get_d2(VALUE obj)
{
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);
  return rb_float_new((double)drotmg->d2 );
}
static VALUE rblas_drotmg_set_d2(VALUE obj, VALUE d2)
{
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);
  drotmg->d2 = NUM2DBL(d2);
  return d2;
}


static VALUE rblas_drotmg_get_x1(VALUE obj)
{
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);
  return rb_float_new((double)drotmg->x1 );
}
static VALUE rblas_drotmg_set_x1(VALUE obj, VALUE x1)
{
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);
  drotmg->x1 = NUM2DBL(x1);
  return x1;
}

static VALUE rblas_drotmg_get_y1(VALUE obj)
{
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);
  return rb_float_new((double)drotmg->y1 );
}
static VALUE rblas_drotmg_set_y1(VALUE obj, VALUE y1)
{
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);
  drotmg->y1 = NUM2DBL(y1);
  return y1;
}


//Creates a drotmg class,
//if c_array is not NULL, then it sets the data field of the Blas class to a copy of c_array.
VALUE rblas_drotmg_new_instance( VALUE self, double a, double b, double c, double s)
{
  char argc = 4;
  VALUE argv[] = { rb_float_new(a), rb_float_new(b),  rb_float_new(c), rb_float_new(s)  };

  return  rb_class_new_instance( argc, argv, self); //new instance of this class.
}

static VALUE rblas_drotmg_get_param(VALUE obj)
{
  int i;
  VALUE param_array;
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);  
  
  param_array = rb_ary_new2(5); 
  for(i=0; i<5; i++)
  {
    rb_ary_store(param_array, i, rb_float_new(drotmg->param[i]));
  }
  
  return param_array;
}

static VALUE rblas_drotmg_to_a(VALUE obj)
{
  int i;
  VALUE result_array, param_array;
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);  
  result_array = rb_ary_new2(5); 
  rb_ary_store(result_array, 0, rb_float_new(drotmg->d1)); 
  rb_ary_store(result_array, 1, rb_float_new(drotmg->d2)); 
  rb_ary_store(result_array, 2, rb_float_new(drotmg->x1));
  rb_ary_store(result_array, 3, rb_float_new(drotmg->y1));
  
  param_array = rb_ary_new2(5); 
  for(i=0; i<5; i++)
  {
    rb_ary_store(param_array, i, rb_float_new(drotmg->param[i]));
  }
  rb_ary_store(result_array, 4, param_array);
  
  return result_array;
}

static VALUE rblas_drotmg_to_s(VALUE obj)
{
  char buffer[16];
  VALUE s;
  int i;
  dRotmg *drotmg;

  Data_Get_Struct(obj, dRotmg, drotmg);

  s = rb_str_new2("[ ");

  sprintf(buffer, "%E", drotmg->d1);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );

  sprintf(buffer, "%E", drotmg->d2);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );

  sprintf(buffer, "%E", drotmg->x1);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );

  sprintf(buffer, "%E", drotmg->y1);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );

  s = rb_str_new2("[ ");
  sprintf(buffer, "%E", drotmg->param[0]);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", [", 2 );    

  for(i=1; i<3; i++)
  {
    sprintf(buffer, "%E", drotmg->param[i]);  
    rb_str_cat(s, buffer, strlen(buffer) );
    rb_str_cat(s, ", ", 2 );    
  }
  sprintf(buffer, "%E", drotmg->param[i]);  
  rb_str_cat(s, buffer, strlen(buffer) );
  
  rb_str_cat(s, "]] ]", 1 );  

  return s;
}


void Init_drotmg(VALUE myModule, VALUE parent_class) 
{
  // Create a Ruby class in this module.
  // rb_cObject is defined in ruby.h

  myClass = rb_define_class_under(myModule, "Drotmg", parent_class);
  rb_define_alloc_func(myClass, rblas_drotmg_alloc); 
  rb_define_method(myClass, "initialize", rblas_drotmg_initialize, -1);
  //rb_define_method(myClass, "initialize_copy", rblas_init_copy, 1);

  rb_define_method(myClass, "d1", rblas_drotmg_get_d1, 0); 
  rb_define_method(myClass, "d1=", rblas_drotmg_set_d1, 1); 
  rb_define_method(myClass, "d2", rblas_drotmg_get_d2, 0); 
  rb_define_method(myClass, "d2=", rblas_drotmg_set_d2, 1); 
  rb_define_method(myClass, "x1", rblas_drotmg_get_x1, 0); 
  rb_define_method(myClass, "x1=", rblas_drotmg_set_x1, 1); 
  rb_define_method(myClass, "y1", rblas_drotmg_get_y1, 0); 
  rb_define_method(myClass, "y=", rblas_drotmg_set_y1, 1); 
  rb_define_method(myClass, "param", rblas_drotmg_get_param, 0); 

  rb_define_method(myClass, "to_a", rblas_drotmg_to_a, 0);
  rb_define_method(myClass, "to_s", rblas_drotmg_to_s, 0);

}

