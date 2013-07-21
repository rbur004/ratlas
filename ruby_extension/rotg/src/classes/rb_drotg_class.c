#include "rb_rotg.h"

/*
FORMAT
  {S,D,C,Z}rotg.new(a, b)

Arguments

  a                    real*4 | real*8 | complex*8 | complex*16
                       the first element of the input vector.

  b                    real*4 | real*8 | complex*8 | complex*16
                       the second element of the input vector.  

                      
Returns {S,D,C,Z}rotg object with
r                     real*4 | real*8 | complex*8 | complex*16
                      the rotated element r.

z                     real*4 | real*8 | complex*8 | complex*16
                      SROTG and DROTG, the reconstruction element z.  
                      For CROTG and ZROTG, b.
                      
c                     real*4 | real*8
                      the first rotation element, that is, the cosine of the angle of rotation.
                      
s                     real*4 | real*8 | complex*8 | complex*16
                      the second rotation element, that is, the sine of the angle of rotation.


Description
  The _ROTG subroutines construct a Givens plane rotation that eliminates the
  second element of a two-element vector and can be used to introduce zeros
  selectively into a matrix.

  Using a and b to represent elements of an input real vector, the Srotg and
  Drotg functions calculate the elements c and s of  an orthogonal matrix
  such that:

   c*a + s*b = r
  -s*a + c*b = 0

  Using a and b to represent elements of an input complex vector, the Crotg
  and Zrotg functions calculate the elements real c and complex s of an
  orthogonal matrix such that:

              c*a + s*b = r
  -conjugate(s)*a + c*b = 0

  A real Givens plane rotation is constructed for values a and b by computing
  values for r, c, s, and z, as follows:

  r=p * (a**(2)+b**(2))**(1/2)

  p = SIGN(a)  if  |a| > |b|
  p = SIGN(b)  if  |a|<=|b|

  c = a/r if r is not equal to 0
  c = 1 if r = 0

  s = b/r if r is not equal to 0
  s = 0 if r = 0

  z = s if |a| > |b|
  z = 1/c if |a|<=|b|, c is not      equal to 0, and r is not equal to 0.
  z = 1 if |a|<=|b|, c = 0, and      r is not equal to 0.
  z = 0 if r = 0

  SROTG and DROTG can use the reconstruction element z to store the rotation
  elements for future use. The quantities c and s are reconstructed from z as
  follows:

  For |z| = 1, c = 0.0  and  s = 1.0

  For |z| < 1, c = (1-z**(2))**(1/2) and s = z

  For |z| > 1, c = 1/z and s = (1-c**(2))**(1/2)

  A complex Givens plane rotation is constructed for values a and b by com-
  puting values for real c, complex s and complex r, as follows:

  p=(|a|**(2)+|b|**(2))**(1/2)

  q = a/|a|

  r = qp  if  |a| is not equal to 0.
  r = b  if  |a| is equal to 0.

  c = |a|/p if |a| is not equal to 0
  c = 0 if |a| is equal to 0

  s = q*conjugate(b)/p if |a| is not equal to 0
  s = (1.0,0.0) if |a| is equal to 0

  The absolute value used in the above definitions corresponds to the strict
  definition of the absolute value of a complex number.

  The arguments c and s are passed to the _ROT subroutines.
*/
#include <ruby.h>
#include "rb_blas.h"


static VALUE myClass;


static void rblas_drotg_mark(dRotg *drotg)
{
  //Empty
}

static void rblas_drotg_free(dRotg *drotg)
{
  if(drotg)
  {
    free(drotg);
  }
}

static VALUE rblas_drotg_alloc(VALUE klass)
{
  dRotg *drotg = calloc(1, sizeof(dRotg));
  return Data_Wrap_Struct(klass, rblas_drotg_mark, rblas_drotg_free, drotg); //wraps a c struct in VALUE
  //recover with Data_Get_Struct()
}

static VALUE rblas_drotg_initialize(int argc, VALUE *argv, VALUE obj)
{
  dRotg *drotg;
  VALUE a, b, c, s;


  rb_scan_args(argc, argv, "22", &a, &b, &c, &s);

  Data_Get_Struct(obj, dRotg, drotg);
  drotg->class_id = myClass; //Remember how we came into this world.

  if(a == Qnil) drotg->r = 0.0; else drotg->r = NUM2DBL(a);
  if(b == Qnil) drotg->z = 0.0; else drotg->z = NUM2DBL(b);
  if(c == Qnil) drotg->c = 0.0; else drotg->c = NUM2DBL(c);
  if(s == Qnil) drotg->z = 0.0; else drotg->s = NUM2DBL(s);
  if(argc == 2) //if not nil then assume have to calculte values for c & s
    cblas_drotg(&drotg->r, &drotg->z, &drotg->c, &drotg->s);

  return obj;
}

static VALUE rblas_drotg_get_r(VALUE obj)
{
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);
  return rb_float_new( drotg->r );
}
static VALUE rblas_drotg_set_r(VALUE obj, VALUE r)
{
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);
  drotg->r = NUM2DBL(r);
  return r;
}

static VALUE rblas_drotg_get_z(VALUE obj)
{
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);
  return rb_float_new( drotg->z );
}
static VALUE rblas_drotg_set_z(VALUE obj, VALUE z)
{
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);
  drotg->z = NUM2DBL(z);
  return z;
}

static VALUE rblas_drotg_get_c(VALUE obj)
{
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);
  return rb_float_new( drotg->c );
}
static VALUE rblas_drotg_set_c(VALUE obj, VALUE c)
{
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);
  drotg->c = NUM2DBL(c);
  return c;
}

static VALUE rblas_drotg_get_s(VALUE obj)
{
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);
  return rb_float_new( drotg->s );
}
static VALUE rblas_drotg_set_s(VALUE obj, VALUE s)
{
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);
  drotg->s = NUM2DBL(s);
  return s;
}


//Creates a drotg class,
//if c_array is not NULL, then it sets the data field of the Blas class to a copy of c_array.
VALUE rblas_drotg_new_instance( VALUE self, double a, double b, double c, double s)
{
  char argc = 4;
  VALUE argv[] = { rb_float_new(a), rb_float_new(b),  rb_float_new(c), rb_float_new(s)  };
  
  return  rb_class_new_instance( argc, argv, self); //new instance of this class.
}

static VALUE rblas_drotg_to_a(VALUE obj)
{
  VALUE result_array;
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);  
  result_array = rb_ary_new2(4); 
  rb_ary_store(result_array, 0, rb_float_new(drotg->r)); 
  rb_ary_store(result_array, 1, rb_float_new(drotg->z)); 
  rb_ary_store(result_array, 2, rb_float_new(drotg->c));
  rb_ary_store(result_array, 3, rb_float_new(drotg->s));
  
  return result_array;
}

static VALUE rblas_drotg_to_s(VALUE obj)
{
  char buffer[16];
  VALUE s;
  dRotg *drotg;

  Data_Get_Struct(obj, dRotg, drotg);
  
  s = rb_str_new2("[ ");
  
  sprintf(buffer, "%E", drotg->r);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );
  
  sprintf(buffer, "%E", drotg->z);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );
  
  sprintf(buffer, "%E", drotg->c);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );
  
  sprintf(buffer, "%E", drotg->s);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, "]", 1 );  
  
  return s;
}

void Init_drotg(VALUE myModule, VALUE parent_class) 
{
  // Create a Ruby class in this module.
  // rb_cObject is defined in ruby.h
  
  myClass = rb_define_class_under(myModule, "Drotg", parent_class);
  rb_define_alloc_func(myClass, rblas_drotg_alloc); 
  rb_define_method(myClass, "initialize", rblas_drotg_initialize, -1);
  //rb_define_method(myClass, "initialize_copy", rblas_init_copy, 1);
  
  rb_define_method(myClass, "r", rblas_drotg_get_r, 0); 
  rb_define_method(myClass, "r=", rblas_drotg_set_r, 1); 
  rb_define_method(myClass, "z", rblas_drotg_get_z, 0); 
  rb_define_method(myClass, "z=", rblas_drotg_set_z, 1); 
  rb_define_method(myClass, "c", rblas_drotg_get_c, 0); 
  rb_define_method(myClass, "c=", rblas_drotg_set_c, 1); 
  rb_define_method(myClass, "s", rblas_drotg_get_s, 0); 
  rb_define_method(myClass, "s=", rblas_drotg_set_s, 1); 
  
  rb_define_method(myClass, "to_a", rblas_drotg_to_a, 0);
  rb_define_method(myClass, "to_s", rblas_drotg_to_s, 0);
}

