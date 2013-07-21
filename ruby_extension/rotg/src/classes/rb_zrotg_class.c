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


static void rblas_zrotg_mark(zRotg *zrotg)
{
  //Empty
}

static void rblas_zrotg_free(zRotg *zrotg)
{
  if(zrotg)
  {
    free(zrotg);
  }
}

static VALUE rblas_zrotg_alloc(VALUE klass)
{
  zRotg *zrotg = calloc(1, sizeof(zRotg));
  return Data_Wrap_Struct(klass, rblas_zrotg_mark, rblas_zrotg_free, zrotg); //wraps a c struct in VALUE
  //recover with Data_Get_Struct()
}

static VALUE rblas_zrotg_initialize(int argc, VALUE *argv, VALUE obj)
{
  zRotg *zrotg;
  VALUE a, b, c, s;
  
  
  rb_scan_args(argc, argv, "22", &a, &b, &c, &s);
  
  Data_Get_Struct(obj, zRotg, zrotg);
  zrotg->class_id = myClass; //Remember how we came into this world.

  zrotg->r[0] = NUM2DBL(rb_funcall( a,  rb_intern("real"),  0 ) );
  zrotg->r[1] = NUM2DBL(rb_funcall( a,  rb_intern("image"),  0 ) );
  zrotg->z[0] = NUM2DBL(rb_funcall( b,  rb_intern("real"),  0 ) );
  zrotg->z[1] = NUM2DBL(rb_funcall( b,  rb_intern("image"),  0 ) );
  if(c == Qnil) zrotg->c = 0; else zrotg->c = NUM2DBL(c);
  if(s == Qnil) 
  {
		zrotg->s[0] = 0.0;
		zrotg->s[1] = 0.0; 
	}
	else
	{
	  zrotg->s[0] = NUM2DBL(rb_funcall( s,  rb_intern("real"),  0 ) );
	  zrotg->s[1] = NUM2DBL(rb_funcall( s,  rb_intern("image"),  0 ) );
	}
	if( argc == 2) //Only calculate the c & s values if they are not given
  	cblas_zrotg(zrotg->r, zrotg->z, &zrotg->c, zrotg->s);

  return obj;
}

static VALUE rblas_zrotg_get_r(VALUE obj)
{
  zRotg *zrotg;
  VALUE complex_result[2];
  
  Data_Get_Struct(obj, zRotg, zrotg);
  complex_result[0] = rb_float_new(zrotg->r[0]);
  complex_result[1] = rb_float_new(zrotg->r[1]);
  
//? (void) rb_require( "complex" );
//?  (void) rb_require( "Math" );
  return rb_class_new_instance(2, complex_result, rb_const_get(rb_mMath, rb_intern("Complex")));
}
static VALUE rblas_zrotg_set_r(VALUE obj, VALUE r)
{
  zRotg *zrotg;

  Data_Get_Struct(obj, zRotg, zrotg);
  zrotg->r[0] = NUM2DBL(rb_funcall( r,  rb_intern("real"),  0 ) );
  zrotg->r[1] = NUM2DBL(rb_funcall( r,  rb_intern("image"),  0 ) );
  return r;
}

static VALUE rblas_zrotg_get_z(VALUE obj)
{
    zRotg *zrotg;
    VALUE complex_result[2];

    Data_Get_Struct(obj, zRotg, zrotg);
    complex_result[0] = rb_float_new(zrotg->z[0]);
    complex_result[1] = rb_float_new(zrotg->z[1]);

  //? (void) rb_require( "complex" );
  //?  (void) rb_require( "Math" );
    return rb_class_new_instance(2, complex_result, rb_const_get(rb_mMath, rb_intern("Complex")));
}
static VALUE rblas_zrotg_set_z(VALUE obj, VALUE z)
{
  zRotg *zrotg;

  Data_Get_Struct(obj, zRotg, zrotg);
  zrotg->z[0] = NUM2DBL(rb_funcall( z,  rb_intern("real"),  0 ) );
  zrotg->z[1] = NUM2DBL(rb_funcall( z,  rb_intern("image"),  0 ) );
  return z;
}

static VALUE rblas_zrotg_get_c(VALUE obj)
{
  zRotg *zrotg;

  Data_Get_Struct(obj, zRotg, zrotg);
  return rb_float_new((double)zrotg->c );
}
static VALUE rblas_zrotg_set_c(VALUE obj, VALUE c)
{
  zRotg *zrotg;
 
  Data_Get_Struct(obj, zRotg, zrotg);
  zrotg->c = NUM2DBL(c);
  return c;
}

static VALUE rblas_zrotg_get_s(VALUE obj)
{
    zRotg *zrotg;
    VALUE complex_result[2];

    Data_Get_Struct(obj, zRotg, zrotg);
    complex_result[0] = rb_float_new(zrotg->s[0]);
    complex_result[1] = rb_float_new(zrotg->s[1]);

  //? (void) rb_require( "complex" );
  //?  (void) rb_require( "Math" );
    return rb_class_new_instance(2, complex_result, rb_const_get(rb_mMath, rb_intern("Complex")));
}
static VALUE rblas_zrotg_set_s(VALUE obj, VALUE s)
{
  zRotg *zrotg;

  Data_Get_Struct(obj, zRotg, zrotg);
  zrotg->s[0] = NUM2DBL(rb_funcall( s,  rb_intern("real"),  0 ) );
  zrotg->s[1] = NUM2DBL(rb_funcall( s,  rb_intern("image"),  0 ) );
  return s;
}

static VALUE rblas_zrotg_to_a(VALUE obj)
{
  VALUE result_array;
  
  result_array = rb_ary_new2(4); 
  rb_ary_store(result_array, 0, rblas_zrotg_get_r(obj)); 
  rb_ary_store(result_array, 1, rblas_zrotg_get_z(obj)); 
  rb_ary_store(result_array, 2, rblas_zrotg_get_c(obj));
  rb_ary_store(result_array, 3, rblas_zrotg_get_s(obj));
  
  return result_array;
}

static char * z_to_s(char *buffer, double *z)
{
  if(z[1] < 0)
    sprintf(buffer, "%E - i%E", z[0],-z[1]);
  else
    sprintf(buffer, "%E + i%E", z[0],z[1]);
  return buffer;
}

static VALUE rblas_zrotg_to_s(VALUE obj)
{
  char buffer[16];
  VALUE s;
  zRotg *zrotg;

  Data_Get_Struct(obj, zRotg, zrotg);
  
  s = rb_str_new2("[ ");
  
  z_to_s(buffer, zrotg->r);
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );
  
  z_to_s(buffer, zrotg->z);
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );
  
  sprintf(buffer, "%E", zrotg->c);  
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, ", ", 2 );
  
  z_to_s(buffer, zrotg->s);
  rb_str_cat(s, buffer, strlen(buffer) );
  rb_str_cat(s, "]", 1 );  
  
  return s;
}

void Init_zrotg(VALUE myModule, VALUE parent_class) 
{
  // Create a Ruby class in this module.
  // rb_cObject is defined in ruby.h
  
  myClass = rb_define_class_under(myModule, "Zrotg", parent_class);
  rb_define_alloc_func(myClass, rblas_zrotg_alloc); 
  rb_define_method(myClass, "initialize", rblas_zrotg_initialize, -1);
  //rb_define_method(myClass, "initialize_copy", rblas_init_copy, 1);
  
  rb_define_method(myClass, "r", rblas_zrotg_get_r, 0); 
  rb_define_method(myClass, "r=", rblas_zrotg_set_r, 1); 
  rb_define_method(myClass, "z", rblas_zrotg_get_z, 0); 
  rb_define_method(myClass, "z=", rblas_zrotg_set_z, 1); 
  rb_define_method(myClass, "c", rblas_zrotg_get_c, 0); 
  rb_define_method(myClass, "c=", rblas_zrotg_set_c, 1); 
  rb_define_method(myClass, "s", rblas_zrotg_get_s, 0); 
  rb_define_method(myClass, "s=", rblas_zrotg_set_s, 1); 
  
  rb_define_method(myClass, "to_a", rblas_zrotg_to_a, 0);
  rb_define_method(myClass, "to_s", rblas_zrotg_to_s, 0);
  
}

