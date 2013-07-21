#ifndef LAPACK_C2RB_H
#define LAPACK_C2RB_H
#include <clapack.h>

void Sub_Init_lapack(VALUE myClass);

void Sub_Init_lapack(VALUE myClass);

VALUE rb_lapack_determinant(VALUE self);

VALUE rb_lapack_xgesv_mod(VALUE self, VALUE result, VALUE ipiv);
VALUE rb_lapack_xgesv( VALUE self, VALUE result, VALUE ipiv );

VALUE rb_lapack_xgetrf_mod(VALUE self, VALUE ipiv);

VALUE rb_lapack_xgetrs_mod(VALUE self,  VALUE result, VALUE ipiv, VALUE transpose );
VALUE rb_lapack_xgetrs( VALUE self, VALUE result, VALUE ipiv, VALUE transpose );

#endif