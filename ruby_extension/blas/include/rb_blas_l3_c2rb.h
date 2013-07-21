#ifndef RB_BLAS_l3_C2RB_H
#define RB_BLAS_l3_C2RB_H

void Sub_Init_blas_l3(VALUE myClass);

VALUE rb_blas_xgemm(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xhemm(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xher2k(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xherk(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xsymm(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xsyr2k(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xsyrk(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xtrmm(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xtrsm(int argc, VALUE *argv, VALUE self);

#endif //RB_BLAS_l3_C2RB_H