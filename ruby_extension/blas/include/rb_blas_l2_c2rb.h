#ifndef RB_BLAS_l2_C2RB_H
#define RB_BLAS_l2_C2RB_H

VALUE rb_blas_xtrmv(int argc, VALUE *argv, VALUE self);
void Sub_Init_blas_l2(VALUE myClass);

VALUE rb_blas_mv(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_r(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_r2(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_rc(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_sv(int argc, VALUE *argv, VALUE self);

#endif //RB_BLAS_l2_C2RB_H