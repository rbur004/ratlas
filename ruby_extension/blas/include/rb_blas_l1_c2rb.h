#ifndef RB_BLAS_l1_C2RB_H
#define RB_BLAS_l1_C2RB_H

void Sub_Init_blas_l1(VALUE myClass);

VALUE rb_atlas_xset_mod(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_ixamax(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xasum(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xaxpy(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xaxpy_add(VALUE self, VALUE vector_y);
VALUE rb_blas_xaxpy_mod(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xaxpy_sub(VALUE self, VALUE vector_y);
VALUE rb_blas_xcopy_mod(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xdot(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xdotc(int argc, VALUE *argv, VALUE self);
//rb_blas_xdotsds(int argc, VALUE *argv, VALUE self);
//rb_blas_xdotu(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xnrm2(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xrot_mod(int argc, VALUE *argv, VALUE self);
//rb_blas_xrotg(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xaxpy_mod_actual( int argc, VALUE *argv, VALUE self, int neg);
VALUE rb_blas_xrotm_mod(int argc, VALUE *argv, VALUE self);
//rb_blas_xrotmg(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xscal(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xscal_mod(int argc, VALUE *argv, VALUE self);
VALUE rb_blas_xswap_mod(int argc, VALUE *argv, VALUE self);

#endif //RB_BLAS_l1_C2RB_H