#ifndef RB_BLAS_EXTRAS_H
#define RB_BLAS_EXTRAS_H

VALUE rb_blas_abs(VALUE self);
VALUE rb_blas_eql(VALUE self, VALUE matrix_y);
VALUE rb_blas_within_bound(VALUE self, VALUE expected, VALUE bound);
VALUE rb_blas_times(VALUE self, VALUE multiplier);

#endif //RB_BLAS_EXTRAS_H