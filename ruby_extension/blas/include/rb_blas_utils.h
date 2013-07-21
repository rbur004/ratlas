#ifndef BLAS_UTIL_H
#define BLAS_UTIL_H

inline void rb_blas_get_member(void *dest, Matrix *m, int row, int col);
inline void rb_blas_set_member(Matrix *m, void *src , int row, int col);
VALUE rb_blas_member_to_value(Matrix *matrix, int row, int col);
void rb_blas_value_to_member(Matrix *matrix, int row, int col, VALUE v);
char * rb_blas_member_to_s(char *buffer, Matrix *matrix, int row, int col);

//void rb_blas_lib_caller(VALUE self, int argc, VALUE *argv, JumpTable *jt);

#endif //BLAS_UTIL_H