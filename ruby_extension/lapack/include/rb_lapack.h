#ifndef RB_LAPACK_H
#define RB_LAPACK_H
#include <cblas.h>
#include <clapack.h>
#include "rb_blas.h"
#include "rb_lapack_c2rb.h"

VALUE Init_lapack(VALUE myModule, VALUE parent_class);
void Init_SingleLapack(VALUE myModule, VALUE parent_class);
void Init_DoubleLapack(VALUE myModule, VALUE parent_class);

#endif
