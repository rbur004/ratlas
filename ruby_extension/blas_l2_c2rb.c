#include "rb_blas.h"
#include "rb_lapack.h"
#include "blas_l2_functions/hack.c"

#ifdef XXXXX

/*Ruby will pass in 
  The Blas class (self) already holds the matrix datatype and dimensions, so these don't need to be passed in.
  the first argument is an array of options argument (or nil) for the UPLO, TRANS and DIAG options
  The second argument is the scalar alpha.
  the third and forth arguments are the vector X and incX
  the fourth, fifth and sixth arguments are optional, and are the scalar beta, the vector Y and incY 

  The function to use is selected from the matrix data_type as an index into mv.
  The Vector Y and its associated arguments are not used for operations on TB, TP and TR matrices.

	Purpose   
    =======
    xGEMV and xGBMV performs one of the matrix-vector operations
       y = alpha*A*x + Beta*y 
    or y = alpha*A.transpose*x + Beta*y 
    or y = alpha*A.conjugate_transpose*x + Beta*y
    where x and y are n element vectors, alpha and beta are scalars, and A is a general matrix.
    
    xHEMV, xHBMV, xHPMV, xSYMV, xSBMV and xSPMV performs one of the matrix-vector operations
       y = alpha*A*x + Beta*y
    where x and y are n element vectors, alpha and beta are scalars, and A is a general matrix.
 
    xTRMV, xTBMV and xTPMV  performs one of the matrix-vector operations   
         x := A*x
    or   x := A.transpose*x  
    or   x := A.conjugate_transpose*x 
    where x is an n element vector and  A is an n by n unit, or non-unit,   
    upper or lower triangular matrix.   
*/

JumpTable mv[] = 
{
  {GB, HAS_TRANSA | NEEDS_M | NEEDS_KU | NEEDS_KL,                        sgbmv, dgbmv, cgbmv, zgbmv},
  {GE, HAS_TRANSA |NEEDS_M,                                               sgemv, dgemv, cgemv, zgemv},
  {HB, HAS_UPLO | NEEDS_K,                                                NULL, NULL, chbmv, zhbmv},
  {HE, HAS_UPLO,                                                          NULL, NULL, cgbmv, zgbmv},
  {HP, HAS_UPLO | NO_LDA,                                                 NULL, NULL, cgbmv, zgbmv},
  {SB, HAS_UPLO | NEEDS_K,                                                sgbmv, dgbmv, NULL, NULL},
  {SP, HAS_UPLO | NO_LDA,                                                 sgbmv, dgbmv, NULL, NULL},
  {SY, HAS_UPLO,                                                          sgbmv, dgbmv, NULL, NULL},
  {TB, HAS_UPLO | HAS_TRANSA | HAS_DIAG | NEEDS_K | NO_Y_VECTOR,          sgbmv, dgbmv, cgbmv, zgbmv},
  {TP, HAS_UPLO | HAS_TRANSA | HAS_DIAG | NO_Y_VECTOR | NO_LDA,           sgbmv, dgbmv, cgbmv, zgbmv},
  {TR, HAS_UPLO | HAS_TRANSA | HAS_DIAG | NO_Y_VECTOR,                    sgbmv, dgbmv, cgbmv, zgbmv},
  {US, 0,                                                                 NULL, NULL, NULL, NULL}
};

VALUE rb_blas_mv(int argc, VALUE *argv, VALUE self)
{
  Matrix *da, *dx;
  int incx;
  char *uplo, *trans, *diag;
  int n;
  float c_f[2], s_f[2];
  double c_d[2], s_d[2];
  char error_msg[64];
  VALUE vector_x,  n_value,  incx_value;
  VALUE uplo_value, trans_value, diag_value;

  rb_scan_args(argc, argv, "15", &vector_x, &n_value,  &incx_value, 
                &uplo_value, &trans_value, &diag_value);

  Data_Get_Struct(self, Matrix, da);
  Data_Get_Struct(vector_x, Matrix, dx);

  if(incx_value == Qnil)
    incx = 1;
  else
    incx = NUM2INT(incx_value);

  if(n_value == Qnil)
    n = da->nrows;
  else
    n = NUM2INT(n_value);

  if(dx == NULL || dx->nrows != 1)
  { sprintf(error_msg, "x is not a Vector");
    rb_raise(rb_eRuntimeError, error_msg);
  }

  if(n > dx->nrows / incx )
  { sprintf(error_msg, "Too few elements inx for incx given");
    rb_raise(rb_eRuntimeError, error_msg);
  }
  
  if(uplo_value == Qnil)
    uplo = "L";
  else
    uplo = NUM2INT(uplo_value) == CblasUpper ? "U":"L";
    
  if(diag_value == Qnil)
    diag = "N";
  else
    diag = NUM2INT(diag_value) == CblasUnit ? "U":"N";
    
  if(trans_value == Qnil)
    trans = 'N';
  else
    trans = NUM2INT(trans_value) == CblasTrans ? "T": 
              (NUM2INT(trans_value) == CblasConjTrans ? "C":"N");
  

  switch(dx->data_type)
  {
  case Single_t: //s
    cblas_strmv(n , (float *)dx->data, incx, (float *)dy->data, incy, c_f[0], s_f[0]); 
    break;
  case Double_t: //d
    if(c_value == Qnil)
      c_d[0] = 1.0;
    else
      c_d[0] = NUM2DBL(c_value);
    if(s_value == Qnil)
      s_d[0] = 1.0;
    else
      s_d[0] = NUM2DBL(s_value);
    cblas_dtrmv(n , (double *)dx->data, incx, (double *)dy->data, incy, c_d[0], s_d[0]); 
    break;
  case Complex_t: //c
    if(c_value == Qnil)
    {
      c_f[0] = 1.0;
      c_f[1] = 1.0;
    }
    else
    {
      c_f[0] = NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, c_value) );
      c_f[1] = NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, c_value ) );
    }
    if(s_value == Qnil)
    {
      s_f[0] = 1.0;
      s_f[1] = 1.0;
    }
    else
    {
      s_f[0] = NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, s_value) );
      s_f[1] = NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, s_value ) );
    }
    cblas_crot(n , dx->data, incx, dy->data, incy, c_f[0], s_f[0]); 
    break;
  case Double_Complex_t: //z
    if(c_value == Qnil)
    {
      c_d[0] = 1.0;
      c_d[1] = 1.0;
    }
    else
    {
      c_d[0] = NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, c_value) );
      c_d[1] = NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, c_value ) );
    }
    if(s_value == Qnil)
    {
      s_d[0] = 1.0;
      s_d[1] = 1.0;
    }
    else
    {
      s_d[0] = NUM2DBL(rb_funcall( rb_intern("Complex"),  rb_intern("real"),  1, s_value) );
      s_d[1] = NUM2DBL(rb_funcall(rb_intern("Complex"),  rb_intern("image"),  1, s_value ) );
    }
    cblas_zrot(n , dx->data, incx, dy->data, incy, c_d[0], s_d[0]); 
    break;
  default:
    sprintf(error_msg, "Invalid data_type (%d) in Matrix", dx->data_type);
    rb_raise(rb_eRuntimeError, error_msg);
    break; //Never reaches here.
  }

  return self;
}

/*Ruby will pass in 
  The Blas class (self) already holds the matrix datatype and dimensions, so these don't need to be passed in.
  the first argument is an array of options argument (or nil) for the UPLO, TRANS and DIAG options
  the second and third arguments are the vector X and incX

  The function to use is selected from the matrix data_type as an index into sv.
 
    xTRSV, xTBSV and xTPSV  performs one of the matrix-vector operations   
         x := A.inverse*x
    or   x := A.transpose.inverse*x  
    or   x := A.conjugate_transpose.inverse*x 
    where x is an n element vector and  A is an n by n unit, or non-unit,   
    upper or lower triangular matrix.   
*/
JumpTable sv[] = 
{
  {GB, 0,                                                         NULL, NULL, NULL, NULL},
  {GE, 0,                                                         NULL, NULL, NULL, NULL},
  {HB, 0,                                                         NULL, NULL, NULL, NULL},
  {HE, 0,                                                         NULL, NULL, NULL, NULL},
  {HP, 0,                                                         NULL, NULL, NULL, NULL},
  {SB, 0,                                                         NULL, NULL, NULL, NULL},
  {SP, 0,                                                         NULL, NULL, NULL, NULL},
  {SY, 0,                                                         NULL, NULL, NULL, NULL},
  {TB, HAS_UPLO | HAS_TRANSA | HAS_DIAG | NEEDS_K | NO_Y_VECTOR,  stbsv, dtbsv, ctbsv, zgbsv},
  {TP, HAS_UPLO | HAS_TRANSA | HAS_DIAG | NO_Y_VECTOR | NO_LDA,   stpsv, dtpsv, ctpsv, ztpsv},
  {TR, HAS_UPLO | HAS_TRANSA | HAS_DIAG | NO_Y_VECTOR,            strsv, dtrsv, ctrsv, ztrsv},
  {US, 0,                                                         NULL, NULL, NULL, NULL}
};

VALUE rb_blas_sv(int argc, VALUE *argv, VALUE self)
{
  return self;
}

/*Ruby will pass in 
  The Blas class (self) already holds the matrix datatype and dimensions, so these don't need to be passed in.
  the first argument is an array of options argument (or nil) for the UPLO, TRANS and DIAG options
  the second and third arguments are the vector X and incX
  the forth and fifth arguments are the vector Y and incY which aren't used by H and S matrices

  The function to use is selected from the matrix data_type as an index into sv.
 
    xGER, xGERU  performs one of the matrix-vector operations   
         x := alpha*X*Y.transpose + A
    where x and y are n element vectors 

    xHER, xHPR performs one of the matrix-vector operations   
         x := alpha*X*X.conjugate_transpose + A

    xSPR, xSYR performs one of the matrix-vector operations   
         x := alpha*X*X.transpose + A
*/
JumpTable r[] = 
{
  {GB, 0,                                                         NULL, NULL, NULL, NULL},
  {GE, NEEDS_M,                                                   sger, dger, cgeru, zgeru},
  {HB, 0,                                                         NULL, NULL, NULL, NULL},
  {HE, UPLO | NO_Y_VECTOR,                                        NULL, NULL, cher, zher},
  {HP, UPLO | NO_Y_VECTOR,                                        NULL, NULL, chpr, zhpr},
  {SB, 0,                                                         NULL, NULL, NULL, NULL},
  {SP, UPLO | NO_Y_VECTOR,                                        cspr, zspr, NULL, NULL},
  {SY, UPLO | NO_Y_VECTOR,                                        csyr, zsyr, NULL, NULL},
  {TB, 0,                                                         NULL, NULL, NULL, NULL},
  {TP, 0,                                                         NULL, NULL, NULL, NULL},
  {TR, 0,                                                         NULL, NULL, NULL, NULL},
  {US, 0,                                                         NULL, NULL, NULL, NULL}
};

VALUE rb_blas_r(int argc, VALUE *argv, VALUE self)
{
  return self;
}

/*Ruby will pass in 
  The Blas class (self) already holds the matrix datatype and dimensions, so these don't need to be passed in.
  the first argument is an array of options argument (or nil) for the UPLO, TRANS and DIAG options
  the second and third arguments are the vector X and incX
  the forth and fifth arguments are the vector Y and incY which aren't used by H and S matrices

  The function to use is selected from the matrix data_type as an index into sv.
 
    xGERC  performs one of the matrix-vector operations   
         x := alpha*X*Y.conjugate_transpose + A
    where x and y are n element vectors 

*/
JumpTable rc[] = 
{
  {GB, 0,                                       NULL, NULL, NULL, NULL},
  {GE, NEEDS_M,                                 NULL, NULL, cgerc, zgerc},
  {HB, 0,                                       NULL, NULL, NULL, NULL},
  {HE, 0,                                       NULL, NULL, NULL, NULL},
  {HP, 0,                                       NULL, NULL, NULL, NULL},
  {SB, 0,                                       NULL, NULL, NULL, NULL},
  {SP, 0,                                       NULL, NULL, NULL, NULL},
  {SY, 0,                                       NULL, NULL, NULL, NULL},
  {TB, 0,                                       NULL, NULL, NULL, NULL},
  {TP, 0,                                       NULL, NULL, NULL, NULL},
  {TR, 0,                                       NULL, NULL, NULL, NULL},
  {US, 0,                                       NULL, NULL, NULL, NULL}
};

VALUE rb_blas_rc(int argc, VALUE *argv, VALUE self)
{
  return self;
}

/*Ruby will pass in 
  The Blas class (self) already holds the matrix datatype and dimensions, so these don't need to be passed in.
  the first argument is an array of options argument (or nil) for the UPLO, TRANS and DIAG options
  the second and third arguments are the vector X and incX
  the forth and fifth arguments are the vector Y and incY which aren't used by H and S matrices

  The function to use is selected from the matrix data_type as an index into sv.
 
    x and y are n element vectors 

    xHER2, xHPR2 performs one of the matrix-vector operations   
         x := alpha*X*X.conjugate_transpose + y*(alpha*x).conjugate_transpose + A

    xSPR2, xSYR2 performs one of the matrix-vector operations   
         x := alpha*X*X.transpose + alpha*y*x.transpose + A
*/

JumpTable r2[] = 
{
  {GB, 0,                                           NULL, NULL, NULL, NULL},
  {GE, 0,                                           NULL, NULL, NULL, NULL},
  {HB, 0,                                           NULL, NULL, NULL, NULL},
  {HE, UPLO,                                        NULL, NULL, cher2, zher2},
  {HP, UPLO,                                        NULL, NULL, chpr2, zhpr2},
  {SB, 0,                                           NULL, NULL, NULL, NULL},
  {SP, UPLO,                                        cspr2, zspr2, NULL, NULL},
  {SY, UPLO,                                        csyr2, zsyr2, NULL, NULL},
  {TB, 0,                                           NULL, NULL, NULL, NULL},
  {TP, 0,                                           NULL, NULL, NULL, NULL},
  {TR, 0,                                           NULL, NULL, NULL, NULL},
  {US, 0,                                           NULL, NULL, NULL, NULL}
};

VALUE rb_blas_r2(int argc, VALUE *argv, VALUE self)
{
  return self;
}

#endif xxxx

void Sub_Init_blas_l2(VALUE myClass)
{
//"gemv", "gbmv", "hemv","hbmv", "symv", "sbmv", "spmv", "trmv", "tbmv", "tpmv"
  rb_define_method(myClass, "mv", rb_blas_mv, -1);

//"trsv", "tbsv", "tpsv"
  rb_define_method(myClass, "sv", rb_blas_sv, -1);

// "ger", "geru", "her", "hpr", "syr", "spr" 
  rb_define_method(myClass, "r", rb_blas_r, -1);

//"gerc"
  rb_define_method(myClass, "rc", rb_blas_rc, -1);
  
// "her2", "hpr2", "syr2", "spr2"
  rb_define_method(myClass, "r2", rb_blas_r2, -1);
}
