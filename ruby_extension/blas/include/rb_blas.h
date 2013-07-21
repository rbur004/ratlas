#ifndef RB_BLAS_H
#define RB_BLAS_H
#include <ruby.h>
#include <cblas.h>

enum DATA_TYPE   {Double_t=1, Single_t=2, Integer_t=3, Complex_t=4, Double_Complex_t=5 };
#ifdef RB_BLAS_C
const char *DATA_TYPE_Text[] =   {"", "Double_t", "Single_t", "Integer_t", "Complex_t", "Double_Complex_t" };
#else /*RB_BLAS_C*/
extern  const char *DATA_TYPE_Text[];
#endif /*RB_BLAS_C*/

enum MATRIX_TYPE {GB=0, GE, HB, HE, HP, SB, SP, SY, TB, TP, TR, US, LU };
#ifdef RB_BLAS_C
const char *MATRIX_TYPE_Text[] = {"GB", "GE", "HB", "HE", "HP", "SB", "SP", "SY", "TB", "TP", "TR", "US", "LU" };
#else /*RB_BLAS_C*/
extern  const char *MATRIX_TYPE_Text[];
#endif /*RB_BLAS_C*/

enum OPTION_FLAGS { HAS_UPLO=1, HAS_TRANSA=2, HAS_DIAG=4, HAS_SIDE=8, HAS_TRANSB=16, 
                    NEEDS_K=32,  NEEDS_M=64, NEEDS_KU=128, NEEDS_KL=256,
                    NO_Y_VECTOR=512, NO_LDA=1024 };

enum FUNCT_TYPES { F1=HAS_TRANSA | NEEDS_M | NEEDS_KU | NEEDS_KL,
                   F2=HAS_TRANSA |NEEDS_M, 
                   F3=HAS_UPLO | NEEDS_K,
                   F4=HAS_UPLO,
                   F5=HAS_UPLO | NO_LDA,
                   F6=HAS_UPLO | HAS_TRANSA | HAS_DIAG | NEEDS_K | NO_Y_VECTOR,
                   F7=HAS_UPLO | HAS_TRANSA | HAS_DIAG | NO_Y_VECTOR | NO_LDA,
                   F8=HAS_UPLO | HAS_TRANSA | HAS_DIAG | NO_Y_VECTOR
                 };

typedef struct Matrix
{
  VALUE   class_id;     //the ruby class id, which can then be used to create new instances.
  void    *unaligned;   
  void    *data;        //One dimensional array to hold the Matrix. Allocated at creation time. 
  int     nrows;
  int     ncols;
  int     kl;           //Banded Matrix number of subdiagonals
  int     ku;           //Banded Matrix number of superdiagonals
  int     cblas_order;  //CBLAS_ORDER. From cblas.h CBlasRowMajor for Row first and CBlasRowMinor for Columns first.
  int     matrix_type;  //ignored at the moment. Intended to select the fortran lib call
  size_t  data_size;    //Filled in at creation time with sizeof(<data_type>)
  int     data_type;    //From the DATA_TYPE enum above. Used to select the fortran lib call
  struct Matrix  *ipiv;         //If this is an LU matrix, this is the associated pivot matrix.
  void    *last_result; //Bit of a Hack. Used to retain Scalar results in Native format to use in the next call.
} Matrix;


VALUE new_subclass_init_helper(int matrix_type, int data_type, int argc, VALUE *argv, VALUE obj);
VALUE new_row_major_matrix_helper(int matrix_type, int data_type, int argc, VALUE *argv, VALUE obj);
VALUE new_column_major_matrix_helper(int matrix_type, int data_type, int argc, VALUE *argv, VALUE obj);
VALUE Init_blas(VALUE myModule);
VALUE rblas_new_instance( VALUE self, void *c_array, int nrows, int ncols,  int data_type, int matrix_type, int cblas_order);
VALUE rblas_integer_rblas_new_instance(void *data, int nrows, int ncols, int matrix_type, int cblas_order);


void Init_IntegerBlas(VALUE myModule, VALUE parent_class);
void Init_DoubleBlas(VALUE myModule, VALUE parent_class);
void Init_SingleBlas(VALUE myModule, VALUE parent_class);
void Init_ComplexBlas(VALUE myModule, VALUE parent_class);
void Init_DoubleComplexBlas(VALUE myModule, VALUE parent_class);

#include "rb_blas_utils.h"
#include "rb_blas_extras.h"
#include "rb_blas_l1_c2rb.h"
#include "rb_blas_l2_c2rb.h"
#include "rb_blas_l3_c2rb.h"
#include "rb_rotg.h"

#endif

