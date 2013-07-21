#include "rb_blas.h"

//Checks that each member of the matrix is < bound.
//Used for testing of results against known outcome.
VALUE rb_blas_within_bound(VALUE self, VALUE bound)
{
  //Definitely a hack. We want to be able to test the results with the same syntax as scalars
  //i.e. if (error = (result - expected)).abs > bound ...
  Matrix *m;
  int r,c;

  Data_Get_Struct(self, Matrix, m);

  for(r = 0; r < m->nrows; r++)
  {
    for(c = 0; c < m->ncols; c++)
    { //use ruby to do this, as it will happily handle the complex data types.
      if (TYPE(rb_funcall(rb_funcall(rb_blas_member_to_value(m, r, c),  rb_intern("abs"), 0), rb_intern("<="), 1, bound)) == TYPE(Qfalse)) return Qfalse;
    } 
  }
  return Qtrue;
}
