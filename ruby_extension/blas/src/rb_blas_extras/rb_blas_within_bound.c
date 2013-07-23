#include "rb_blas.h"

//Checks that each member of the matrix is < bound.
//Used for testing of results against known outcome.
VALUE rb_blas_within_bound(VALUE self, VALUE expected, VALUE bound)
{
  //Definitely a hack. We want to be able to test the results with the same syntax as scalars
  //i.e. if (error = (result - expected)).abs > bound ...
  Matrix *m, *n;
  int r,c;

  Data_Get_Struct(self, Matrix, m);
  Data_Get_Struct(expected, Matrix, n);

	if(m->nrows != n->nrows || m->ncols != n->ncols)
	{
		return Qfalse;
	}
  for(r = 0; r < m->nrows; r++)
  {
    for(c = 0; c < m->ncols; c++)
    { //use ruby to do this, as it will happily handle the complex data types.
      if( TYPE(rb_funcall(rb_blas_member_to_value(m, r, c),  rb_intern("within_bound"), 2, rb_blas_member_to_value(n, r, c), bound)) == TYPE(Qfalse) ) 
				return Qfalse;
    } 
  }
  return Qtrue;
}
