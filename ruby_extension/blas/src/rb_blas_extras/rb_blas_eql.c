#include "rb_blas.h"

VALUE rb_blas_eql(VALUE self, VALUE y)
{
  Matrix *m, *n;
  int r,c;

  Data_Get_Struct(self, Matrix, m);
  Data_Get_Struct(y, Matrix, n);
	if(m->nrows != n->nrows || m->ncols != n->ncols)
	{
		return Qfalse;
	}
		
  for(r = 0; r < m->nrows; r++)
  {
    for(c = 0; c < m->ncols; c++)
    { //use ruby to do this, as it will happily handle the complex data types.
      if (TYPE(rb_funcall(rb_blas_member_to_value(m, r, c),rb_intern("=="), 1, rb_blas_member_to_value(n, r, c))) == TYPE(Qfalse)) return Qfalse;
    } 
  }

  return Qtrue;
}
