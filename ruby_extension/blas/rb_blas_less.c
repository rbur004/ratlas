static VALUE rblas_less(VALUE obj, VALUE bound)
{
  //Definitely a hack. We want to be able to test the results with the same syntax as scalars
  //i.e. if (error = (result - expected)).abs > bound ...
  Matrix *m;
  int r,c;

  Data_Get_Struct(obj, Matrix, m);

  for(r = 0; r < m->nrows; r++)
  {
    for(c = 0; c < m->ncols; c++)
    { //use ruby to do this, as it will happily handle the complex data types.
      if (TYPE(rb_funcall(member_to_value(m, r, c),rb_intern("<"), 1, bound)) == TYPE(Qfalse)) return Qfalse;
    } 
  }
  return Qtrue;
}
