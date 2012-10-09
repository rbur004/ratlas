static VALUE rblas_abs(VALUE obj)
{
  //Definitely a hack. We want to be able to test the results with the same syntax as scalars
  //i.e. if (error = (result - expected)).abs > bound ...
  Matrix *m, *abs_m;
  int r,c;
  VALUE abs_matrix;

  Data_Get_Struct(obj, Matrix, m);
  //argv[0] = rblas_new_instance(dy->class_id, dy->data, dy->nrows, dy->ncols, dy->data_type);
  abs_matrix = rblas_new_instance( m->class_id, NULL, m->nrows, m->ncols, m->data_type, m->matrix_type, m->cblas_order);
  Data_Get_Struct(abs_matrix, Matrix, abs_m);

  for(r = 0; r < m->nrows; r++)
  {
    for(c = 0; c < m->ncols; c++)
    { //use ruby to do this, as it will happily handle the complex data types.
      value_to_member(abs_m, r, c, rb_funcall(member_to_value(m, r, c),  rb_intern("abs"), 0) );
    } 
  }

  return abs_matrix;
}

