#ifndef RB_ROTG_H
#define RB_ROTG_H

typedef struct Rotg //dummy base class. Needed for no real reason that I can currently think of. 
{
  VALUE class_id;
} Rotg;

typedef struct sRotg //single precision version
{
  VALUE class_id;
  float r, z, c, s;
} sRotg;

typedef struct dRotg //double precision version
{
  VALUE class_id;
  double r, z, c, s;
} dRotg;

typedef struct cRotg //single precision imaginary version
{
  VALUE class_id;
  float r[2], z[2], c, s[2];
} cRotg;

typedef struct zRotg //double precision imaginary version
{
  VALUE class_id;
  double r[2], z[2], c, s[2];
} zRotg;

typedef struct sRotmg //single precision modified plane rotation
{
  VALUE class_id;
  float d1, d2, x1, y1, param[5]; 
} sRotmg;

typedef struct dRotmg //double precision modified plane rotation
{
  VALUE class_id;
  double d1, d2, x1, y1, param[5]; 
} dRotmg;

VALUE Init_Rotg(VALUE myModule);

void Init_srotg(VALUE myModule, VALUE parent_class) ;
void Init_drotg(VALUE myModule, VALUE parent_class) ;
void Init_crotg(VALUE myModule, VALUE parent_class) ;
void Init_zrotg(VALUE myModule, VALUE parent_class) ;
void Init_srotgm(VALUE myModule, VALUE parent_class);
void Init_drotgm(VALUE myModule, VALUE parent_class) ;

#endif //RB_ROTG_H