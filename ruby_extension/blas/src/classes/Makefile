include ../../../conf.mk

LIBDIR=../../lib
LIBOFILES = ${LIBDIR}/complex_blas.o ${LIBDIR}/double_blas.o ${LIBDIR}/double_complex_blas.o ${LIBDIR}/integer_blas.o ${LIBDIR}/rb_blas.o ${LIBDIR}/single_blas.o 
INCLUDEDIR=../../include
INCLUDEFILES= $(INCLUDEDIR)/rb_blas.h $(INCLUDEDIR)/rb_blas_extras.h $(INCLUDEDIR)/rb_blas_l1_c2rb.h $(INCLUDEDIR)/rb_blas_l2_c2rb.h $(INCLUDEDIR)/rb_blas_l3_c2rb.h $(INCLUDEDIR)/rb_blas_utils.h

all: ${LIBOFILES} $(INCLUDEFILES)

clean:
	rm -f ${LIBOFILES}

$(LIBDIR)/%.o: %.c
	$(CC)  $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $< 

#.c.o:
#	$(CC) $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<
