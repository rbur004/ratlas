include ../../../conf.mk

LIBDIR=../../lib
LIBOFILES = ${LIBDIR}/blas_l2_c2rb.o ${LIBDIR}/rb_blas_mv.o ${LIBDIR}/rb_blas_r.o ${LIBDIR}/rb_blas_r2.o ${LIBDIR}/rb_blas_rc.o ${LIBDIR}/rb_blas_sv.o 
INCLUDEDIR=../../include
INCLUDEFILES= $(INCLUDEDIR)/rb_blas.h $(INCLUDEDIR)/rb_blas_extras.h $(INCLUDEDIR)/rb_blas_l1_c2rb.h $(INCLUDEDIR)/rb_blas_l2_c2rb.h $(INCLUDEDIR)/rb_blas_l3_c2rb.h $(INCLUDEDIR)/rb_blas_utils.h

all: ${LIBOFILES} $(INCLUDEFILES)

clean:
	rm -f ${LIBOFILES}

$(LIBDIR)/%.o: %.c
	$(CC)  $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $< 




