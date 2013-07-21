include ../../../conf.mk

LIBDIR=../../lib
LIBOFILES = ${LIBDIR}/rb_blas_util.o ${LIBDIR}/rb_blas_abs.o ${LIBDIR}/rb_blas_eql.o ${LIBDIR}/rb_blas_within_bound.o ${LIBDIR}/rb_blas_times.o
INCLUDEDIR=../../include
INCLUDEFILES= $(INCLUDEDIR)/rb_blas.h $(INCLUDEDIR)/rb_blas_extras.h $(INCLUDEDIR)/rb_blas_l1_c2rb.h $(INCLUDEDIR)/rb_blas_l2_c2rb.h $(INCLUDEDIR)/rb_blas_l3_c2rb.h $(INCLUDEDIR)/rb_blas_utils.h

all: ${LIBOFILES} $(INCLUDEFILES)

clean:
	rm -f ${LIBOFILES}

$(LIBDIR)/%.o: %.c
	$(CC)  $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $< 




