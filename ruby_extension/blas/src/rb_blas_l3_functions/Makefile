include ../../../conf.mk

LIBDIR=../../lib
LIBOFILES =  ${LIBDIR}/blas_l3_c2rb.o ${LIBDIR}/rb_blas_xgemm.o ${LIBDIR}/rb_blas_xhemm.o ${LIBDIR}/rb_blas_xher2k.o ${LIBDIR}/rb_blas_xherk.o ${LIBDIR}/rb_blas_xsymm.o ${LIBDIR}/rb_blas_xsyr2k.o ${LIBDIR}/rb_blas_xsyrk.o ${LIBDIR}/rb_blas_xtrmm.o ${LIBDIR}/rb_blas_xtrsm.o 
INCLUDEDIR=../../include
INCLUDEFILES= $(INCLUDEDIR)/rb_blas.h $(INCLUDEDIR)/rb_blas_extras.h $(INCLUDEDIR)/rb_blas_l1_c2rb.h $(INCLUDEDIR)/rb_blas_l2_c2rb.h $(INCLUDEDIR)/rb_blas_l3_c2rb.h $(INCLUDEDIR)/rb_blas_utils.h

all: ${LIBOFILES} $(INCLUDEFILES)

clean:
	rm -f ${LIBOFILES}

$(LIBDIR)/%.o: %.c
	$(CC)  $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $< 




