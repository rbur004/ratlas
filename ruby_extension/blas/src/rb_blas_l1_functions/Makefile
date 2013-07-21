include ../../../conf.mk

LIBDIR=../../lib

LIBOFILES = ${LIBDIR}/blas_l1_c2rb.o ${LIBDIR}/rb_atlas_xset_mod.o ${LIBDIR}/rb_blas_ixamax.o ${LIBDIR}/rb_blas_xasum.o ${LIBDIR}/rb_blas_xaxpy.o ${LIBDIR}/rb_blas_xaxpy_add.o ${LIBDIR}/rb_blas_xaxpy_mod.o ${LIBDIR}/rb_blas_xaxpy_sub.o ${LIBDIR}/rb_blas_xcopy_mod.o ${LIBDIR}/rb_blas_xdot.o ${LIBDIR}/rb_blas_xdotc.o ${LIBDIR}/rb_blas_xdotsds.o ${LIBDIR}/rb_blas_xdotu.o ${LIBDIR}/rb_blas_xnrm2.o ${LIBDIR}/rb_blas_xrot_mod.o ${LIBDIR}/rb_blas_xrotg.o ${LIBDIR}/rb_blas_xrotm_mod.o ${LIBDIR}/rb_blas_xrotmg.o ${LIBDIR}/rb_blas_xscal.o ${LIBDIR}/rb_blas_xscal_mod.o ${LIBDIR}/rb_blas_xswap_mod.o 
INCLUDEDIR=../../include
INCLUDEFILES= $(INCLUDEDIR)/rb_blas.h $(INCLUDEDIR)/rb_blas_extras.h $(INCLUDEDIR)/rb_blas_l1_c2rb.h $(INCLUDEDIR)/rb_blas_l2_c2rb.h $(INCLUDEDIR)/rb_blas_l3_c2rb.h $(INCLUDEDIR)/rb_blas_utils.h

all: ${LIBOFILES} $(INCLUDEFILES)

clean:
	rm -f ${LIBOFILES}

$(LIBDIR)/%.o: %.c
	$(CC)  $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $< 
