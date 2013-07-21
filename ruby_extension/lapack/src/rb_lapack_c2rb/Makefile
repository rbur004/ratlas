include ../../../conf.mk

LIBDIR=../../lib
LIBOFILES = ${LIBDIR}/rb_lapack_c2rb.o ${LIBDIR}/rb_lapack_determinant.o ${LIBDIR}/rb_lapack_xgesv.o ${LIBDIR}/rb_lapack_xgetrf.o ${LIBDIR}/rb_lapack_xgetrs.o 
INCLUDEDIR=../../include
INCLUDEFILES= $(INCLUDEDIR)/rb_lapack.h $(INCLUDEDIR)/rb_lapack_c2rb.h 

all: ${LIBOFILES} $(INCLUDEFILES)

clean:
	rm -f ${LIBOFILES}

$(LIBDIR)/%.o: %.c
	$(CC)  $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $< 
