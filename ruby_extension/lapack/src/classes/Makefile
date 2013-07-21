include ../../../conf.mk

LIBDIR=../../lib
LIBOFILES = ${LIBDIR}/double_lapack.o ${LIBDIR}/rb_lapack.o ${LIBDIR}/single_lapack.o
INCLUDEDIR=../../include
INCLUDEFILES= $(INCLUDEDIR)/rb_lapack.h $(INCLUDEDIR)/rb_lapack_c2rb.h 

all: ${LIBOFILES} $(INCLUDEFILES)

clean:
	rm -f ${LIBOFILES}

$(LIBDIR)/%.o: %.c
	$(CC)  $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $< 
