include ../../../conf.mk

LIBDIR=../../lib
LIBOFILES = ${LIBDIR}/rb_crotg_class.o ${LIBDIR}/rb_drotg_class.o ${LIBDIR}/rb_drotmg_class.o ${LIBDIR}/rb_rotg.o ${LIBDIR}/rb_srotg_class.o ${LIBDIR}/rb_srotmg_class.o ${LIBDIR}/rb_zrotg_class.o 
INCLUDEDIR=../../include
INCLUDEFILES= $(INCLUDEDIR)/rb_rotg.h 

all: ${LIBOFILES} $(INCLUDEFILES)

clean:
	rm -f ${LIBOFILES}

$(LIBDIR)/%.o: %.c
	$(CC)  $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) -o $@ -c $< 
