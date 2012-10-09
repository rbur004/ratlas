INCLUDE	=	-I. -I.. -I/usr/lib/ruby/1.8/i686-darwin9.3.0 -I/Users/rbur004/Documents/src/ratlas/trunk/ruby_extension -I/Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/vecLib.framework/Versions/A/Headers/ 
CFLAGS   =  -fno-common -g -O2 -pipe -fno-common $(INCLUDE) 


../librotg.a: rb_rotg.o rb_crotg_class.o rb_drotg_class.o rb_srotg_class.o rb_zrotg_class.o rb_srotmg_class.o rb_drotmg_class.o
	AR -r $@ $+