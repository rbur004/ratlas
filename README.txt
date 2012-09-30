Ratlas is a Ruby library interfacing to the Atlas Blas and Lapack linear algebra c libraries. The Atlas libraries must be installed before this code will work. There are also variants of the clapack libraries with different calling and error handling. This code was written against the source forge release.

See the Atlas FAQ for details. http://math-atlas.sourceforge.net/faq.html

Atlas is released under a modified BSD license. To be compatible, this code is released under the same terms and condition as Atlas. 

Conventions:

Blas and Lapack are Fortran libraries. Atlas provides a C calling layer, which also allows both row and column ordered matrices. The library function's naming convention uses the first letter to indicate the data type the call uses.
	i.e. S - Real (or float), D - Double, C - Complex (float), Z - Double Complex.
	     SDS - Real Source and Result, but use Double's in calculations.
	
The next two letters indicate the Matrix type.
	i.e.	GE-General,		GB - General Banded
	   		HE-Hermitian,	HB - Hermitian Banded,	HP - Hermitian Packed
				SY-Symmetric,	SB - Symmetric Banded,	SP - Symmetric Packed
				TR-Triangular,TB - Triangular Banded,	TP - Triangular Packed
				
In this Ruby implementation, these two prefixes are not included in the name. The C library function called is determined by the Ruby object. e.g. The level one cblas library has functions cblas_saxpy(...), cblas_daxpy(...), cblas_caxpy(...) and cblas_zaxpy(...). The Ruby library only has a single method call, axpy(...). 

The cblas library also requires a number of arguments to be passed to calls to specify the matrix dimensions and layout. The Ruby calls reorder the arguments so that the common default values can be determined from inspection  of the matrix objects and only need to be included in the argument list if they are being overridden. 

The Ruby matrix object is stored in a C data structure exposed in Ruby as a family of Blas objects. The elements are only translated into Ruby numeric objects if they are directly accessed from the Ruby code. This includes a field to store the last scalar result so it may be used in the next call without having to go through a translation process. 

Complex numbers are stored as pairs of floats, and double complex numbers as pairs of doubles immediately adjacent in memory. This is the format used in cblas. When accessed in Ruby, both types become Ruby Complex objects.

There is no internal difference between a Vector and a Matrix. The Vector is stored as a single row Matrix.

There is a parent class Blas, which children classes FloatBlas, DoubleBlas, ComplexBlas, DoubleComplexBlas, GBFloatBlas, GBDoubleBlas, GBComplexBlas, GBDoubleComplexBlas, HEFloatBlas, HBDoubleBlas, ...

The Ruby interfaces to the cLapack library build on these Blas objects. Hence there are FloatLapack, DoubleLapack, ComplexLapack, ... objects which inherit from the Blas objects of the same date type and matrix storage type.


Building something that works:

First you need the Atlas cblas and clapack libraries installed, which you can get from source forge or they may well be installed already.

Next you need to run extconf.rb in the ruby_extension folder. It generates the Makefile. You might need to modify it to help it find the Atlas libraries (its only a few lines long.).

Then you type make and make install.

Testing:
The tests folder contains a large number of ruby based tests that have been copied from the gnu gsl projects cblas code. 



 			
