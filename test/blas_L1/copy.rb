require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestCopy < TestBlas

	def initialize 
	  super
	  ibm_examples
    gsl
  end
  
  def test_copy(x,incx,y,incy,x_expected,y_expected, error_bound, test_message, n = nil)
    #puts y
    if n == nil
      y.copy!(x, incy,incx);
    else
      y.copy!(x, incy,incx,n);
    end
    #puts x
    #puts y
    print_on_error( "copy! #{test_message}", x, x_expected, error_bound) 
    print_on_error( "copy! #{test_message}", y, y_expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    test_copy(SingleBlas[ 0.898 ], 1, SingleBlas[ 0.699 ], -1, SingleBlas[ 0.898 ], SingleBlas[ 0.898 ], @flteps, "gsl Test 76")
    test_copy(DoubleBlas[ 0.002 ], 1, DoubleBlas[ -0.921 ], -1,  DoubleBlas[ 0.002 ], DoubleBlas[ 0.002 ],  @dbleps, "gsl Test 77")
    #test_copy(ComplexBlas[ Complex(-0.166, 0.639) ], 1, ComplexBlas[ Complex(0.863, 0.613) ], -1, ComplexBlas[ Complex(-0.166, 0.639) ], ComplexBlas[Complex(-0.166, 0.639)],  Complex(@flteps,@flteps), "gsl Test 78")
    #test_copy(DoubleComplexBlas[ Complex(0.315, -0.324 ) ], 1, DoubleComplexBlas[ Complex(-0.312, -0.748) ], -1, DoubleComplexBlas[Complex(0.315, -0.324 )],   DoubleComplexBlas[Complex(0.315, -0.324 ) ], Complex(@dbleps,@dbleps), "gsl Test 79")

    test_copy(SingleBlas[  0.222 ], -1, SingleBlas[ 0.522 ], 1, SingleBlas[  0.222 ], SingleBlas[  0.222 ], @flteps, "gsl Test 80")
    test_copy(DoubleBlas[ 0.021 ], -1, DoubleBlas[  0.898 ], 1,  DoubleBlas[ 0.021 ], DoubleBlas[ 0.021 ],  @dbleps, "gsl Test 81")
    #test_copy(ComplexBlas[ Complex(0.376, 0.229) ], -1, ComplexBlas[ Complex(0.143, -0.955) ], 1, ComplexBlas[ Complex(0.376, 0.229) ], ComplexBlas[Complex(0.376, 0.229)],  Complex(@flteps,@flteps), "gsl Test 82")
    #test_copy(DoubleComplexBlas[ Complex(-0.265, -0.84) ], -1, DoubleComplexBlas[ Complex( -0.156, 0.939) ], 1, DoubleComplexBlas[Complex(-0.265, -0.84)],   DoubleComplexBlas[Complex(-0.265, -0.84) ], Complex(@dbleps,@dbleps), "gsl Test 83")
 
    test_copy(SingleBlas[ 0.074 ], -1, SingleBlas[ -0.802 ], -1, SingleBlas[ 0.074 ], SingleBlas[ 0.074 ], @flteps, "gsl Test 84")
    test_copy(DoubleBlas[ -0.374 ], -1, DoubleBlas[ -0.161  ], -1,  DoubleBlas[ -0.374 ], DoubleBlas[-0.374 ],  @dbleps, "gsl Test 85")
    #test_copy(ComplexBlas[ Complex(0.084, 0.778) ], -1, ComplexBlas[ Complex( 0.31, -0.797) ], -1, ComplexBlas[ Complex( 0.084, 0.778) ], Complex(0.084, 0.778),  Complex(@flteps,@flteps), "gsl Test 86")
    #test_copy(DoubleComplexBlas[ Complex(0.831, -0.282) ], -1, DoubleComplexBlas[ Complex(  -0.62, 0.32 ) ], -1, DoubleComplexBlas[Complex( 0.831, -0.282)],   DoubleComplexBlas[Complex(0.831, -0.282) ], Complex(@dbleps,@dbleps), "gsl Test 87")
  end

  def ibm_examples
    #    Example 1
    #        This example shows input vector x and output vector y with positive strides.
    #        Call Statement and Input:
    #
    #                    N   X  INCX  Y  INCY
    #                    |   |   |    |   |
    #        CALL SCOPY( 5 , X , 1  , Y , 2  )
    #
    #        X        =  (1.0, 2.0, 3.0, 4.0, 5.0)
    #
    #        Output:
    #
    #        Y        =  (1.0, . , 2.0, . , 3.0, . , 4.0, . , 5.0)
    #
    test_copy(SingleBlas[ 1.0, 2.0, 3.0, 4.0, 5.0 ], 1, 
              SingleBlas[ 29.0, 33.0 , 28.0, 33.0 , 27.0, 33.0 , 26.0, 33.0 , 25.0 ], 2, 
              SingleBlas[ 1.0, 2.0, 3.0, 4.0, 5.0 ], 
              SingleBlas[1.0, 33.0 , 2.0, 33.0 , 3.0, 33.0 , 4.0, 33.0 , 5.0], @flteps, "ibm example copy 01", 5)
    #    Example 2
    #        This example shows how to obtain a reverse copy of the input vector x by
    #        specifying strides with the same absolute value, but with opposite signs, 
    #        for input vector x and output vector y. For y, which has a negative stride,
    #        results are stored beginning at element Y(5).
    #        Call Statement and Input:
    #
    #                    N   X  INCX  Y   INCY
    #                    |   |   |    |    |
    #        CALL SCOPY( 5 , X , 1  , Y , -1  )
    #
    #        X        =  (1.0, 2.0, 3.0, 4.0, 5.0)
    #
    #        Output:
    #
    #        Y        =  (5.0, 4.0, 3.0, 2.0, 1.0)
    #
    test_copy(SingleBlas[ 1.0, 2.0, 3.0, 4.0, 5.0 ], 1, 
              SingleBlas[ 0.0, 0.0, 0.0, 0.0, 0,0 ], -1, 
              SingleBlas[ 1.0, 2.0, 3.0, 4.0, 5.0 ], 
              SingleBlas[5.0, 4.0, 3.0, 2.0, 1.0], @flteps, "ibm example copy 02", 5)
    #    Example 3
    #        This example shows an input vector, x, with 0 stride. Vector x is treated 
    #        like a vector of length n, all of whose elements are the same as the single
    #        element in x. This is a technique for replicating an element of a vector.
    #        Call Statement and Input:
    #
    #                    N   X  INCX  Y  INCY
    #                    |   |   |    |   |
    #        CALL SCOPY( 5 , X , 0  , Y , 1  )
    #
    #        X        =  (13.0)
    #
    #        Output:
    #
    #        Y        =  (13.0, 13.0, 13.0, 13.0, 13.0)
    #
    test_copy(SingleBlas[ 13.0 ], 0, 
              SingleBlas[ 0.0, 0.0, 0.0, 0.0, 0,0 ], 1, 
              SingleBlas[ 13.0 ], 
              SingleBlas[13.0, 13.0, 13.0, 13.0, 13.0], @flteps, "ibm example copy 03", 5)
    #    Example 4
    #        This example shows input vector x and output vector y, containing complex
    #        numbers and having positive strides.
    #        Call Statement and Input:
    #
    #                    N   X  INCX  Y  INCY
    #                    |   |   |    |   |
    #        CALL CCOPY( 4 , X , 1  , Y , 2  )
    #
    #        X        =  ((1.0, 1.0), (2.0, 2.0), (3.0, 3.0), (4.0, 4.0))
    #
    #        Output:
    #
    #        Y        =  ((1.0, 1.0), . , (2.0, 2.0), . , (3.0, 3.0), . ,
    #                     (4.0, 4.0))
=begin
    test_copy(ComplexBlas[ Complex(1.0, 1.0), Complex(2.0, 2.0), Complex(3.0, 3.0), Complex(4.0, 4.0) ], 1, 
              ComplexBlas[ Complex( 0.0, 0.0),Complex( 0.0, 0.0),Complex( 0.0, 0.0),Complex( 0.0, 0.0),Complex( 0.0, 0.0),Complex( 0.0, 0.0),Complex( 0.0, 0.0) ], 2,
              ComplexBlas[ Complex(1.0, 1.0), Complex(2.0, 2.0), Complex(3.0, 3.0), Complex(4.0, 4.0)  ], 
              ComplexBlas[ Complex(1.0, 1.0),Complex( 0.0, 0.0), Complex(2.0, 2.0),Complex( 0.0, 0.0), Complex(3.0, 3.0),Complex( 0.0, 0.0), Complex(4.0, 4.0) ],  
              Complex(@flteps,@flteps), "ibm example copy 04", 4)
=end
  end
  
end

TestCopy.new
