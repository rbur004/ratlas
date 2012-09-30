require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestCopy < TestBlas

	def initialize 
	  super
    gsl
  end
  
  def test_copy(x,incx,y,incy,x_expected,y_expected, error_bound, test_message)
    x.copy!(y, incx,incy);
    print_on_error( "copy #{test_message}", x, x_expected, error_bound) 
    print_on_error( "copy #{test_message}", y, y_expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    test_copy(SingleBlas[ 0.898 ], 1, SingleBlas[ 0.699 ], -1, SingleBlas[ 0.898 ], SingleBlas[ 0.898 ], @flteps, "gsl Test 76")
    test_copy(DoubleBlas[ 0.002 ], 1, DoubleBlas[ -0.921 ], -1,  DoubleBlas[ 0.002 ], DoubleBlas[ 0.002 ],  @dbleps, "gsl Test 77")
    test_copy(ComplexBlas[ Complex(-0.166, 0.639) ], 1, ComplexBlas[ Complex(0.863, 0.613) ], -1, ComplexBlas[ Complex(-0.166, 0.639) ], ComplexBlas[Complex(-0.166, 0.639)],  Complex(@flteps,@flteps), "gsl Test 78")
    test_copy(DoubleComplexBlas[ Complex(0.315, -0.324 ) ], 1, DoubleComplexBlas[ Complex(-0.312, -0.748) ], -1, DoubleComplexBlas[Complex(0.315, -0.324 )],   DoubleComplexBlas[Complex(0.315, -0.324 ) ], Complex(@dbleps,@dbleps), "gsl Test 79")

    test_copy(SingleBlas[  0.222 ], -1, SingleBlas[ 0.522 ], 1, SingleBlas[  0.222 ], SingleBlas[  0.222 ], @flteps, "gsl Test 80")
    test_copy(DoubleBlas[ 0.021 ], -1, DoubleBlas[  0.898 ], 1,  DoubleBlas[ 0.021 ], DoubleBlas[ 0.021 ],  @dbleps, "gsl Test 81")
    test_copy(ComplexBlas[ Complex(0.376, 0.229) ], -1, ComplexBlas[ Complex(0.143, -0.955) ], 1, ComplexBlas[ Complex(0.376, 0.229) ], ComplexBlas[Complex(0.376, 0.229)],  Complex(@flteps,@flteps), "gsl Test 82")
    test_copy(DoubleComplexBlas[ Complex(-0.265, -0.84) ], -1, DoubleComplexBlas[ Complex( -0.156, 0.939) ], 1, DoubleComplexBlas[Complex(-0.265, -0.84)],   DoubleComplexBlas[Complex(-0.265, -0.84) ], Complex(@dbleps,@dbleps), "gsl Test 83")
 
    test_copy(SingleBlas[ 0.074 ], -1, SingleBlas[ -0.802 ], -1, SingleBlas[ 0.074 ], SingleBlas[ 0.074 ], @flteps, "gsl Test 84")
    test_copy(DoubleBlas[ -0.374 ], -1, DoubleBlas[ -0.161  ], -1,  DoubleBlas[ -0.374 ], DoubleBlas[-0.374 ],  @dbleps, "gsl Test 85")
    test_copy(ComplexBlas[ Complex(0.084, 0.778) ], -1, ComplexBlas[ Complex( 0.31, -0.797) ], -1, ComplexBlas[ Complex( 0.084, 0.778) ], Complex(0.084, 0.778),  Complex(@flteps,@flteps), "gsl Test 86")
    test_copy(DoubleComplexBlas[ Complex(0.831, -0.282) ], -1, DoubleComplexBlas[ Complex(  -0.62, 0.32 ) ], -1, DoubleComplexBlas[Complex( 0.831, -0.282)],   DoubleComplexBlas[Complex(0.831, -0.282) ], Complex(@dbleps,@dbleps), "gsl Test 87")
  end

end

TestCopy.new
