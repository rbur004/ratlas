require 'ratlas'
include RAtlas
require '../testblas.rb'
require 'complex' 
include Math

class TestAsum < TestBlas

	def initialize 
	  super
    gsl
  end
  
  def test_amax(x,incx, expected, error_bound, test_message)
    result = x.amax(incx);
    print_on_error( "amax #{test_message}", result, expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    test_amax(SingleBlas[ -0.388 ], -1, 0, 0, "gsl Test 52")
    test_amax(DoubleBlas[  0.247 ], -1, 0,  0, "gsl Test 53")
    test_amax(ComplexBlas[ Complex(0.704, 0.665) ], -1, 0,  0, "gsl Test 54")
    test_amax(DoubleComplexBlas[ Complex(-0.599, -0.758 ) ], -1, 0, 0, "gsl Test 55")

    #2D Vector, incx 1
    test_amax(SingleBlas[ 0.909, 0.037 ], 1, 0, 0, "gsl Test 56")
    test_amax(DoubleBlas[ 0.271, -0.426 ], 1, 1,  0, "gsl Test 57")
    test_amax(ComplexBlas[ Complex(-0.648, 0.317), Complex(0.62, 0.392) ], 1,  1,  0, "gsl Test 58")
    test_amax(DoubleComplexBlas[ Complex( -0.789, 0.352 ) , Complex(0.562, 0.697)], 1, 1, 0, "gsl Test 59")

    #2D Vector, incx -1
    #gsl cblas returns 0 when incx is negative. Hence the gsl tests need modifying for the Atlas results!
    test_amax(SingleBlas[ 0.487, 0.918 ], -1, 1, 0, "gsl Test 60")
    test_amax(DoubleBlas[ 0.537, 0.826 ], -1, 1,  0, "gsl Test 61")
    test_amax(ComplexBlas[ Complex(0.993, 0.172), Complex(-0.825, 0.873) ], -1,  1,  0, "gsl Test 62")
    test_amax(DoubleComplexBlas[ Complex(  0.913, -0.436) , Complex( -0.134, 0.129)], -1, 0, 0, "gsl Test 63")

  end

end

TestAsum.new
