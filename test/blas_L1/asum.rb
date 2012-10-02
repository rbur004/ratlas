require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestAsum < TestBlas

	def initialize 
	  super
    gsl
  end
  
  def test_asum(x,incx, expected, error_bound, test_message)
    result = x.asum(incx);
    print_on_error( "asum #{test_message}", result, expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    #Atlas cblas treats -incx as walking through the vector in reverse.
    #gsl cblas returns 0 when incx is negative. Hence the gsl tests need modifying for the Atlas results!
    
    test_asum(SingleBlas[ 0.239 ], -1, 0.239, @flteps, "gsl Test 40") 
    test_asum(DoubleBlas[  -0.413 ], -1, 0.413,  @dbleps, "gsl Test 41") 
    test_asum(ComplexBlas[ Complex(0.1, 0.017) ], -1, 0.117,  @flteps, "gsl Test 42")
    test_asum(DoubleComplexBlas[ Complex(-0.651, 0.079) ], -1, 0.73, @dbleps, "gsl Test 43")

    #2D Vector, incx 1
    test_asum(SingleBlas[ 0.899, -0.72 ], 1, 1.619, @flteps, "gsl Test 44")
    test_asum(DoubleBlas[ 0.271, -0.012 ], 1, 0.283,  @dbleps, "gsl Test 45")
    test_asum(ComplexBlas[ Complex(-0.567, -0.645), Complex(0.098, 0.256) ], 1,  1.566,  @flteps, "gsl Test 46")
    test_asum(DoubleComplexBlas[ Complex( -0.046, -0.671) , Complex(-0.323, 0.785)], 1, 1.825, @dbleps, "gsl Test 47")

    #2D Vector, incx -1
    #Again Atlas cblas returns values, where gsl cblas returns 0
    test_asum(SingleBlas[ 0.169, 0.833 ], -1, 1.002, @flteps, "gsl Test 48")
    test_asum(DoubleBlas[ -0.586, -0.486 ], -1, 1.072,  @dbleps, "gsl Test 49")
    test_asum(ComplexBlas[ Complex(-0.314, -0.318), Complex(-0.835, -0.807) ], -1,  2.274,  @flteps, "gsl Test 50")
    test_asum(DoubleComplexBlas[ Complex( -0.927, 0.152) , Complex(-0.554, -0.844)], -1, 2.477, @dbleps, "gsl Test 51")

  end

end

TestAsum.new
