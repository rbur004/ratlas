require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestNrm2 < TestBlas

	def initialize
	  super
    gsl
  end
  
  def test_nrm2(x,incx, expected, error_bound, test_message)
    result = x.nrm2(incx);
    print_on_error( "nrm2 #{test_message}", result, expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    #gsl cblas returns 0 when incx is negative. Hence the gsl tests need modifying for the Atlas results!
    test_nrm2(SingleBlas[ 0.317 ], -1, 0.317, @flteps, "gsl Test 28")
    test_nrm2(DoubleBlas[ 0.071 ], -1, 0.071,  @dbleps, "gsl Test 29")
    test_nrm2(ComplexBlas[ Complex(0.776, 0.983) ], -1, 1.25238370895386,  @flteps, "gsl Test 30")
    test_nrm2(DoubleComplexBlas[ Complex(0.549, -0.354) ], -1, 0.653235792038373, @dbleps, "gsl Test 31")

    #2D Vector, incx 1
    test_nrm2(SingleBlas[ 0.14, -0.632 ], 1, 0.647320631527, @flteps, "gsl Test 32")
    test_nrm2(DoubleBlas[ 0.696, -0.804 ], 1, 1.06340584915,  @dbleps, "gsl Test 33")
    test_nrm2(ComplexBlas[ Complex(0.281, -0.063), Complex(0.367, 0.232) ], 1,  0.521001919382,  @flteps, "gsl Test 34")
    test_nrm2(DoubleComplexBlas[ Complex( -0.359, -0.76) , Complex(-0.906, -0.108)], 1, 1.24055672986, @dbleps, "gsl Test 35")

    #2D Vector, incx -1
    #gsl cblas returns 0 when incx is negative. Hence the gsl tests need modifying for the Atlas results!
    test_nrm2(SingleBlas[ 0.918, -0.126 ], -1, 0.926606714725494, @flteps, "gsl Test 36")
    test_nrm2(DoubleBlas[ 0.217, -0.588 ], -1, 0.626763910894685,  @dbleps, "gsl Test 37")
    test_nrm2(ComplexBlas[ Complex(0.31, 0.059), Complex(-0.442, 0.987) ], -1,  1.12654960155487,  @flteps, "gsl Test 38")
    test_nrm2(DoubleComplexBlas[ Complex( 0.609, 0.615) , Complex(-0.143, -0.957)], -1, 1.29823110423376, @dbleps, "gsl Test 39")

  end

end

TestNrm2.new
