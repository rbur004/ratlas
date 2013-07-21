require 'ratlas'
include RAtlas
require 'complex' 
include Math
require_relative '../testblas.rb'

class TestNrm2 < TestBlas

	def initialize
	  super
   # gsl
    ibm_examples
  end
  
  def test_nrm2(x,incx, expected, error_bound, test_message, n = nil)
    result = n == nil ? x.nrm2(incx) : x.nrm2(incx,n);
    #puts "result = #{result} expected = #{expected} +/- #{error_bound}"
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
  
  def ibm_examples
        #Example 1
        #    This example shows a vector, x, whose elements must be 
        #    scaled to prevent overflow.
        #
        #                   N   X  INCX
        #                   |   |   |
        #    DNORM = DNRM2( 6 , X , 1  )
        #
        #    X      = (0.68056D+200, 0.25521D+200, 0.34028D+200,
        #              0.85071D+200, 0.25521D+200, 0.85071D+200)
        #
        #    Output:
        #
        #    DNORM    =  0.1469D+201
        #
        test_nrm2(DoubleBlas[ 0.68056E200, 0.25521E200, 0.34028E200, 0.85071E200, 0.25521E200, 0.85071E200 ], 1, 0.1469E201,  0.0001E201, "ibm examples nrm2-01")
        #Example 2
        #    This example shows a vector, x, whose elements must be scaled 
        #    to prevent destructive underflow.
        #    Function Reference and Input:
        #
        #                   N   X  INCX
        #                   |   |   |
        #    DNORM = DNRM2( 4 , X , 2  )
        #
        #    X     = (0.10795D-200, . , 0.10795D-200, . , 0.10795D-200,
        #             . , 0.10795D-200)
        #
        #    Output:
        #
        #    DNORM    =  0.21590D-200
        #
        test_nrm2(DoubleBlas[ 0.10795E-200, 0.0 , 0.10795E-200, 0.0 , 0.10795E-200, 0.0 , 0.10795E-200 ], 2, 0.21590E-200, 0.00001E-200, "ibm examples nrm2-02", 4)
        #Example 3
        #    This example shows a vector, x, with a stride of 0. The result in SNORM is:
        #    Math Graphic
        #    Function Reference and Input:
        #
        #                   N   X  INCX
        #                   |   |   |
        #    SNORM = SNRM2( 4 , X , 0  )
        #
        #    X        =  (4.0)
        #
        #    Output:
        #
        #    SNORM    =  8.0
        #
        puts "nrm2-03: I get 0 as the result with my blas library unless I set incx to 1?"
        test_nrm2(SingleBlas[ 4.0 ], 0, 8.0, @flteps, "ibm examples nrm2-03", 4)
        #Example 4
        #    This example shows a vector, x, containing complex numbers,
        #    and whose elements must be scaled to prevent overflow.
        #    Function Reference and Input:
        #
        #                     N   X  INCX
        #                     |   |   |
        #    DZNORM = DZNRM2( 3 , X , 1  )
        #
        #    X     = ((0.68056D+200, 0.25521D+200), (0.34028D+200, 0.85071D+200),
        #             (0.25521D+200, 0.85071D+200))
        #
        #    Output:
        #
        #    DZNORM   =  0.1469D+201
        #
        test_nrm2(DoubleComplexBlas[ Complex( 0.68056E200, 0.25521E200) , Complex(0.34028E200, 0.85071E200), Complex(0.25521E200, 0.85071E200)], 1, 0.1469E201, 0.0001E201, "ibm examples nrm2-04", 3)
        #Example 5
        #    This example shows a vector, x, containing complex numbers, 
        #    and whose elements must be scaled to prevent destructive underflow.
        #    Function Reference and Input:
        #
        #                     N   X  INCX
        #                     |   |   |
        #    DZNORM = DZNRM2( 2 , X , 2  )
        #
        #    X     = ((0.10795D-200, 0.10795D-200), . ,
        #             (0.10795D-200, 0.10795D-200))
        #
        #    Output:
        #
        #    DZNORM   =  0.2159D-200
        test_nrm2(DoubleComplexBlas[ Complex( 0.10795E-200, 0.10795E-200) , Complex(0.34028E200, 0.85071E200), Complex(0.10795E-200, 0.10795E-200)], 2, 0.2159E-200, 0.0001E-200, "ibm examples nrm2-05", 2)
    
  end

end

TestNrm2.new
