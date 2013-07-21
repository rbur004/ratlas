require 'ratlas'
include RAtlas
require 'complex' 
include Math
require_relative '../testblas.rb'

class TestAsum < TestBlas

	def initialize 
	  super
    gsl
    ibm_examples
  end
  
  def test_asum(x,incx, n, expected, error_bound, test_message)
    result = n == nil ? x.asum(incx) : x.asum(incx, n)
    #puts result
    print_on_error( "asum #{test_message}", result, expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    #Atlas cblas treats -incx as walking through the vector in reverse.
    #gsl cblas returns 0 when incx is negative. Hence the gsl tests need modifying for the Atlas results!
    
    test_asum(SingleBlas[ 0.239 ], -1, nil, 0.239, @flteps, "gsl Test 40") 
    test_asum(DoubleBlas[  -0.413 ], -1, nil, 0.413,  @dbleps, "gsl Test 41") 
    test_asum(ComplexBlas[ Complex(0.1, 0.017) ], -1, nil, 0.117,  @flteps, "gsl Test 42")
    test_asum(DoubleComplexBlas[ Complex(-0.651, 0.079) ], -1, nil, 0.73, @dbleps, "gsl Test 43")

    #2D Vector, incx 1
    test_asum(SingleBlas[ 0.899, -0.72 ], 1, nil, 1.619, @flteps, "gsl Test 44")
    test_asum(DoubleBlas[ 0.271, -0.012 ], 1, nil, 0.283,  @dbleps, "gsl Test 45")
    test_asum(ComplexBlas[ Complex(-0.567, -0.645), Complex(0.098, 0.256) ], 1, nil,  1.566,  @flteps, "gsl Test 46")
    test_asum(DoubleComplexBlas[ Complex( -0.046, -0.671) , Complex(-0.323, 0.785)], 1, nil, 1.825, @dbleps, "gsl Test 47")

    #2D Vector, incx -1
    #Again Atlas cblas returns values, where gsl cblas returns 0
    test_asum(SingleBlas[ 0.169, 0.833 ], -1, nil, 1.002, @flteps, "gsl Test 48")
    test_asum(DoubleBlas[ -0.586, -0.486 ], -1, nil, 1.072,  @dbleps, "gsl Test 49")
    test_asum(ComplexBlas[ Complex(-0.314, -0.318), Complex(-0.835, -0.807) ], -1, nil,  2.274,  @flteps, "gsl Test 50")
    test_asum(DoubleComplexBlas[ Complex( -0.927, 0.152) , Complex(-0.554, -0.844)], -1, nil, 2.477, @dbleps, "gsl Test 51")

  end
  
  def ibm_examples
    #    Example 1
    #  This example shows a vector, x, with a stride of 1.
    #  Function Reference and Input:
    #
    #                N   X  INCX
    #                |   |   |
    #  SUMM = SASUM( 7 , X , 1  )
    #
    #  X        =  (1.0, -3.0, -6.0, 7.0, 5.0, 2.0, -4.0)
    #
    #  Output:
    #
    #  SUMM     =  28.0
    #
    test_asum(SingleBlas[ 1.0, -3.0, -6.0, 7.0, 5.0, 2.0, -4.0 ], 1, 7, 28.0, @flteps, "ibm example asum-01") 
    #    Example 2
    #  This example shows a vector, x, with a stride greater than 1.
    #  Function Reference and Input:
    #
    #                N   X  INCX
    #                |   |   |
    #  SUMM = SASUM( 4 , X , 2  )
    #
    #  X        =  (1.0, . , -6.0, . , 5.0, . , -4.0)
    #
    #  Output:
    #
    #  SUMM     =  16.0
    #
    test_asum(SingleBlas[ 1.0, 10.0 , -6.0, 10.0 , 5.0, 10.0 , -4.0 ], 2, 4, 16.0, @flteps, "ibm example asum-02") 
    #    Example 3
    #  This example shows a vector, x, with negative stride. Processing begins at element X(7), which is -4.0.
    #  Function Reference and Input:
    #
    #                N   X   INCX
    #                |   |    |
    #  SUMM = SASUM( 4 , X , -2  )
    #
    #  X        =  (1.0, . , -6.0, . , 5.0, . , -4.0)
    #
    #  Output:
    #
    #  SUMM     =  16.0
    #
    test_asum(SingleBlas[ 1.0, 10.0 , -6.0, 10.0 , 5.0, 10.0 , -4.0 ], -2, 4, 16.0, @flteps, "ibm example asum-03") 
    #    Example 4
    #  This example shows a vector, x, with a stride of 0. The result in SUMM is nx1.
    #  Function Reference and Input:
    #
    #                N   X  INCX
    #                |   |   |
    #  SUMM = SASUM( 7 , X , 0  )
    #
    #  X        =  (-2.0, . , . , . , . , . , .)
    #
    #  Output:
    #
    #  SUMM     =  14.0
    #
    #IBM example says the answer is 14, which is wrong!
    test_asum(SingleBlas[ -2.0, -3.0, -6.0, 7.0, 5.0, 2.0, -4.0 ], 0, 7, 0.0, @flteps, "ibm example asum-04") 
    #    Example 5
    #  This example shows a vector, x, containing complex numbers and having a stride of 1.
    #  Function Reference and Input:
    #
    #                 N   X  INCX
    #                 |   |   |
    #  SUMM = SCASUM( 5 , X , 1  )
    #
    #  X        =  ((1.0, 2.0), (-3.0, 4.0), (5.0, -6.0 ), (-7.0, -8.0),
    #               (9.0, 10.0))
    #
    #  Output:
    #
    #  SUMM     =  55.0
    test_asum(ComplexBlas[ Complex(1.0, 2.0), Complex(-3.0, 4.0), Complex(5.0, -6.0 ), Complex(-7.0, -8.0), Complex(9.0, 10.0) ], 1, 5, 55.0, @flteps, "ibm example asum-05") 

  
  end

end

TestAsum.new
