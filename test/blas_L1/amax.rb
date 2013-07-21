require 'ratlas'
include RAtlas
require_relative '../testblas.rb'
require 'complex' 
include Math

class TestAsum < TestBlas

	def initialize 
	  super
    gsl
    ibm_examples
  end
  
  def test_amax(x,incx, expected, error_bound, test_message)
    result = x.amax(incx);
    #puts result
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
  
  def ibm_examples
    #From IBM Blas documentations examples. 
    # Example 1
    #     This example shows a vector, x, with a stride of 1.
    #     Function Reference and Input:
    #
    #                    N   X   INCX
    #                    |   |    |
    #     IMAX = ISAMAX( 9 , X ,  1   )
    #
    #     X        =  (1.0, 2.0, 7.0, -8.0, -5.0, -10.0, -9.0, 10.0, 6.0)
    #
    #     Output:
    #
    #     IMAX     =  6
    #
    test_amax(SingleBlas[ 1.0, 2.0, 7.0, -8.0, -5.0, -10.0, -9.0, 10.0, 6.0 ], 1, 5, 0, "ibm Test amax-01")
    # Example 2
    #     This example shows a vector, x, with a stride greater than 1.
    #     Function Reference and Input:
    #
    #                    N   X   INCX
    #                    |   |    |
    #     IMAX = ISAMAX( 5 , X ,  2   )
    #
    #     X        =  (1.0, . , 7.0, . , -5.0, . , -9.0, . , 6.0)
    #
    #     Output:
    #
    #     IMAX     =  4
    #
    test_amax(SingleBlas[ 1.0, 0.0 , 7.0, 0.0 , -5.0, 1000 , -9.0, 0.0 , 6.0 ], 2, 3, 0, "ibm Test amax-02")
    # Example 3
    #     This example shows a vector, x, with a stride of 0.
    #     Function Reference and Input:
    #
    #                    N   X   INCX
    #                    |   |    |
    #     IMAX = ISAMAX( 9 , X ,  0   )
    #
    #     X        =  (1.0, . , . , . , . , . , . , . , .)
    #
    #     Output:
    #
    #     IMAX     =  1
    #
    test_amax(SingleBlas[ 1.0, 0.0 , 7.0, 0.0 , -5.0, 1000 , -9.0, 0.0 , 6.0 ], 0, 0, 0, "ibm Test amax-03")
    # Example 4
    #     This example shows a vector, x, with a negative stride. Processing 
    #     begins at element X(15), which is 2.0.
    #     Function Reference and Input:
    #
    #                    N   X   INCX
    #                    |   |    |
    #     IMAX = ISAMAX( 8 , X , -2   )
    #
    #     X        =  (3.0, . , 5.0, . , -8.0, . , 6.0, . , 8.0, . ,
    #                  4.0, . , 8.0, . , 2.0)
    #
    #     Output:
    #
    #     IMAX     =  7
    #
    test_amax(SingleBlas[ 3.0, 0.0 , 5.0, 0.0 , -8.0, 0.0 , 6.0, 0.0 , 8.0, 0.0 ,
                 4.0, 0.0 , 8.0, 0.0 , 2.0 ], -2, 2, 0, "ibm Test amax-04")
     # Example 5
     #     This example shows a vector, x, containing complex numbers and having a 
     #     stride of 1.
     #     Function Reference and Input:
     #
     #                    N   X   INCX
     #                    |   |    |
     #     IMAX = ICAMAX( 5 , X ,  1   )
     #
     #     X        =  ((9.0 , 2.0) , (7.0 , -8.0) , (-5.0 , -10.0) , (-4.0 , 10.0),
     #                  (6.0 , 3.0))
     #
     #     Output:
     #
     #     IMAX     =  2
     test_amax(ComplexBlas[ Complex(9.0 , 2.0) , Complex(7.0 , -8.0) , Complex(-5.0 , -10.0) , Complex(-4.0 , 10.0), Complex(6.0 , 3.0)], -2, 1, 0, "ibm Test amax-05")

  end

end

TestAsum.new
