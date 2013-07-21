require 'ratlas'
include RAtlas
require 'complex' 
include Math
require '../testblas.rb'

class TestSwap < TestBlas

	def initialize 
	  super
    gsl
    ibm_examples
  end
  
  def test_swap(x,incx,y,incy,x_expected,y_expected, error_bound, test_message)
    x.swap!(y, incx,incy);
    print_on_error( "swap #{test_message}", x, x_expected, error_bound) 
    print_on_error( "swap #{test_message}", y, y_expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    test_swap(SingleBlas[ 0.539 ], 1, SingleBlas[ -0.262 ], -1, SingleBlas[ -0.262 ], SingleBlas[ 0.539 ], @flteps, "gsl Test 88/89")
    test_swap(DoubleBlas[ 0.906 ], 1, DoubleBlas[ 0.373 ], -1,  DoubleBlas[ 0.373 ], DoubleBlas[ 0.906 ],  @dbleps, "gsl Test 90/91")
#    test_swap(ComplexBlas[ Complex(-0.316, -0.529) ], 1, ComplexBlas[ Complex(-0.313, 0.363) ], -1, ComplexBlas[ Complex(-0.313, 0.363) ], ComplexBlas[Complex(-0.316, -0.529)],  Complex(@flteps,@flteps), "gsl Test 92/93")
#    test_swap(DoubleComplexBlas[ Complex(0.512, -0.89) ], 1, DoubleComplexBlas[ Complex(-0.225, -0.511) ], -1, DoubleComplexBlas[Complex(-0.225, -0.511)],   DoubleComplexBlas[Complex(0.512, -0.89) ], Complex(@dbleps,@dbleps), "gsl Test 94/95")

    test_swap(SingleBlas[ 0.336 ], -1, SingleBlas[ -0.431 ], 1, SingleBlas[ -0.431 ], SingleBlas[ 0.336 ], @flteps, "gsl Test 96/97")
    test_swap(DoubleBlas[ 0.764 ], -1, DoubleBlas[ -0.293 ], 1,  DoubleBlas[ -0.293 ], DoubleBlas[0.764 ],  @dbleps, "gsl Test 98/99")
#    test_swap(ComplexBlas[ Complex(-0.239, 0.361) ], -1, ComplexBlas[ Complex(0.149, 0.347) ], 1, ComplexBlas[ Complex(0.149, 0.347) ], ComplexBlas[Complex(-0.239, 0.361)],  Complex(@flteps,@flteps), "gsl Test 100/101")
#    test_swap(DoubleComplexBlas[ Complex(-0.171, -0.936) ], -1, DoubleComplexBlas[ Complex( 0.495, -0.835) ], 1, DoubleComplexBlas[Complex( 0.495, -0.835)],   DoubleComplexBlas[Complex(-0.171, -0.936) ], Complex(@dbleps,@dbleps), "gsl Test 102/103")
 
    test_swap(SingleBlas[ -0.405 ], -1, SingleBlas[ -0.213 ], -1, SingleBlas[ -0.213 ], SingleBlas[ -0.405 ], @flteps, "gsl Test 104/105")
    test_swap(DoubleBlas[ -0.761 ], -1, DoubleBlas[ -0.585 ], -1,  DoubleBlas[ -0.585 ], DoubleBlas[-0.761 ],  @dbleps, "gsl Test 106/107")
 #   test_swap(ComplexBlas[ Complex(0.853, 0.146) ], -1, ComplexBlas[ Complex( 0.009, -0.178) ], -1, ComplexBlas[ Complex( 0.009, -0.178) ], Complex(0.853, 0.146),  Complex(@flteps,@flteps), "gsl Test 108/109")
 #   test_swap(DoubleComplexBlas[ Complex(-0.228, 0.386) ], -1, DoubleComplexBlas[ Complex(  0.988, -0.084) ], -1, DoubleComplexBlas[Complex( 0.988, -0.084)],   DoubleComplexBlas[Complex(-0.228, 0.386) ], Complex(@dbleps,@dbleps), "gsl Test 110/111")

    test_swap(SingleBlas[ 0.539, 0.906 ], 1, SingleBlas[ -0.262, 0.373 ], 1, SingleBlas[ -0.262, 0.373 ], SingleBlas[0.539, 0.906 ], @flteps, "gsl Test 88/89")
    test_swap(DoubleBlas[ 0.539, 0.906], 1, DoubleBlas[ -0.262, 0.373 ], 1,  DoubleBlas[ -0.262, 0.373 ], DoubleBlas[0.539, 0.906],  @dbleps, "gsl Test 90/91")
  end
  
  def ibm_examples
        #Example 1
        #    This example shows vectors x and y with positive strides.
        #    Call Statement and Input:
        #
        #                N   X  INCX  Y  INCY
        #                |   |   |    |   |
        #    CALL SSWAP( 5 , X , 1  , Y , 2  )
        #
        #    X        =  (1.0, 2.0, 3.0, 4.0, 5.0)
        #    Y        =  (-1.0, . , -2.0, . , -3.0, . , -4.0, . , -5.0)
        #
        #    Output:
        #
        #    X        =  (-1.0, -2.0, -3.0, -4.0, -5.0)
        #    Y        =  (1.0, . , 2.0, . , 3.0, . , 4.0, . , 5.0)
        #
        test_swap(SingleBlas[ 1.0, 2.0, 3.0, 4.0, 5.0 ], 1, 
                  SingleBlas[ -1.0, 1.1 , -2.0, 2.1 , -3.0, 3.1 , -4.0, 4.1 , -5.0 ], 2,
                  SingleBlas[ -1.0, -2.0, -3.0, -4.0, -5.0 ], 
                  SingleBlas[1.0, 1.1 , 2.0, 2.1 , 3.0, 3.1 , 4.0, 4.1 , 5.0 ], 
                  @flteps, "ibm examples swap-01")
        #Example 2
        #    This example shows how to obtain output vectors x and y that are reverse copies of the input vectors y and x. You must specify strides with the same absolute value, but with opposite signs. For y, which has negative stride, processing begins at element Y(5), which is -5.0, and the results of the swap are stored beginning at the same element.
        #    Call Statement and Input:
        #
        #                N   X  INCX  Y   INCY
        #                |   |   |    |    |
        #    CALL SSWAP( 5 , X , 1  , Y , -1  )
        #
        #    X        =  (1.0, 2.0, 3.0, 4.0, 5.0)
        #    Y        =  (-1.0, -2.0, -3.0, -4.0, -5.0)
        #
        #    Output:
        #
        #    X        =  (-5.0, -4.0, -3.0, -2.0, -1.0)
        #    Y        =  (5.0, 4.0, 3.0, 2.0, 1.0)
        #
        test_swap(SingleBlas[ 1.0, 2.0, 3.0, 4.0, 5.0 ], 1, 
                  SingleBlas[ -1.0, -2.0, -3.0, -4.0, -5.0 ], -1,
                  SingleBlas[ -5.0, -4.0, -3.0, -2.0, -1.0 ], 
                  SingleBlas[5.0, 4.0, 3.0, 2.0, 1.0], 
                  @flteps, "ibm examples swap-02")
        #Example 3
        #    This example shows how SSWAP can be used to interchange scalar values in vectors x and y by specifying 0 strides and the number of elements to be processed as 1.
        #    Call Statement and Input:
        #
        #                N   X  INCX  Y  INCY
        #                |   |   |    |   |
        #    CALL SSWAP( 1 , X , 0  , Y , 0  )
        #
        #    X        =  (1.0)
        #    Y        =  (-4.0)
        #
        #    Output
        #
        #    X        =  (-4.0)
        #    Y        =  (1.0)
        #
        test_swap(SingleBlas[ 1.0 ], 0, 
                  SingleBlas[ -4.0 ], 0,
                  SingleBlas[ -4.0 ], 
                  SingleBlas[ 1.0], 
                  @flteps, "ibm examples swap-03")
        #Example 4
        #    This example shows vectors x and y, containing complex numbers and having positive strides.
        #    Call Statement and Input:
        #
        #                N   X  INCX  Y  INCY
        #                |   |   |    |   |
        #    CALL CSWAP( 4 , X , 1  , Y , 2  )
        #
        #    X        =  ((1.0, 6.0), (2.0, 7.0), (3.0, 8.0), (4.0, 9.0))
        #    Y        =  ((-1.0, -1.0), . , (-2.0, -2.0), . , (-3.0, -3.0), . ,
        #                 (-4.0, -4.0))
        #
        #    Output:
        #
        #    X        =  ((-1.0, -1.0), (-2.0, -2.0), (-3.0, -3.0), (-4.0, -4.0))
        #    Y        =  ((1.0, 6.0), . , (2.0, 7.0), . , (3.0, 8.0), . ,
        #                 (4.0, 9.0))
=begin
        test_swap(ComplexBlas[ Complex(1.0, 6.0), Complex(2.0, 7.0), Complex(3.0, 8.0), Complex(4.0, 9.0) ], 1, 
                  ComplexBlas[ Complex(-1.0, -1.0), Complex(-1.1, -1.1) , Complex(-2.0, -2.0), Complex(-2.1, -2.1) , Complex(-3.0, -3.0), Complex(-3.1, -3.1) , Complex(-4.0, -4.0)], 2,
                  ComplexBlas[ Complex(-1.0, -1.0), Complex(-2.0, -2.0), Complex(-3.0, -3.0), Complex(-4.0, -4.0)], 
                  ComplexBlas[ Complex(1.0, 6.0), Complex(-1.1, -1.1), Complex(2.0, 7.0), Complex(-2.1, -2.1), Complex(3.0, 8.0), Complex(-3.1, -3.1), Complex(4.0, 9.0)], 
                  @flteps, "ibm examples swap-04", 4)
=end
  end

end

TestSwap.new
