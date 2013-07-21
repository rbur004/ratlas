require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestDot < TestBlas

  def initialize
    super
    gsl
    ibm_examples
  end
  
  def test_dotsds(x,incx,y,incy,alpha , expected, error_bound, test_message)
    result = x.dotsds(y, alpha,  incx, incy);
    print_on_error( "dotsds #{test_message}", result, expected, error_bound) 
  end
  
  def test_dot(x,incx,y,incy , expected, error_bound, test_message, n = nil)
    result = n == nil ? x.dot( y,  incx, incy) : x.dot( y,  incx, incy, n);
    print_on_error( "dot #{test_message}", result, expected, error_bound) 
  end
  
  def test_dotc(x,incx,y,incy , expected, error_bound, test_message)
    result = x.dotc( y,  incx, incy);
    print_on_error( "dotc #{test_message}", result, expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite. 
=begin
    test_dotsds(SingleBlas[ 0.733 ], 1, SingleBlas[ 0.825 ], -1, 0.0,   0.604725,  @flteps, "gsl Test 1")
    test_dotsds(SingleBlas[ 0.733 ], 1, SingleBlas[ 0.825 ], -1, 0.1,   0.704725,  @flteps, "gsl Test 2")
    test_dotsds(SingleBlas[ 0.733 ], 1, SingleBlas[ 0.825 ], -1, 1.0,   1.604725,  @flteps, "gsl Test 3")
    
    test_dotsds(SingleBlas[ -0.812 ], -1, SingleBlas[ -0.667 ], 1, 0.0,   0.541604,  @flteps, "gsl Test 4")
    test_dotsds(SingleBlas[ -0.812 ], -1, SingleBlas[ -0.667 ], 1, 0.1,   0.641604,  @flteps, "gsl Test 5")
    test_dotsds(SingleBlas[ -0.812 ], -1, SingleBlas[ -0.667 ], 1, 1.0,   1.541604,  @flteps, "gsl Test 6")
    
    test_dotsds(SingleBlas[ 0.481 ], -1, SingleBlas[ 0.523 ], -1, 0.0,   0.251563,  @flteps, "gsl Test 7")
    test_dotsds(SingleBlas[ 0.481 ], -1, SingleBlas[ 0.523 ], -1, 0.1,   0.351563,  @flteps, "gsl Test 8")
    test_dotsds(SingleBlas[ 0.481 ], -1, SingleBlas[ 0.523 ], -1, 1.0,   1.251563,  @flteps, "gsl Test 9")
=end
    test_dot(SingleBlas[ 0.785 ], 1, SingleBlas[ -0.7 ], -1,  -0.5495,  @flteps, "gsl Test 10")
    test_dot(DoubleBlas[ 0.79 ], 1, DoubleBlas[ -0.679 ], -1, -0.53641,  @dbleps, "gsl Test 11")
#    test_dot(ComplexBlas[ Complex(0.474, -0.27) ], 1, ComplexBlas[ Complex(-0.144, -0.392) ], -1, Complex(-0.174096, -0.146928),  Complex(@flteps,@flteps), "gsl Test 12")
#    test_dotc(ComplexBlas[ Complex(0.474, -0.27) ], 1, ComplexBlas[ Complex(-0.144, -0.392) ], -1, Complex(0.037584, -0.224688),  Complex(@flteps,@flteps), "gsl Test 13")
#    test_dot(DoubleComplexBlas[ Complex(-0.87, -0.631) ], 1, DoubleComplexBlas[ Complex(-0.7, -0.224) ], -1, Complex(0.467656, 0.63658),  Complex(@dbleps,@dbleps), "gsl Test 14")
#    test_dotc(DoubleComplexBlas[ Complex(-0.87, -0.631) ], 1, DoubleComplexBlas[ Complex(-0.7, -0.224) ], -1, Complex(0.750344, -0.24682),  Complex(@dbleps,@dbleps), "gsl Test 15")

    test_dot(SingleBlas[ -0.457 ], -1, SingleBlas[ 0.839 ], 1, -0.383423,  @flteps, "gsl Test 16")
    test_dot(DoubleBlas[  0.949 ], -1, DoubleBlas[ -0.873 ], 1, -0.828477,  @dbleps, "gsl Test 17")
#    test_dot(ComplexBlas[ Complex(0.852, -0.045) ], -1, ComplexBlas[ Complex(0.626, -0.164) ], 1, Complex(0.525972, -0.167898), Complex(@flteps,@flteps), "gsl Test 18")
#    test_dotc(ComplexBlas[ Complex(0.852, -0.045) ], -1, ComplexBlas[ Complex(0.626, -0.164) ], 1, Complex(0.540732, -0.111558),  Complex(@flteps,@flteps), "gsl Test 19")
#    test_dot(DoubleComplexBlas[ Complex(-0.786, -0.341) ], -1, DoubleComplexBlas[ Complex(-0.271, -0.896) ], 1, Complex(-0.09253, 0.796667),  Complex(@dbleps,@dbleps), "gsl Test 20")
#    test_dotc(DoubleComplexBlas[ Complex(-0.786, -0.341) ], -1, DoubleComplexBlas[ Complex(-0.271, -0.896) ], 1, Complex(0.518542, 0.611845),  Complex(@dbleps,@dbleps), "gsl Test 21")

    test_dot(SingleBlas[ -0.088 ], -1, SingleBlas[ -0.165 ], -1, 0.01452,  @flteps, "gsl Test 22")
    test_dot(DoubleBlas[  -0.434 ], -1, DoubleBlas[ -0.402 ], -1, 0.174468,  @dbleps, "gsl Test 23")
#    test_dot(ComplexBlas[ Complex(-0.347, 0.899) ], -1, ComplexBlas[ Complex(-0.113, -0.858) ], -1, Complex(0.810553, 0.196139),  Complex(@flteps,@flteps), "gsl Test 24")
#    test_dotc(ComplexBlas[ Complex(-0.347, 0.899) ], -1, ComplexBlas[ Complex(-0.113, -0.858) ], -1, Complex(-0.732131, 0.399313),  Complex(@flteps,@flteps), "gsl Test 25")
#    test_dot(DoubleComplexBlas[ Complex(-0.897, -0.204 ) ], -1, DoubleComplexBlas[ Complex(-0.759, 0.557) ], -1, Complex(0.794451, -0.344793),  Complex(@dbleps,@dbleps), "gsl Test 26")
#    test_dotc(DoubleComplexBlas[ Complex(-0.897, -0.204 ) ], -1, DoubleComplexBlas[ Complex(-0.759, 0.557) ], -1, Complex(0.567195, -0.654465),  Complex(@dbleps,@dbleps), "gsl Test 27")
  end
  
  def ibm_examples
        #    Example 1
        #        This example shows how to compute the dot product of two vectors, 
        #        x and y, having strides of 1.
        #        Function Reference and Input:
        #
        #                     N   X  INCX  Y  INCY
        #                     |   |   |    |   |
        #        DOTT = SDOT( 5 , X , 1  , Y , 1  )
        #
        #        X        =  (1.0, 2.0, -3.0, 4.0, 5.0)
        #        Y        =  (9.0, 8.0, 7.0, -6.0, 5.0)
        #
        #        Output:
        #
        #        DOTT     =  (9.0 + 16.0 - 21.0 - 24.0 + 25.0) = 5.0
        #
        test_dot(SingleBlas[ 1.0, 2.0, -3.0, 4.0, 5.0 ], 1, 
                 SingleBlas[ 9.0, 8.0, 7.0, -6.0, 5.0 ], 1, 
                 5.0,  @flteps, "ibm example dot-01")
        #    Example 2
        #        This example shows how to compute the dot product of a vector, x, 
        #        with a stride of 1, and a vector, y, with a stride greater than 1.
        #        Function Reference and Input:
        #
        #                     N   X  INCX  Y  INCY
        #                     |   |   |    |   |
        #        DOTT = SDOT( 5 , X , 1  , Y , 2  )
        #
        #        X        =  (1.0, 2.0, -3.0, 4.0, 5.0)
        #        Y        =  (9.0, . , 7.0, . , 5.0, . , -3.0, . , 1.0)
        #
        #        Output:
        #
        #        DOTT     =  (9.0 + 14.0 - 15.0 - 12.0 + 5.0) = 1.0
        #
        test_dot(SingleBlas[ 1.0, 2.0, -3.0, 4.0, 5.0 ], 1, 
                 SingleBlas[ 9.0, 1.0 , 7.0, 2.0 , 5.0, 3.0 , -3.0, 4.0 , 1.0 ], 2, 
                 1.0,  @flteps, "ibm example dot-02", 5)
        #    Example 3
        #        This example shows how to compute the dot product of a vector, x, 
        #        with a negative stride, and a vector, y, with a stride greater than 1.
        #        For x, processing begins at element X(5), which is 5.0.
        #        Function Reference and Input:
        #
        #                     N   X   INCX  Y  INCY
        #                     |   |    |    |   |
        #        DOTT = SDOT( 5 , X , -1  , Y , 2  )
        #
        #        X        =  (1.0, 2.0, -3.0, 4.0, 5.0)
        #        Y        =  (9.0, . , 7.0, . , 5.0, . , -3.0, . , 1.0)
        #
        #        Output:
        #
        #        DOTT     =  (45.0 + 28.0 - 15.0 - 6.0 + 1.0) = 53.0
        #
        test_dot(SingleBlas[ 1.0, 2.0, -3.0, 4.0, 5.0 ], -1, 
                 SingleBlas[ 9.0, 1.0 , 7.0, 2.0 , 5.0, 3.0 , -3.0, 4.0 , 1.0 ], 2, 
                 53.0,  @flteps, "ibm example dot-03", 5)
        #    Example 4
        #        This example shows how to compute the dot product of a vector, x, 
        #        with a stride of 0, and a vector, y, with a stride of 1. The result 
        #        in DOTT is x1(y1+…...+yn).
        #        Function Reference and Input:
        #
        #                     N   X  INCX  Y  INCY
        #                     |   |   |    |   |
        #        DOTT = SDOT( 5 , X , 0  , Y , 1  )
        #
        #        X        =  (1.0, . , . , . , .)
        #        Y        =  (9.0, 8.0, 7.0, -6.0, 5.0)
        #
        #        Output:
        #
        #        DOTT     =  (1.0) × (9.0 + 8.0 + 7.0 - 6.0 + 5.0) = 23.0
        #
        test_dot(SingleBlas[ 1.0, 2.0, -3.0, 4.0, 5.0 ], 0, 
                 SingleBlas[ 9.0, 8.0, 7.0, -6.0, 5.0 ], 1, 
                 23.0,  @flteps, "ibm example dot-04", 5)
        #    Example 5
        #        This example shows how to compute the dot product of two vectors,
        #        x and y, with strides of 0. The result in DOTT is nx1y1.
        #        Function Reference and Input:
        #
        #                     N   X  INCX  Y  INCY
        #                     |   |   |    |   |
        #        DOTT = SDOT( 5 , X , 0  , Y , 0  )
        #
        #        X        =  (1.0, . , . , . , .)
        #        Y        =  (9.0, . , . , . , .)
        #
        #        Output:
        #
        #        DOTT     =  (5) × (1.0) × (9.0) = 45.0
        #
        test_dot(SingleBlas[ 1.0, 2.0, -3.0, 4.0, 5.0 ], 0, 
                 SingleBlas[ 9.0, 8.0, 7.0, -6.0, 5.0 ], 0, 
                 45.0,  @flteps, "ibm example dot-05", 5)
        #    Example 6
        #        This example shows how to compute the dot product of two vectors, 
        #        x and y, containing complex numbers, where x has a stride of 1, and y 
        #        has a stride greater than 1.
        #        Function Reference and Input:
        #
        #                      N   X  INCX  Y  INCY
        #                      |   |   |    |   |
        #        DOTT = CDOTU( 3 , X , 1  , Y , 2  )
        #
        #        X        =  ((1.0, 2.0), (3.0, -4.0), (-5.0, 6.0))
        #        Y        =  ((10.0, 9.0), . , (-6.0, 5.0), . , (2.0, 1.0))
        #
        #        Output:
        #
        #        DOTT     =  ((10.0 - 18.0 - 10.0) - (18.0 - 20.0 + 6.0),
        #                     (9.0 + 15.0 - 5.0) + (20.0 + 24.0 + 12.0))
        #                 =  (-22.0, 75.0)
        #
=begin
        test_dot(ComplexBlas[ Complex(1.0, 2.0), Complex(3.0, -4.0), Complex(-5.0, 6.0) ], 1, 
                 ComplexBlas[ Complex(10.0, 9.0), Complex(10.0, 9.0) , Complex(-6.0, 5.0), Complex(10.0, 9.0) , Complex(2.0, 1.0) ], 2, 
                 Complex(-22.0, 75.0),  @flteps, "ibm example dot-06", 3)
        #    Example 7
        #        This example shows how to compute the dot product of the conjugate of a
        #        vector, x, with vector y, both containing complex numbers, where x 
        #        has a stride of 1, and y has a stride greater than 1.
        #        Function Reference and Input:
        #
        #                      N   X  INCX  Y  INCY
        #                      |   |   |    |   |
        #        DOTT = CDOTC( 3 , X , 1  , Y , 2  )
        #
        #        X        =  ((1.0, 2.0), (3.0, -4.0), (-5.0, 6.0))
        #        Y        =  ((10.0, 9.0), . , (-6.0, 5.0), . , (2.0, 1.0))
        #
        #        Output:
        #
        #        DOTT     =  ((10.0 - 18.0 - 10.0) + (18.0 - 20.0 + 6.0),
        #                     (9.0  +  15.0  -  5.0)  -  (20.0  +  24.0  +  12.0))
        #                 =  (-14.0, -37.0)
        test_dotc(ComplexBlas[ Complex(1.0, 2.0), Complex(3.0, -4.0), Complex(-5.0, 6.0) ], 1, 
                 ComplexBlas[ Complex(10.0, 9.0), Complex(10.0, 9.0) , Complex(-6.0, 5.0), Complex(10.0, 9.0) , Complex(2.0, 1.0) ], 2, 
                 Complex(-14.0, -37.0),  @flteps, "ibm example dot-07", 3)
=end

  end
  
end

TestDot.new
