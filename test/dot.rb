require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestDot < TestBlas

  def initialize
    super
    gsl
  end
  
  def test_dotsds(x,incx,y,incy,alpha , expected, error_bound, test_message)
    result = x.dotsds(y, alpha,  incx, incy);
    print_on_error( "dotsds #{test_message}", result, expected, error_bound) 
  end
  
  def test_dot(x,incx,y,incy , expected, error_bound, test_message)
    result = x.dot( y,  incx, incy);
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
    test_dot(ComplexBlas[ Complex(0.474, -0.27) ], 1, ComplexBlas[ Complex(-0.144, -0.392) ], -1, Complex(-0.174096, -0.146928),  Complex(@flteps,@flteps), "gsl Test 12")
    test_dotc(ComplexBlas[ Complex(0.474, -0.27) ], 1, ComplexBlas[ Complex(-0.144, -0.392) ], -1, Complex(0.037584, -0.224688),  Complex(@flteps,@flteps), "gsl Test 13")
    test_dot(DoubleComplexBlas[ Complex(-0.87, -0.631) ], 1, DoubleComplexBlas[ Complex(-0.7, -0.224) ], -1, Complex(0.467656, 0.63658),  Complex(@dbleps,@dbleps), "gsl Test 14")
    test_dotc(DoubleComplexBlas[ Complex(-0.87, -0.631) ], 1, DoubleComplexBlas[ Complex(-0.7, -0.224) ], -1, Complex(0.750344, -0.24682),  Complex(@dbleps,@dbleps), "gsl Test 15")

    test_dot(SingleBlas[ -0.457 ], -1, SingleBlas[ 0.839 ], 1, -0.383423,  @flteps, "gsl Test 16")
    test_dot(DoubleBlas[  0.949 ], -1, DoubleBlas[ -0.873 ], 1, -0.828477,  @dbleps, "gsl Test 17")
    test_dot(ComplexBlas[ Complex(0.852, -0.045) ], -1, ComplexBlas[ Complex(0.626, -0.164) ], 1, Complex(0.525972, -0.167898), Complex(@flteps,@flteps), "gsl Test 18")
    test_dotc(ComplexBlas[ Complex(0.852, -0.045) ], -1, ComplexBlas[ Complex(0.626, -0.164) ], 1, Complex(0.540732, -0.111558),  Complex(@flteps,@flteps), "gsl Test 19")
    test_dot(DoubleComplexBlas[ Complex(-0.786, -0.341) ], -1, DoubleComplexBlas[ Complex(-0.271, -0.896) ], 1, Complex(-0.09253, 0.796667),  Complex(@dbleps,@dbleps), "gsl Test 20")
    test_dotc(DoubleComplexBlas[ Complex(-0.786, -0.341) ], -1, DoubleComplexBlas[ Complex(-0.271, -0.896) ], 1, Complex(0.518542, 0.611845),  Complex(@dbleps,@dbleps), "gsl Test 21")

    test_dot(SingleBlas[ -0.088 ], -1, SingleBlas[ -0.165 ], -1, 0.01452,  @flteps, "gsl Test 22")
    test_dot(DoubleBlas[  -0.434 ], -1, DoubleBlas[ -0.402 ], -1, 0.174468,  @dbleps, "gsl Test 23")
    test_dot(ComplexBlas[ Complex(-0.347, 0.899) ], -1, ComplexBlas[ Complex(-0.113, -0.858) ], -1, Complex(0.810553, 0.196139),  Complex(@flteps,@flteps), "gsl Test 24")
    test_dotc(ComplexBlas[ Complex(-0.347, 0.899) ], -1, ComplexBlas[ Complex(-0.113, -0.858) ], -1, Complex(-0.732131, 0.399313),  Complex(@flteps,@flteps), "gsl Test 25")
    test_dot(DoubleComplexBlas[ Complex(-0.897, -0.204 ) ], -1, DoubleComplexBlas[ Complex(-0.759, 0.557) ], -1, Complex(0.794451, -0.344793),  Complex(@dbleps,@dbleps), "gsl Test 26")
    test_dotc(DoubleComplexBlas[ Complex(-0.897, -0.204 ) ], -1, DoubleComplexBlas[ Complex(-0.759, 0.557) ], -1, Complex(0.567195, -0.654465),  Complex(@dbleps,@dbleps), "gsl Test 27")
  end
  
end

TestDot.new
