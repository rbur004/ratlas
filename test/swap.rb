require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestSwap < TestBlas

	def initialize 
	  super
    gsl
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
    test_swap(ComplexBlas[ Complex(-0.316, -0.529) ], 1, ComplexBlas[ Complex(-0.313, 0.363) ], -1, ComplexBlas[ Complex(-0.313, 0.363) ], ComplexBlas[Complex(-0.316, -0.529)],  Complex(@flteps,@flteps), "gsl Test 92/93")
    test_swap(DoubleComplexBlas[ Complex(0.512, -0.89) ], 1, DoubleComplexBlas[ Complex(-0.225, -0.511) ], -1, DoubleComplexBlas[Complex(-0.225, -0.511)],   DoubleComplexBlas[Complex(0.512, -0.89) ], Complex(@dbleps,@dbleps), "gsl Test 94/95")

    test_swap(SingleBlas[ 0.336 ], -1, SingleBlas[ -0.431 ], 1, SingleBlas[ -0.431 ], SingleBlas[ 0.336 ], @flteps, "gsl Test 96/97")
    test_swap(DoubleBlas[ 0.764 ], -1, DoubleBlas[ -0.293 ], 1,  DoubleBlas[ -0.293 ], DoubleBlas[0.764 ],  @dbleps, "gsl Test 98/99")
    test_swap(ComplexBlas[ Complex(-0.239, 0.361) ], -1, ComplexBlas[ Complex(0.149, 0.347) ], 1, ComplexBlas[ Complex(0.149, 0.347) ], ComplexBlas[Complex(-0.239, 0.361)],  Complex(@flteps,@flteps), "gsl Test 100/101")
    test_swap(DoubleComplexBlas[ Complex(-0.171, -0.936) ], -1, DoubleComplexBlas[ Complex( 0.495, -0.835) ], 1, DoubleComplexBlas[Complex( 0.495, -0.835)],   DoubleComplexBlas[Complex(-0.171, -0.936) ], Complex(@dbleps,@dbleps), "gsl Test 102/103")
 
    test_swap(SingleBlas[ -0.405 ], -1, SingleBlas[ -0.213 ], -1, SingleBlas[ -0.213 ], SingleBlas[ -0.405 ], @flteps, "gsl Test 104/105")
    test_swap(DoubleBlas[ -0.761 ], -1, DoubleBlas[ -0.585 ], -1,  DoubleBlas[ -0.585 ], DoubleBlas[-0.761 ],  @dbleps, "gsl Test 106/107")
    test_swap(ComplexBlas[ Complex(0.853, 0.146) ], -1, ComplexBlas[ Complex( 0.009, -0.178) ], -1, ComplexBlas[ Complex( 0.009, -0.178) ], Complex(0.853, 0.146),  Complex(@flteps,@flteps), "gsl Test 108/109")
    test_swap(DoubleComplexBlas[ Complex(-0.228, 0.386) ], -1, DoubleComplexBlas[ Complex(  0.988, -0.084) ], -1, DoubleComplexBlas[Complex( 0.988, -0.084)],   DoubleComplexBlas[Complex(-0.228, 0.386) ], Complex(@dbleps,@dbleps), "gsl Test 110/111")
  end

end

TestSwap.new
