require 'ratlas'
include RAtlas
require 'complex' 
include Math
require '../testblas.rb'

class TestAxpy < TestBlas

	def initialize 
	  super
	  simple_sums
    gsl
  end
  
  def test_axpy(x,incx,y,incy,alpha , x_expected, y_expected, error_bound, test_message)
    x.axpy!(y, alpha, incx, incy)
    print_on_error( "axpy #{test_message}", x, x_expected, error_bound) 
    print_on_error( "axpy #{test_message}", y, y_expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    test_axpy(SingleBlas[ 0.018 ], 1, SingleBlas[-0.417], -1, 0.0, SingleBlas[ 0.018 ], SingleBlas[-0.417], @flteps, "gsl Test 64")
    test_axpy(DoubleBlas[ 0.071 ], 1, DoubleBlas[-0.888], -1, 0.0, DoubleBlas[ 0.071 ], DoubleBlas[-0.888], @dbleps, "gsl Test 65")
#    test_axpy(ComplexBlas[ Complex(-0.542, -0.362) ], 1, ComplexBlas[ Complex(-0.459, -0.433) ], -1, Complex(1.0, 0.0),  ComplexBlas[ Complex(-0.542, -0.362) ], ComplexBlas[ Complex(-1.001, -0.795) ], Complex(@flteps,@flteps), "gsl Test 66")
#    test_axpy(DoubleComplexBlas[ Complex(0.003, -0.514) ], 1, DoubleComplexBlas[ Complex(-0.529, 0.743) ], -1, Complex(-1.0, 0.0),  DoubleComplexBlas[ Complex(0.003, -0.514) ], DoubleComplexBlas[ Complex(-0.529, 0.743) ], Complex(@dbleps,@dbleps), "gsl Test 67")
 
    test_axpy(SingleBlas[ 0.771 ], -1, SingleBlas[0.507], 1, 0.1, SingleBlas[ 0.771 ], SingleBlas[0.5841], @flteps, "gsl Test 68")
    test_axpy(DoubleBlas[ 0.029 ], -1, DoubleBlas[-0.992], 1, -0.3, DoubleBlas[ 0.029 ], DoubleBlas[-1.0007], @dbleps, "gsl Test 69")
#    test_axpy(ComplexBlas[ Complex(0.194, -0.959) ], -1, ComplexBlas[ Complex(0.096, 0.032) ], 1, Complex(-0.3, 0.1),  ComplexBlas[ Complex(0.194, -0.959) ], ComplexBlas[ Complex(0.1337, 0.3391) ], Complex(@flteps,@flteps), "gsl Test 70")
#    test_axpy(DoubleComplexBlas[ Complex(0.776, -0.671) ], -1, DoubleComplexBlas[ Complex(0.39, 0.404) ], 1, Complex(0.0, 1.0),  DoubleComplexBlas[ Complex(0.776, -0.671) ], DoubleComplexBlas[ Complex(1.061, 1.18) ], Complex(@dbleps,@dbleps), "gsl Test 71")
 
    test_axpy(SingleBlas[ 0.647 ], -1, SingleBlas[0.016], -1, 1.0, SingleBlas[ 0.647 ], SingleBlas[0.663], @flteps, "gsl Test 72")
    test_axpy(DoubleBlas[ -0.558 ], -1, DoubleBlas[0.308], -1, -1.0, DoubleBlas[ -0.558 ], DoubleBlas[0.866], @dbleps, "gsl Test 73")
#    test_axpy(ComplexBlas[ Complex(0.899, -0.624) ], -1, ComplexBlas[ Complex(0.155, -0.33) ], -1, Complex(-0.3, 0.1),  ComplexBlas[ Complex(0.899, -0.624) ], ComplexBlas[ Complex(-0.0523, -0.0529) ], Complex(@flteps,@flteps), "gsl Test 74")
#    test_axpy(DoubleComplexBlas[ Complex(-0.451, 0.768) ], -1, DoubleComplexBlas[ Complex(0.007, 0.732) ], -1, Complex(0.0, 1.0),  DoubleComplexBlas[ Complex(-0.451, 0.768) ], DoubleComplexBlas[ Complex(-0.761, 0.281) ], Complex(@dbleps,@dbleps), "gsl Test 75")
  end

  def simple_sums
    #Don't use the built in array checking as it requires these to work first.
    r = SingleBlas[ 0.018 ] - SingleBlas[ 0.018 ]
    print_on_error( "1D -", r[0,0], 0.0 , @flteps)
    r = SingleBlas[ 0.018 ] + SingleBlas[ 0.018 ]
    print_on_error( "1D +", r[0,0], 0.036, @flteps)
    r = SingleBlas[ 0.018 ] * 2
    print_on_error( "1D *", r[0,0], 0.036, @flteps)
    
    #2D
    r = SingleBlas[ 0.018, 0.16 ] - SingleBlas[ 0.018, 0.16 ]
    print_on_error( "2D - 0", r[0], 0.0, @flteps)
    print_on_error( "2D - 1", r[1], 0.0, @flteps)
    r = SingleBlas[ 0.018, 0.16 ] + SingleBlas[ 0.018, 0.16 ]
    print_on_error( "2D + 0", r[0], 0.036, @flteps)
    print_on_error( "2D + 0", r[1], 0.32, @flteps)
    r = SingleBlas[ 0.018, 0.16 ] * 2
    print_on_error( "2D * 0", r[0], 0.036, @flteps)
    print_on_error( "2D * 0", r[1], 0.32, @flteps)
  end
end 

TestAxpy.new
