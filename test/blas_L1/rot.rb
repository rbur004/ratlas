require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestRot < TestBlas

	def initialize 
	  super
    gsl
  end
  
  def test_rot(x,incx,y,incy,c,s,x_expected,y_expected, error_bound, test_message)
    x.rot!(y, incx, incy, c, s);
    print_on_error( "rot x #{test_message}", x, x_expected, error_bound) 
    print_on_error( "rot y #{test_message}", y, y_expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    test_rot(SingleBlas[ -0.314 ], 1, SingleBlas[ -0.406 ], -1, 0.0, 0.0, SingleBlas[ 0.0 ], SingleBlas[ 0.0 ], @flteps, "gsl Test 558/559")
    test_rot(SingleBlas[ -0.314 ], 1, SingleBlas[ -0.406 ], -1, 0.866025403784, 0.5, SingleBlas[ -0.474932 ], SingleBlas[ -0.194606 ], @flteps, "gsl Test 560/561")
    test_rot(SingleBlas[ -0.314 ], 1, SingleBlas[ -0.406 ], -1, 0.0, -1.0, SingleBlas[ 0.406 ], SingleBlas[ -0.314 ], @flteps, "gsl Test 562/563")
    test_rot(SingleBlas[ -0.314 ], 1, SingleBlas[ -0.406 ], -1, -1.0, 0.0,  SingleBlas[ 0.314 ], SingleBlas[ 0.406 ],  @flteps, "gsl Test 564/565")

    test_rot(DoubleBlas[ -0.493 ], 1, DoubleBlas[ -0.014 ], -1, 0.0, 0.0,  DoubleBlas[ 0.0 ], DoubleBlas[ 0.0 ],  @dbleps, "gsl Test 566/567")
    test_rot(DoubleBlas[ -0.493 ], 1, DoubleBlas[ -0.014 ], -1, 0.866025403784, 0.5,  DoubleBlas[ -0.433950524066 ], DoubleBlas[ 0.234375644347 ],  @dbleps, "gsl Test 568/569")
    test_rot(DoubleBlas[ -0.493 ], 1, DoubleBlas[ -0.014 ], -1, 0.0, -1.0,  DoubleBlas[ 0.014 ], DoubleBlas[ -0.493 ],  @dbleps, "gsl Test 570/571")
    test_rot(DoubleBlas[ -0.493 ], 1, DoubleBlas[ -0.014 ], -1, -1.0, 0.0,  DoubleBlas[ 0.493 ], DoubleBlas[ 0.014 ],  @dbleps, "gsl Test 572/573")
    
    test_rot(SingleBlas[ -0.808 ], -1, SingleBlas[ -0.511 ], 1, 0.0, 0.0, SingleBlas[ 0.0 ], SingleBlas[ 0.0 ], @flteps, "gsl Test 574/575")
    test_rot(SingleBlas[ -0.808 ], -1, SingleBlas[ -0.511 ], 1, 0.866025403784, 0.5, SingleBlas[ -0.955249 ], SingleBlas[ -0.038539 ], @flteps, "gsl Test 576/577")
    test_rot(SingleBlas[ -0.808 ], -1, SingleBlas[ -0.511 ], 1, 0.0, -1.0, SingleBlas[ 0.511 ], SingleBlas[ -0.808 ], @flteps, "gsl Test 578/579")
    test_rot(SingleBlas[ -0.808 ], -1, SingleBlas[ -0.511 ], 1, -1.0, 0.0,  SingleBlas[ 0.808 ], SingleBlas[ 0.511 ],  @flteps, "gsl Test 580/581")

    test_rot(DoubleBlas[ -0.176 ], -1, DoubleBlas[ -0.165 ], 1, 0.0, 0.0,  DoubleBlas[ 0.0 ], DoubleBlas[ 0.0 ],  @dbleps, "gsl Test 582/583")
    test_rot(DoubleBlas[ -0.176 ], -1, DoubleBlas[ -0.165 ], 1, 0.866025403784, 0.5,  DoubleBlas[ -0.234920471066 ], DoubleBlas[ -0.0548941916244 ],  @dbleps, "gsl Test 584/585")
    test_rot(DoubleBlas[ -0.176 ], -1, DoubleBlas[ -0.165 ], 1, 0.0, -1.0,  DoubleBlas[ 0.165 ], DoubleBlas[ -0.176 ],  @dbleps, "gsl Test 586/587")
    test_rot(DoubleBlas[ -0.176 ], -1, DoubleBlas[ -0.165 ], 1, -1.0, 0.0,  DoubleBlas[ 0.176 ], DoubleBlas[ 0.165 ],  @dbleps, "gsl Test 588/589")
    
    test_rot(SingleBlas[ -0.201 ], -1, SingleBlas[ 0.087 ], -1, 0.0, 0.0, SingleBlas[ 0.0 ], SingleBlas[ 0.0 ], @flteps, "gsl Test 590/591")
    test_rot(SingleBlas[ -0.201 ], -1, SingleBlas[ 0.087 ], -1, 0.866025403784, 0.5, SingleBlas[ -0.130571 ], SingleBlas[ 0.175844 ], @flteps, "gsl Test 592/593")
    test_rot(SingleBlas[ -0.201 ], -1, SingleBlas[ 0.087 ], -1, 0.0, -1.0, SingleBlas[ -0.087 ], SingleBlas[ -0.201 ], @flteps, "gsl Test 594/595")
    test_rot(SingleBlas[ -0.201 ], -1, SingleBlas[ 0.087 ], -1, -1.0, 0.0,  SingleBlas[ 0.201 ], SingleBlas[ -0.087 ],  @flteps, "gsl Test 596/597")

    test_rot(DoubleBlas[ -0.464 ], -1, DoubleBlas[ 0.7 ], -1, 0.0, 0.0,  DoubleBlas[ 0.0 ], DoubleBlas[ 0.0 ],  @dbleps, "gsl Test 598/599")
    test_rot(DoubleBlas[ -0.464 ], -1, DoubleBlas[ 0.7 ], -1, 0.866025403784, 0.5,  DoubleBlas[ -0.051835787356 ], DoubleBlas[ 0.838217782649 ],  @dbleps, "gsl Test 600/601")
    test_rot(DoubleBlas[ -0.464 ], -1, DoubleBlas[ 0.7 ], -1, 0.0, -1.0,  DoubleBlas[ -0.7 ], DoubleBlas[ -0.464 ],  @dbleps, "gsl Test 602/603")
    test_rot(DoubleBlas[ -0.464 ], -1, DoubleBlas[ 0.7 ], -1, -1.0, 0.0,  DoubleBlas[ 0.464 ], DoubleBlas[ -0.7 ],  @dbleps, "gsl Test 604/605")
  end

end

TestRot.new

