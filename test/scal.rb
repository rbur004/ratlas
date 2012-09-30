require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestScal < TestBlas

	def initialize 
	  super
    gsl
  end
  
  def test_scal(x,incx, alpha, x_expected, error_bound, test_message)
    x.scal(alpha, incx);
    print_on_error( "scal #{test_message}", x, x_expected, error_bound) 
  end
  
  def gsl
    #Simple tests based on the values used in the gsl cblas test suite.
    test_scal(SingleBlas[ 0.651 ], -1, 0.0, SingleBlas[ 0.651 ], @flteps, "gsl Test 112")
    test_scal(SingleBlas[ 0.651 ], -1, 0.1, SingleBlas[ 0.651 ], @flteps, "gsl Test 113")
    test_scal(SingleBlas[ 0.651 ], -1, 1.0, SingleBlas[ 0.651 ], @flteps, "gsl Test 114")

    test_scal(DoubleBlas[  0.686 ], -1, 0.0, DoubleBlas[ 0.686 ],  @dbleps, "gsl Test 115")
    test_scal(DoubleBlas[  0.686 ], -1, 0.1, DoubleBlas[ 0.686 ],  @dbleps, "gsl Test 116")
    test_scal(DoubleBlas[  0.686 ], -1, 1.0, DoubleBlas[ 0.686 ],  @dbleps, "gsl Test 117")

    test_scal(ComplexBlas[ Complex(0.986, -0.775) ], -1, Complex(0.0, 0.0) , ComplexBlas[Complex(0.986, -0.775)],  Complex(@flteps,@flteps), "gsl Test 118")
    test_scal(ComplexBlas[ Complex(0.986, -0.775) ], -1, Complex(0.1,0.0) , ComplexBlas[Complex(0.986, -0.775)],  Complex(@flteps,@flteps), "gsl Test 119")
    test_scal(ComplexBlas[ Complex(0.986, -0.775) ], -1, Complex(1.0, 0.0) , ComplexBlas[Complex(0.986, -0.775)],  Complex(@flteps,@flteps), "gsl Test 120")
    test_scal(ComplexBlas[ Complex(0.986, -0.775) ], -1, Complex(0.0, 0.1) , ComplexBlas[Complex(0.986, -0.775)],  Complex(@flteps,@flteps), "gsl Test 121")
    test_scal(ComplexBlas[ Complex(0.986, -0.775) ], -1, Complex(0.1, 0.2) , ComplexBlas[Complex(0.986, -0.775)],  Complex(@flteps,@flteps), "gsl Test 122")
    test_scal(ComplexBlas[ Complex(0.986, -0.775) ], -1, Complex(1.0, 0.3) , ComplexBlas[Complex(0.986, -0.775)],  Complex(@flteps,@flteps), "gsl Test 123")

    test_scal(DoubleComplexBlas[ Complex(0.454, -0.478 ) ], -1, Complex(0.0, 0.0),  DoubleComplexBlas[Complex(0.454, -0.478 )], Complex(@dbleps,@dbleps), "gsl Test 124")
    test_scal(DoubleComplexBlas[ Complex(0.454, -0.478 ) ], -1, Complex(0.1, 0.0),  DoubleComplexBlas[Complex(0.454, -0.478 )], Complex(@dbleps,@dbleps), "gsl Test 125")
    test_scal(DoubleComplexBlas[ Complex(0.454, -0.478 ) ], -1, Complex(1.0, 0.0),  DoubleComplexBlas[Complex(0.454, -0.478 )], Complex(@dbleps,@dbleps), "gsl Test 126")
    test_scal(DoubleComplexBlas[ Complex(0.454, -0.478 ) ], -1, Complex(0.0, 0.1),  DoubleComplexBlas[Complex(0.454, -0.478 )], Complex(@dbleps,@dbleps), "gsl Test 127")
    test_scal(DoubleComplexBlas[ Complex(0.454, -0.478 ) ], -1, Complex(0.1, 0.2),  DoubleComplexBlas[Complex(0.454, -0.478 )], Complex(@dbleps,@dbleps), "gsl Test 128")
    test_scal(DoubleComplexBlas[ Complex(0.454, -0.478 ) ], -1, Complex(1.0, 0.3),  DoubleComplexBlas[Complex(0.454, -0.478 )], Complex(@dbleps,@dbleps), "gsl Test 129")

    #2D Vector, incx 1
    test_scal(SingleBlas[ 0.389, -0.236 ], 1, 0.0, SingleBlas[ 0.0, 0.0 ], @flteps, "gsl Test 130")
    test_scal(SingleBlas[ 0.389, -0.236 ], 1, 0.1, SingleBlas[ 0.0389, -0.0236 ], @flteps, "gsl Test 131")
    test_scal(SingleBlas[ 0.389, -0.236 ], 1, 1.0, SingleBlas[ 0.389, -0.236 ], @flteps, "gsl Test 132")

    test_scal(DoubleBlas[  -0.429, -0.183 ], 1, 0.0, DoubleBlas[ 0.0, 0.0 ],  @dbleps, "gsl Test 133")
    test_scal(DoubleBlas[  -0.429, -0.183 ], 1, 0.1, DoubleBlas[ -0.0429, -0.0183 ],  @dbleps, "gsl Test 134")
    test_scal(DoubleBlas[  -0.429, -0.183 ], 1, 1.0, DoubleBlas[ -0.429, -0.183 ],  @dbleps, "gsl Test 135")

    test_scal(ComplexBlas[ Complex(-0.603, 0.239), Complex(0.339, -0.58) ], 1, Complex(0.0, 0.0) , ComplexBlas[Complex(0.0, 0.0),Complex(0.0, 0.0)],  Complex(@flteps,@flteps), "gsl Test 136")
    test_scal(ComplexBlas[ Complex(-0.603, 0.239), Complex(0.339, -0.58) ], 1, Complex(0.1,0.0) , ComplexBlas[Complex(-0.0603, 0.0239), Complex(0.0339, -0.058)],  Complex(@flteps,@flteps), "gsl Test 137")
    test_scal(ComplexBlas[ Complex(-0.603, 0.239), Complex(0.339, -0.58) ], 1, Complex(1.0, 0.0) , ComplexBlas[Complex(-0.603, 0.239), Complex(0.339, -0.58)],  Complex(@flteps,@flteps), "gsl Test 138")
    test_scal(ComplexBlas[ Complex(-0.603, 0.239), Complex(0.339, -0.58) ], 1, Complex(0.0, 0.1) , ComplexBlas[Complex(-0.0239, -0.0603), Complex(0.058, 0.0339)],  Complex(@flteps,@flteps), "gsl Test 139")
    test_scal(ComplexBlas[ Complex(-0.603, 0.239), Complex(0.339, -0.58) ], 1, Complex(0.1, 0.2) , ComplexBlas[Complex(-0.1081, -0.0967), Complex(0.1499, 0.0098)],  Complex(@flteps,@flteps), "gsl Test 140")
    test_scal(ComplexBlas[ Complex(-0.603, 0.239), Complex(0.339, -0.58) ], 1, Complex(1.0, 0.3) , ComplexBlas[Complex(-0.6747, 0.0581), Complex(0.513, -0.4783)],  Complex(@flteps,@flteps), "gsl Test 141")

    test_scal(DoubleComplexBlas[ Complex(-0.956, 0.613 ), Complex( 0.443, 0.503) ], 1, Complex(0.0, 0.0),  DoubleComplexBlas[Complex(0.0, 0.0 ), Complex(0.0, 0.0 ) ], Complex(@dbleps,@dbleps), "gsl Test 142")
    test_scal(DoubleComplexBlas[Complex(-0.956, 0.613 ), Complex( 0.443, 0.503) ], 1, Complex(0.1, 0.0),  DoubleComplexBlas[Complex(-0.0956, 0.0613 ), Complex( 0.0443, 0.0503 )], Complex(@dbleps,@dbleps), "gsl Test 143")
    test_scal(DoubleComplexBlas[ Complex(-0.956, 0.613 ), Complex( 0.443, 0.503) ], 1, Complex(1.0, 0.0),  DoubleComplexBlas[Complex(-0.956, 0.613 ), Complex( 0.443, 0.503)], Complex(@dbleps,@dbleps), "gsl Test 144")
    test_scal(DoubleComplexBlas[ Complex(-0.956, 0.613 ), Complex( 0.443, 0.503) ], 1, Complex(0.0, 0.1),  DoubleComplexBlas[Complex(-0.0613, -0.0956 ), Complex( -0.0503, 0.0443)], Complex(@dbleps,@dbleps), "gsl Test 145")
    test_scal(DoubleComplexBlas[ Complex(-0.956, 0.613 ), Complex( 0.443, 0.503) ], 1, Complex(0.1, 0.2),  DoubleComplexBlas[Complex(-0.2182, -0.1299), Complex( -0.0563, 0.1389 )], Complex(@dbleps,@dbleps), "gsl Test 146")
    test_scal(DoubleComplexBlas[ Complex(-0.956, 0.613 ), Complex( 0.443, 0.503) ], 1, Complex(1.0, 0.3),  DoubleComplexBlas[Complex(-1.1399, 0.3262), Complex( 0.2921, 0.6359 )], Complex(@dbleps,@dbleps), "gsl Test 147")

    #2D Vector, incx -1
    test_scal(SingleBlas[ 0.629, -0.419 ], -1, 0.0, SingleBlas[ 0.629, -0.419 ], @flteps, "gsl Test 148")
    test_scal(SingleBlas[ 0.629, -0.419 ], -1, 0.1, SingleBlas[ 0.629, -0.419 ], @flteps, "gsl Test 149")
    test_scal(SingleBlas[ 0.629, -0.419 ], -1, 1.0, SingleBlas[ 0.629, -0.419 ], @flteps, "gsl Test 150")

    test_scal(DoubleBlas[   0.398, -0.656 ], -1, 0.0, DoubleBlas[  0.398, -0.656 ],  @dbleps, "gsl Test 151")
    test_scal(DoubleBlas[   0.398, -0.656 ], -1, 0.1, DoubleBlas[  0.398, -0.656 ],  @dbleps, "gsl Test 152")
    test_scal(DoubleBlas[   0.398, -0.656 ], -1, 1.0, DoubleBlas[  0.398, -0.656 ],  @dbleps, "gsl Test 153")

    test_scal(ComplexBlas[ Complex(0.736, 0.331), Complex(-0.318, 0.622) ], -1, Complex(0.0, 0.0) , ComplexBlas[Complex(0.736, 0.331), Complex(-0.318, 0.622)],  Complex(@flteps,@flteps), "gsl Test 154")
    test_scal(ComplexBlas[ Complex(0.736, 0.331), Complex(-0.318, 0.622) ], -1, Complex(0.1,0.0) , ComplexBlas[Complex(0.736, 0.331), Complex(-0.318, 0.622)],  Complex(@flteps,@flteps), "gsl Test 155")
    test_scal(ComplexBlas[ Complex(0.736, 0.331), Complex(-0.318, 0.622) ], -1, Complex(1.0, 0.0) , ComplexBlas[Complex(0.736, 0.331), Complex(-0.318, 0.622)],  Complex(@flteps,@flteps), "gsl Test 156")
    test_scal(ComplexBlas[ Complex(0.736, 0.331), Complex(-0.318, 0.622) ], -1, Complex(0.0, 0.1) , ComplexBlas[Complex(0.736, 0.331), Complex(-0.318, 0.622)],  Complex(@flteps,@flteps), "gsl Test 157")
    test_scal(ComplexBlas[ Complex(0.736, 0.331), Complex(-0.318, 0.622) ], -1, Complex(0.1, 0.2) , ComplexBlas[Complex(0.736, 0.331), Complex(-0.318, 0.622)],  Complex(@flteps,@flteps), "gsl Test 158")
    test_scal(ComplexBlas[ Complex(0.736, 0.331), Complex(-0.318, 0.622) ], -1, Complex(1.0, 0.3) , ComplexBlas[Complex(0.736, 0.331), Complex(-0.318, 0.622)],  Complex(@flteps,@flteps), "gsl Test 159")

    test_scal(DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], -1, Complex(0.0, 0.0),  DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], Complex(@dbleps,@dbleps), "gsl Test 160")
    test_scal(DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], -1, Complex(0.1, 0.0),  DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], Complex(@dbleps,@dbleps), "gsl Test 161")
    test_scal(DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], -1, Complex(1.0, 0.0),  DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], Complex(@dbleps,@dbleps), "gsl Test 162")
    test_scal(DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], -1, Complex(0.0, 0.1),  DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], Complex(@dbleps,@dbleps), "gsl Test 163")
    test_scal(DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], -1, Complex(0.1, 0.2),  DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], Complex(@dbleps,@dbleps), "gsl Test 164")
    test_scal(DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], -1, Complex(1.0, 0.3),  DoubleComplexBlas[ Complex( 0.521, -0.811), Complex(0.556, -0.147) ], Complex(@dbleps,@dbleps), "gsl Test 165")
  end

end

TestScal.new

