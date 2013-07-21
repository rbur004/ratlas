require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestRot < TestBlas

	def initialize 
	  super
    gsl
    ibm_examples
  end
  
  def test_rot(x,incx,y,incy,c,s,x_expected,y_expected, error_bound, test_message)
    #puts x.class, y.class, incx.class, incy.class, c.class, s.class
    puts x, y, incx, incy, c, s
    x.rot!(y, c, s, incx, incy);
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
  
  def ibm_examples
    #Example 1
    #    This example shows how to apply a real plane rotation to real 
    #    vectors x and y having positive strides.
    #    Call Statement and Input:

    #               N   X  INCX  Y  INCY   C    S
    #               |   |   |    |   |     |    |
    #    CALL SROT( 5 , X , 1  , Y , 2  , 0.5 , S )

    #    X        =  (1.0, 2.0, 3.0, 4.0, 5.0)
    #    Y        =  (-1.0, . , -2.0, . , -3.0, . , -4.0, . , -5.0)

    #    Output:

    #    X        =  (-0.366, -0.732, -1.098, -1.464, -1.830)
    #    Y        =  (-1.366, -2.732, -4.098, -5.464, -6.830)

    #Example 2
    #    This example shows how to apply a real plane rotation to real 
    #    vectors x and y having strides of opposite sign.
    #    Call Statement and Input:

    #               N   X  INCX  Y   INCY   C    S
    #               |   |   |    |    |     |    |
    #    CALL SROT( 5 , X , 1  , Y , -1  , 0.5 , S )

    #    X        =  (1.0, 2.0, 3.0, 4.0, 5.0)
    #    Y        =  (-5.0, -4.0, -3.0, -2.0, -1.0)

    #    Output:

    #    X        =(same as output X in Example 1)
    #    Y        =  (-6.830, -5.464, -4.098, -2.732, -1.366)
    #Example 3
    #    This example shows how scalar values in vectors x and y can be 
    #    processed by specifying 0 strides and the number of elements to be 
    #    processed, n, equal to 1.
    #    Call Statement and Input:

    #               N   X  INCX  Y  INCY   C    S
    #               |   |   |    |   |     |    |
    #    CALL SROT( 1 , X , 0  , Y , 0  , 0.5 , S )

    #    X        =  (1.0)
    #    Y        =  (-1.0)

    #    Output:

    #    X        =  (-0.366)
    #    Y        =  (-1.366)

    #Example 4
    #    This example shows how to apply a complex plane rotation to complex
    #     vectors x and y having positive strides.
    #    Call Statement and Input:

    #               N   X  INCX  Y  INCY   C    S
    #               |   |   |    |   |     |    |
    #    CALL CROT( 3 , X , 1  , Y , 2  , 0.5 , S )

    #    X        =  ((1.0, 2.0), (2.0, 3.0), (3.0, 4.0))
    #    Y        =  ((-1.0, 5.0), . , (-2.0, 4.0), . , (-3.0, 3.0))
    #    S        =  (0.75, 0.50)

    #    Output:

    #    X        =  ((-2.750, 4.250), (-2.500, 3.500), (-2.250, 2.750))
    #    Y        =  ((-2.250, 1.500), . , (-4.000, 0.750), . ,
    #                 (-5.750, 0.000))

    #Example 5
    #    This example shows how to apply a real plane rotation to complex vectors 
    #    x and y having positive strides.
    #    Call Statement and Input:

    #                N   X  INCX  Y  INCY   C    S
    #                |   |   |    |   |     |    |
    #    CALL CSROT( 3 , X , 1  , Y , 2  , 0.5 , S )

    #    X        =  ((1.0, 2.0), (2.0, 3.0), (3.0, 4.0))
    #    Y        =  ((-1.0, 5.0), . , (-2.0, 4.0), . , (-3.0, 3.0))
    
    #    Output:

    #    X        =  ((-0.366, 5.330), (-0.732, 4.964), (-1.098, 4.598))
    #    Y        =  ((-1.366, 0.768), . , (-2.732, -0.598), . ,
    #                 (-4.098, -1.964))

    end

end

TestRot.new

