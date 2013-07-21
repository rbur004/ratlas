require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestRotg < TestBlas
  
	def initialize 
	  super
	  local_tests
    ibm_examples
  end
  
  def test_rot(a,b,x_expected,y_expected, error_bound, test_message)
    #puts x.class, y.class, incx.class, incy.class, c.class, s.class
    puts x, y, incx, incy, c, s
    x.rotg!(y, c, s, incx, incy);
    print_on_error( "rot x #{test_message}", x, x_expected, error_bound) 
    print_on_error( "rot y #{test_message}", y, y_expected, error_bound) 
  end
  
  def test_srotg(a,b,r_expected,z_expected, c_expected, s_expected, error_bound, test_message, c = nil, s = nil)
    r = c == nil ? Srotg.new(a, b) : Srotg.new(a,b,c,s);
    print_on_error( "rot r #{test_message}", r.r, r_expected, error_bound) 
    print_on_error( "rot z #{test_message}", r.z, z_expected, error_bound) 
    print_on_error( "rot c #{test_message}", r.c, c_expected, error_bound) 
    print_on_error( "rot s #{test_message}", r.s, s_expected, error_bound) 
  end

  def local_tests
    #Giving all params means the srotg call is not made.
    #Tests given params can be recalled correctly
    test_srotg( 1.0, 2.0, 1.0, 2.0, 3.0, 4.0, 0.0, "Check giving all params", 3.0, 4.0)
  end
  
  def ibm_examples
        #  Example 1
        #      This example shows the construction of a real Givens plane rotation, 
        #      where r is 0.
        #      Call Statement and Input:
        #  
        #                   A     B    C   S
        #                   |     |    |   |
        #      CALL SROTG( 0.0 , 0.0 , C , S )
        #  
        #      Output:
        #  
        #      A        =  0.0
        #      B        =  0.0
        #      C        =  1.0
        #      S        =  0.0
        #  
        test_srotg( 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, "ibm example rotg-01")
        #  Example 2
        #      This example shows the construction of a real Givens plane rotation, 
        #      where c is 0.
        #      Call Statement and Input:
        #  
        #                   A     B    C   S
        #                   |     |    |   |
        #      CALL SROTG( 0.0 , 2.0 , C , S )
        #  
        #      Output:
        #  
        #      A        =  2.0
        #      B        =  1.0
        #      C        =  0.0
        #      S        =  1.0
        #  
        test_srotg( 0.0, 2.0, 2.0, 1.0, 0.0, 1.0, 0.0, "ibm example rotg-02")
        #  Example 3
        #      This example shows the construction of a real Givens plane rotation, 
        #      where |b| > greater than |a|.
        #      Call Statement and Input:
        #  
        #                   A      B    C   S
        #                   |      |    |   |
        #      CALL SROTG( 6.0 , -8.0 , C , S )
        #  
        #      Output:
        #  
        #      A        =  -10.0
        #      B        =  -1.666
        #      C        =  -0.6
        #      S        =  0.8
        #  
        test_srotg( 6.0, -8.0, -10.0, -1.666, -0.6, 0.8, 0.001, "ibm example rotg-03")
        #  Example 4
        #      This example shows the construction of a real Givens plane rotation, 
        #      where |a| > greater than |b|.
        #      Call Statement and Input:
        #  
        #                   A     B    C   S
        #                   |     |    |   |
        #      CALL SROTG( 8.0 , 6.0 , C , S )
        #  
        #      Output:
        #  
        #      A        =  10.0
        #      B        =  0.6
        #      C        =  0.8
        #      S        =  0.6
        #  
        test_srotg( 8.0, 6.0, 10.0, 0.6, 0.8, 0.6, 0.001, "ibm example rotg-04")
=begin
        #  Example 5
        #      This example shows the construction of a complex Givens plane rotation, 
        #      where |a| = equal to 0.
        #      Call Statement and Input:
        #  
        #                  A   B   C   S
        #                  |   |   |   |
        #      CALL CROTG( A , B , C , S )
        #       
        #      A        =  (0.0, 0.0)
        #      B        =  (1.0, 0.0)
        #       
        #  
        #      Output:
        #  
        #      A        =  (1.0, 0.0)
        #      C        =  0.0
        #      S        =  (1.0, 0.0)
        #  
        puts "example 5 crotg()"
        sr5 = Crotg.new(Complex(0.0, 0.0), Complex(1.0, 0.0))
        puts sr5.r, sr5.z, sr5.c, sr5.s
        #  Example 6
        #      This example shows the construction of a complex Givens plane rotation,
        #      where |a| ≠ not equal to 0.
        #      Call Statement and Input:
        #  
        #                  A   B   C   S
        #                  |   |   |   |
        #      CALL CROTG( A , B , C , S )
        #       
        #      A        =  (3.0, 4.0)
        #      B        =  (4.0, 6.0)
        #  
        #      Output:
        #  
        #      A        =  (5.26, 7.02)
        #      C        =  0.57
        #      S        =  (0.82, -0.05)
        puts "example 6 zrotg()"
        sr6 = Zrotg.new(Complex(3.0, 4.0), Complex(4.0, 6.0))
        puts sr6.r, sr6.z, sr6.c, sr6.s
=end
  end
end

TestRotg.new


