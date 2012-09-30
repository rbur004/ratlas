require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestConstructors < TestBlas
  def initialize
    super
    test_constructors
  end
  
  def validate(matrix, array, error_bound, testmsg)
    matrix.each_by_row_with_index do |v, r,c|
     # print "#{testmsg}: v #{v}, r #{r}, c #{c}, array[r][c]  #{array[r][c]}, error_bound #{error_bound}\n"
      print_on_error( "Test Constructors #{testmsg} by row: matrix[#{r},#{c}] #{v} != array[#{r},#{c}] #{array[r][c]}\n" , v,  array[r][c], error_bound)
    end
    
    matrix.each_by_col_with_index do |v, r,c|
     # print "#{testmsg}: v #{v}, r #{r}, c #{c}, array[r][c]  #{array[r][c]}, error_bound #{error_bound}\n"
      print_on_error( "Test Constructors #{testmsg} by column: matrix[#{r},#{c}] #{v} != array[#{r},#{c}] #{array[r][c]}\n" , v,  array[r][c], error_bound)
    end
    
    #Reverse the test to see if we can index the Matrix members correctly.
    array.each_with_index do |row, r|
      row.each_with_index do |v, c|
        print_on_error( "Test Constructors #{testmsg}:  array[#{r},#{c}] #{v} != matrix.get(#{r},#{c}) #{matrix.get(r,c)}\n" , v,   matrix.get(r,c), error_bound)
      end
    end
    
    #Reverse the test to see if we can index the Matrix members correctly using [] notation.
    array.each_with_index do |row, r|
      row.each_with_index do |v, c|
        print_on_error(  "Test Constructors #{testmsg}:  array[#{r},#{c}] #{v} != matrix[#{r},#{c}] #{matrix[r,c]}\n" , v,   matrix[r,c], error_bound)
      end
    end
  end
  
  def test_constructors
    validate(Blas.rows(Blas::GE, Blas::S, [3.6,4.5], [1.2,1.3]), [[3.6,4.5], [1.2,1.3]], @flteps, "Test 1")
    validate(SingleBlas.rows([3.6,4.5],[1.2,1.3]), [[3.6,4.5],[1.2,1.3]], @flteps, "Test 2")
    validate(SingleBlas.columns([3.6,4.5],[1.2,1.3]), [[3.6,1.2],[4.5,1.3]], @flteps, "Test 3")
    validate(SingleBlas[ [3.6,4.5],[1.2,1.3] ], [ [3.6,4.5],[1.2,1.3] ], @flteps, "Test 4")

    validate(Blas.rows(Blas::GE, Blas::D, [3.6,4.5],[1.2,1.3]), [[3.6,4.5],[1.2,1.3]], @dbleps, "Test 5")
    validate(DoubleBlas.rows([3.6,4.5],[1.2,1.3]), [[3.6,4.5],[1.2,1.3]], @dbleps, "Test 6")
    validate(DoubleBlas.columns([3.6,4.5],[1.2,1.3]), [[3.6,1.2],[4.5,1.3]], @dbleps, "Test 7")
    validate(DoubleBlas[ [3.6,4.5],[1.2,1.3] ], [ [3.6,4.5],[1.2,1.3] ], @dbleps, "Test 8")

    validate(Blas.rows(Blas::GE, Blas::C, [Complex(3.6,2.8), Complex(4.5,-8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]), [[Complex(3.6,2.8), Complex(4.5,-8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]], Complex(@flteps, @flteps), "Test 9")
    validate(ComplexBlas.rows([Complex(3.6,2.8), Complex(4.5,-8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]), [[Complex(3.6,2.8), Complex(4.5,-8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]], Complex(@flteps, @flteps), "Test 10")
    validate(ComplexBlas.columns([Complex(3.6,2.8), Complex(4.5,-8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]), [[Complex(3.6,2.8), Complex(1.2,4.4)],[Complex(4.5,-8.2), Complex(1.3,9.1)]], Complex(@flteps, @flteps), "Test 11")
    validate(ComplexBlas[ [Complex(3.6,2.8), Complex(4.5,-8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)] ], [ [Complex(3.6,2.8), Complex(4.5,-8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)] ], Complex(@flteps, @flteps), "Test 12")
 
    validate(Blas.rows(Blas::GE, Blas::Z, [Complex(3.6,2.8), Complex(4.5,8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]), [[Complex(3.6,2.8), Complex(4.5,8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]], Complex(@dbleps, @dbleps), "Test 13")
    validate(DoubleComplexBlas.rows([Complex(3.6,2.8), Complex(4.5,8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]), [[Complex(3.6,2.8), Complex(4.5,8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]], Complex(@dbleps, @dbleps), "Test 14")
    validate(DoubleComplexBlas.columns([Complex(3.6,2.8), Complex(4.5,8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]), [[Complex(3.6,2.8), Complex(1.2,4.4) ],[Complex(4.5,8.2), Complex(1.3,9.1)]], Complex(@dbleps, @dbleps), "Test 15")
    validate(DoubleComplexBlas[ [Complex(3.6,2.8), Complex(4.5,8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)] ], [ [Complex(3.6,2.8), Complex(4.5,8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)] ], Complex(@dbleps, @dbleps), "Test 16")

    #Test the complex constructors being passed real values, with no complex component.
    validate(Blas.rows(Blas::GE, Blas::C, [3.6,4.5], [1.2,1.3]), [[3.6,4.5], [1.2,1.3]], @flteps, "Test 17")
    validate(ComplexBlas.rows([3.6,4.5],[1.2,1.3]), [[3.6,4.5],[1.2,1.3]], @flteps, "Test 18")
    validate(ComplexBlas.columns([3.6,4.5],[1.2,1.3]), [[3.6,1.2],[4.5,1.3]], @flteps, "Test 19")
    validate(ComplexBlas[ [3.6,4.5],[1.2,1.3] ], [ [3.6,4.5],[1.2,1.3] ], @flteps, "Test 20")

    validate(Blas.rows(Blas::GE, Blas::Z, [3.6,4.5],[1.2,1.3]), [[3.6,4.5],[1.2,1.3]], @dbleps, "Test 21")
    validate(DoubleComplexBlas.rows([3.6,4.5],[1.2,1.3]), [[3.6,4.5],[1.2,1.3]], @dbleps, "Test 22")
    validate(DoubleComplexBlas.columns([3.6,4.5],[1.2,1.3]), [[3.6,1.2],[4.5,1.3]], @dbleps, "Test 23")
    validate(DoubleComplexBlas[ [3.6,4.5],[1.2,1.3] ], [ [3.6,4.5],[1.2,1.3] ], @dbleps, "Test 24")

    validate(Blas.rows(Blas::GE, Blas::I, [3,4],[1,2]), [[3,4],[1,2]], 0, "Test 25")
    validate(IntegerBlas.rows([3,4],[1,2]), [[3,4],[1,2]], 0, "Test 26")
    validate(IntegerBlas.columns([3,4],[1,2]), [[3,1],[4,2]], 0, "Test 27")
    validate(IntegerBlas[ [3,4],[1,2] ], [ [3,4],[1,2] ], 0,"Test 28")

  end
end

TestConstructors.new