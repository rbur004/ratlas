require 'ratlas'
include RAtlas
require 'complex' 
include Math
require 'testblas.rb'

class TestAxpy < TestBlas

	def initialize 
	  super
	  solve_these
  end
  
  def solve_these
    a = DoubleLapack[ [1,2,3], [0,1,4], [1,2,1] ]
    puts a
    print "Determinent #{a.determinant} == -2?\n\n"
    
    a = DoubleLapack[ [1,2,3], [0,1,4], [1,2,1] ]
    p = IntegerBlas.new(3)
    puts a
    puts p
    puts "********xgetrf!************"
    a.xgetrf!(p)
    puts a
    puts p
    
    puts "determinant of LU matrix #{a.determinant}\n"
    
    a = DoubleLapack[ [1,2,3], [1,2,3], [1,2,1] ]
    puts a
    print "Determinent of singular matrix #{a.determinant} == 0?\n\n"
    
  end
end

TestAxpy.new