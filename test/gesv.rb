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
    a = DoubleLapack[ [3,2,-1], [6,6,2], [3,-2,1], ]
    b = DoubleLapack[1,12,11]
    p = IntegerBlas.new(3)
    puts a
    puts b
    puts p
    puts "********xgesv!************"
    a.xgesv!(b,p)
    puts a
    puts b
    puts p
    puts "********xgesv************"
    a = DoubleLapack[ [3,2,-1], [6,6,2], [3,-2,1], ]
    b = DoubleLapack[1,12,11]
    p = IntegerBlas.new(3)
    r = a.xgesv(b,p)
    puts a
    puts b
    puts r
    puts p
    
    
  end
end

TestAxpy.new