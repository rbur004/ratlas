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
    a = DoubleLapack[ [3,2,-1], [6,6,2], [3,-2,1] ]
    p = IntegerBlas.new(3)
    puts a
    puts p
    puts "********xgetrf!************"
    a.xgetrf!(p)
    puts a
    puts p
    
    puts "********xgetrs!************"
    b = DoubleLapack[1,12,11]
    r = a.xgetrs!(b, p, Blas::NoTrans)
    puts a
    puts b
    puts r
    
    puts "********xgetrs************"
    b = DoubleLapack[1,12,11]
    r = a.xgetrs(b, p, Blas::NoTrans)
    puts a
    puts b
    puts r
  end
end

TestAxpy.new