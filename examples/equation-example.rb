require 'ratlas'
include RAtlas
require 'complex' 
include Math

class Blas
  def x
    self[0,0]
  end
  def y
    self[0,1]
  end
  def z
    self[0,2]
  end
end
=begin
p = Blas.new( 2,2, Blas::Double_t)
print p, "\n"

z = SingleBlas.new(1,2)
print z, "\n"

b = Blas.columns(1, Blas::Double_t, [3.6,4.5],[1.2,1.3])
print b, "\n"
b = Blas.rows(1, Blas::Single_t, [3.6,4.5],[1.2,1.3])
print b, "\n"
=end
c = SingleBlas.rows([3.6,4.5])
print c, "\n"
c = SingleBlas.columns([3.6,4.5],[1.2,1.3])
print c, "\n"
c = ComplexBlas.rows([3.6,4.5],[1.2,1.3])
print c, "\n"
c = ComplexBlas[[Complex(3.6,2.8), Complex(4.5,-8.2)],[Complex(1.2,4.4), Complex(1.3,9.1)]]
print c, "\n"


a = DoubleBlas[[2.0,3.0],[4.0,1.0]]
p = Blas.new( 2,2, Blas::Integer_t)

print a, "\n"
print p, "\n"

=begin

#a.xgesv(b,p)
#Gives result Vector[0.99, 0.54]

#print a, "\n"
#print b, "\n"
#print p, "\n"

print b.x, " ", b.y, "\n\n"

v1 = Blas[Blas::Float_t, 1.0, 4.5]
v2 = Blas[Blas::Float_t, 1.0, 4.6]
v3 = v2 - v1
print v1, "\n", v2, "\n", v3, "\n\n"

v4 = v1.scal(10.0)
print v4, "\n", v1, "\n\n"

v4.scal!(10.0)
print v4, "\n\n"

puts v1.dot(v2)
puts (v1 * v1)
puts v1 * 10.0 

puts Blas[Blas::Double_t, 1.0, 4.5, 2.3, 4.6, 0.9].max

v5 = Blas.new(10)
puts v5
puts v5.set!(2.0, 5, 2)

puts v1, v2
v1.swap!(v2)
puts v1, v2

puts v1.nrm2
=end