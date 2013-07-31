require 'ratlas'
include RAtlas
require 'complex' 
include Math
class Numeric
  def within_bound(value, bound)
    puts "in Numeric"
    puts self.class, value.class, bound.class
    return (self - value).abs <= bound
  end
end

a = SingleBlas[ 1.0, 2.0, 7.0, -8.0, -5.0, -10.0, -9.0, 10.0, 6.0 ]
b = SingleBlas[ 1.0, 2.0, 7.0, -7.0, -5.0, -10.0, -9.0, 10.0, 6.0 ]
bound = 1.0

puts "a[#{a.ncols},#{a.nrows}] = ", a, a.class
puts "b[#{b.ncols},#{b.nrows}] = ", b, b.class
puts "a-b = ", a-b
puts
puts "(a - b) < #{bound} == #{a.within_bound(b,bound)}"
puts "(a - b) < #{bound/2} == #{a.within_bound(b,bound/2)}"

puts "Complex"
a = ComplexBlas[ Complex(1.0, 2.0), Complex(7.0, -8.0), Complex(-5.0, -10.0), Complex(-9.0, 10.0) ]
a.each do |c|
  puts c.class
end
b = ComplexBlas[ Complex(1.0, 2.0), Complex(7.0, -7.0), Complex(-5.0, -10.0), Complex(-9.0, 10.0) ]
bound = 1.0

puts "a[#{a.ncols},#{a.nrows}] = ", a, a.class
puts "b[#{b.ncols},#{b.nrows}] = ", b, b.class
puts "a-b = ", a-b
puts
puts "(a - b) < #{bound} == #{a.within_bound(b,bound)}"
puts "(a - b) < #{bound/2} == #{a.within_bound(b,bound/2)}"
