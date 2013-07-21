require 'ratlas'
include RAtlas
require 'complex' 
include Math

a = SingleBlas[ 1.0, 2.0, 7.0, -8.0, -5.0, -10.0, -9.0, 10.0, 6.0 ]
b = SingleBlas[ 1.0, 2.0, 7.0, -7.0, -5.0, -10.0, -9.0, 10.0, 6.0 ]
bound = 1.0

puts "a[#{a.ncols},#{a.nrows}] = ", a, a.class
puts "b[#{b.ncols},#{b.nrows}] = ", b, b.class
puts "a-b = ", a-b
puts
puts "(a - b) < #{bound} == #{(a - b).in_bound(bound)}"
puts "(a - b) < #{bound/2} == #{(a - b).in_bound(bound/2)}"
