require 'ratlas'
include RAtlas
require 'complex' 
include Math

a = SingleBlas[ 1.0, 2.0, 7.0, -8.0, -5.0, -10.0, -9.0, 10.0, 6.0 ]
b = SingleBlas[ 1.0, 2.0, 7.0, -8.0, -5.0, -10.0, -9.0, 10.0, 6.0 ]

puts "a[#{a.ncols},#{a.nrows}] = ", a, a.class
puts "b[#{b.ncols},#{b.nrows}] = ", b, b.class
puts
puts "(a < b) == #{a < b }"
puts "(a <= b) == #{a <= b }"


