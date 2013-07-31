require 'ratlas'
include RAtlas
require 'complex' 
include Math

a = ComplexBlas[ Complex(1.0, 2.0), Complex(7.0, -8.0), Complex(-5.0, -10.0), Complex(-9.0, 10.0) ]
puts a.class
#puts a.methods
puts "data type #{a.data_type}, data size #{a.data_size}"
puts "a.nrows #{a.nrows}, a.ncols #{a.ncols}"
a.each_row do |c| 
  puts c.length
end

#puts a.to_s


#a.each do |c|
#  puts c.class
#end
