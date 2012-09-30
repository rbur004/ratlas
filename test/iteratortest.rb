require 'ratlas'
include RAtlas
require 'complex' 
include Math

puts "********* each Vector *********"
a = IntegerBlas[12,10,8]
puts a
a.each { |i| puts i }

puts "********* each 3x3 Matrix ***********"
a = IntegerBlas[[10,20,30],[40,50,60],[70,80,90]]
puts a
a.each { |i| puts i }

puts "********* each_by_row Vector *********"
a = IntegerBlas[12,10,8]
puts a
a.each_by_row { |i| puts i }

puts "********* each_by_row 3x3 Matrix *********"
a = IntegerBlas[[1,2,3],[4,5,6],[7,8,9]]
puts a
a.each_by_row { |i| puts i }

puts "*********  each_by_row_with_index Vector *********"
a = IntegerBlas[12,10,8]
puts a
a.each_by_row_with_index { |i,r,c| puts "a[#{r},#{c}] = #{i}" }

puts "*********  each_by_row_with_index 3x3 Matrix *********"
a = IntegerBlas[[1,2,3],[4,5,6],[7,8,9]]
puts a
a.each_by_row_with_index { |i,r,c| puts "a[#{r},#{c}] = #{i}" }


puts "********* each_row Vector *********"
a = IntegerBlas[12,10,8]
puts a
a.each_row { |i,r| puts "row #{r} => #{i}"}

puts "********* each_row 3x3 Matrix *********"
a = IntegerBlas[[10,20,30],[40,50,60],[70,80,90]]
puts a
a.each_row { |i,r| puts "row #{r} => #{i}"}

puts "********* each_row 3x3 Matrix  .each *********"
a = IntegerBlas[[10,20,30],[40,50,60],[70,80,90]]
puts a
a.each_row do |i,r| 
  i.each { |j| puts "row #{r} => #{j}" }
end



puts "********* each_by_col 3x3 Vector *********"
a = IntegerBlas[12,10,8]
a.each_by_col { |i| puts i }

puts "********* each_by_col 3x3 Matrix *********"
a = IntegerBlas[[10,20,30],[40,50,60],[70,80,90]]
a.each_by_col { |i| puts i }

puts "*********  each_by_col_with_index Vector *********"
a = IntegerBlas[12,10,8]
a.each_by_col_with_index { |i,r,c| puts "a[#{r},#{c}] = #{i}" }

puts "*********  each_by_col_with_index 3x3 Matrix *********"
a = IntegerBlas[[1,2,3],[4,5,6],[7,8,9]]
a.each_by_col_with_index { |i,r,c| puts "a[#{r},#{c}] = #{i}" }

puts "********* each_col Vector *********"
a = IntegerBlas[12,10,8]
a.each_col { |i,c| puts "column #{c} => #{i}"}

puts "********* each_col 3x3 Matrix *********"
a = IntegerBlas[[10,20,30],[40,50,60],[70,80,90]]
a.each_col { |i,c| puts "column #{c} => #{i}"}

