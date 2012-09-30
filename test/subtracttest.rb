require 'ratlas'
include RAtlas
require 'complex' 
include Math

class Blas
  def pp
    s = "[\n"
    self.each_row do |r,i|
      s += "  ["
      r.each { |c| s += "#{sprintf("%2.2f", c)}, "}
      s[-1] = "]"
      s += "\n"
    end
    s += "]"
  end
end


a = SingleBlas[10,12,16]
b = SingleBlas[5,6,3]
puts "a = #{a.pp}"
puts "b = #{b.pp}\nrows = #{b.nrows}\ncol = #{b.ncols}"
puts "a+b = #{(a + b).pp}"
puts "a.axpy(b) = #{a.axpy(b).pp}"
puts 
puts "a = #{a.pp}"
puts "b = #{b.pp}"
puts "a-b = #{(a - b).pp}"
puts 

puts "a = #{a.pp}"
puts "b = #{b.pp}"
puts "a*b = #{(a * b)}"
puts 

a = DoubleBlas[10,12,16]
b = DoubleBlas[5,6,3]
puts "a = #{a.pp}"
puts "b = #{b.pp}"
puts "a+b = #{(a + b).pp}"
puts 
puts "a = #{a.pp}"
puts "b = #{b.pp}"
puts "a-b = #{(a - b).pp}"
puts 

puts "a = #{a.pp}"
puts "b = #{b.pp}"
puts "a*b = #{(a * b)}"
puts 

a = DoubleLapack[10,12,16]
b = DoubleLapack[5,6,3]
puts "a = #{a.pp}"
puts "b = #{b.pp}"
puts "a+b = #{(a + b).pp}"
puts 
puts "a = #{a.pp}"
puts "b = #{b.pp}"
puts "a-b = #{(a - b).pp}"
puts 

puts "a = #{a.pp}"
puts "b = #{b.pp}"
puts "a*b = #{(a * b)}"
puts 
