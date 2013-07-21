require 'ratlas'
include RAtlas

class Lapack
  def pp
    s = "[\n"
    self.each_row do |r,i|
      s += "  ["
      r.each { |c| s += "#{sprintf("%2.1f",c)}, "}
      s += "]\n"
    end
    s += "]"
  end
end

a = DoubleLapack[1,2,3,4,5]
puts a

a.each { | i | puts i }

a.each_by_row { | i | puts i.class }

a.each_row do |r,j|
  puts r.class, j
  r.each { |x| puts x.class }
end

puts a.pp

