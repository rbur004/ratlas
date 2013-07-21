require 'ratlas'
include RAtlas
require 'complex' 
include Math

class Numeric
  def within_bound(bound)
    return self.abs <= bound
  end
end

class Numeric
  def within_bound(bound)
    return self.abs <= bound
  end
end

class TestBlas
  def initialize
    @flteps = 1e-4
    @dbleps = 1e-6
  end
  
  def print_on_error(message, result, expected, bound )
    #Small issue here, for matrix tests Blas subtraction, >, and abs must work.
    if (error = (result - expected)).within_bound(bound) == false
      print "#{message} Got #{result} Expected #{expected} Error = #{error}\n"
    end
  end
end