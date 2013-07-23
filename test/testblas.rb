require 'ratlas'
include RAtlas
require 'complex' 
include Math

class Numeric
  def within_bound(value, bound)
    return (self - value).abs <= bound
  end
end

class Complex
  def within_bound(value, bound)
    return (self - value).abs <= bound
  end
end


class TestBlas
  def initialize
    @flteps = 1e-4
    @dbleps = 1e-6
  end
  
  def print_on_error(message, result, expected, bound )
    #Small issue here, for matrix tests Blas subtraction, >, and abs must work.
    if result.within_bound(expected, bound) == false
      print "#{message} Got #{result} Expected #{expected}\n"
    end
  end
end