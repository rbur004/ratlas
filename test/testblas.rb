require 'ratlas'
include RAtlas
require 'complex' 
include Math

class TestBlas
  def initialize
    @flteps = 1e-4
    @dbleps = 1e-6
  end
  
  def print_on_error(message, result, expected, bound )
    #Small issue here, for matrix tests Blas subtraction, >, and abs must work.
    if (error = (result - expected)).abs > bound
      print "#{message} Error = #{error}\n"
    end
  end
end