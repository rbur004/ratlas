require 'ratlas'
include RAtlas
require '../testblas.rb'
require 'complex' 
include Math

class TestAsum < TestBlas

	def initialize 
	  super
]    ibm_examples
  end
  
  def test_mv(x,incx, expected, error_bound, test_message)
    result = x.amax(incx);
    #puts result
    print_on_error( "mv #{test_message}", result, expected, error_bound) 
  end
  
  def ibm_examples
    
  end

end

