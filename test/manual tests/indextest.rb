require 'ratlas'
include RAtlas

a = SingleBlas[1,2,3,4,5]
puts a[1]

a[2] = 8
puts a

a = SingleLapack[1,2,3,4,5]
puts a[1]
a[2] = 8
puts a
