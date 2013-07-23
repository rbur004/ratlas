require 'mkmf'
puts __dir__
find_header('cblas.h',"/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers") #additional include-path, also checks .h file exists
dir_config("/usr/local") #additional include-path and lib-path
dir_config("/usr/local/atlas`") #additional include-path and lib-path
dir_config("blas","#{__dir__}/blas/include", "#{__dir__}/blas/lib")
dir_config("lapack","#{__dir__}/lapack/include", "#{__dir__}/lapack/lib")
dir_config("rotg","#{__dir__}/rotg/include", "#{__dir__}/rotg/lib") #additional include-path and lib-path
#dir_config('ratlas')

#libraries given here link in reverse order to listing.
#Link order should be -llapack -lcblas -lf77blas -latlas
have_library('atlas')
#have_library('atlasf77blas')
have_library('blas')
have_library('cblas')
have_library('f77lapack')
have_library('lapack')
have_library('clapack')

File.open("conf.mk", "w") do |f_out|
  f_out.puts configuration('ratlas')
  f_out.puts "LIBPATH = #{libpathflag($LIBPATH)}"
end
#create_makefile('ratlas')

  

