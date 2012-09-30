require 'mkmf'
dir_config('ratlas', "/Developer/SDKs/MacOSX10.5.sdk/System/Library/Frameworks/vecLib.framework/Versions/A/Headers/", "/Developer/SDKs/MacOSX10.5.sdk//usr/lib/") #optional include-path, optional lib-path
#dir_config('ratlas', "/usr/local/include/") #optional include-path, optional lib-path
dir_config('ratlas')

#libraries given here link in reverse order to listing.
#Link order should be -llapack -lcblas -lf77blas -latlas
have_library('atlas')
#have_library('atlasf77blas')
have_library('blas')
have_library('cblas')
have_library('f77lapack')
have_library('lapack')
have_library('clapack')
create_makefile('ratlas')
