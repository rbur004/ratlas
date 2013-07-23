require 'ratlas'
include RAtlas
require 'complex' 
include Math
require_relative '../testblas.rb'

class TestGesv < TestBlas

	def initialize 
	  super
	  ibm_examples
	  #solve_these
  end
  
  def test_gesv(a, bx, ipiv, ax_expected, bx_expected, ipiv_expected, error_bound, test_message)
    a.xgesv!(bx,ipiv)
    #puts a, ax_expected
    #puts
    #puts bx, bx_expected
    #puts
    #puts ipiv, ipiv_expected
    print_on_error( "xgesv a #{test_message}", a, ax_expected, error_bound) 
    print_on_error( "xgesv bx #{test_message}", bx, bx_expected, error_bound) 
    print_on_error( "xgesv ipiv #{test_message}", ipiv, ipiv_expected, error_bound) 
  end

  def solve_these
    a = DoubleLapack[ [3,2,-1], [6,6,2], [3,-2,1], ]
    b = DoubleLapack[1,12,11]
    p = IntegerBlas.new(3)
    puts a
    puts b
    puts p
    puts "********xgesv!************"
    a.xgesv!(b,p)
    puts a
    puts b
    puts p
    puts "********xgesv************"
    a = DoubleLapack[ [3,2,-1], [6,6,2], [3,-2,1], ]
    b = DoubleLapack[1,12,11]
    p = IntegerBlas.new(3)
    r = a.xgesv(b,p)
    puts a
    puts b
    puts r
    puts p
  end
  
  def ibm_examples
    #Example 1
    #    This example shows how to solve the system AX = equal to B, where:

    #    Matrix A is the same used as input in Example 1 for DGETRF.
    #    Matrix B is the same used as input in Example 1 for DGETRS.
    #    Call Statement and Input:

    #                N  NRHS  A  LDA  IPVT  BX  LDB  INFO
    #                |   |    |   |    |    |    |    | 
    #    CALL DGESV( 9 , 5 ,  A , 9 , IPVT, BX , 9 , INFO)

    #            | 1.0  1.2  1.4  1.6  1.8  2.0  2.2  2.4  2.6 |
    #            | 1.2  1.0  1.2  1.4  1.6  1.8  2.0  2.2  2.4 |
    #            | 1.4  1.2  1.0  1.2  1.4  1.6  1.8  2.0  2.2 |
    #            | 1.6  1.4  1.2  1.0  1.2  1.4  1.6  1.8  2.0 |
    #    A    =  | 1.8  1.6  1.4  1.2  1.0  1.2  1.4  1.6  1.8 |
    #            | 2.0  1.8  1.6  1.4  1.2  1.0  1.2  1.4  1.6 |
    #            | 2.2  2.0  1.8  1.6  1.4  1.2  1.0  1.2  1.4 |
    #            | 2.4  2.2  2.0  1.8  1.6  1.4  1.2  1.0  1.2 |
    #            | 2.6  2.4  2.2  2.0  1.8  1.6  1.4  1.2  1.0 |
    #     
    #             | 93.0  186.0  279.0  372.0  465.0 |
    #             | 84.4  168.8  253.2  337.6  422.0 |
    #             | 76.6  153.2  229.8  306.4  383.0 |
    #             | 70.0  140.0  210.0  280.0  350.0 |
    #     BX   =  | 65.0  130.0  195.0  260.0  325.0 |
    #             | 62.0  124.0  186.0  248.0  310.0 |
    #             | 61.4  122.8  184.2  245.6  307.0 |
    #             | 63.6  127.2  190.8  254.4  318.0 |
    #             | 69.0  138.0  207.0  276.0  345.0 |

    #    Output:

    #    IPIV     =  (9, 9, 9, 9, 9, 9, 9, 9, 9) 

    #            | 2.6   2.4  2.2  2.0  1.8  1.6  1.4  1.2  1.0 |
    #            | 0.4   0.3  0.6  0.8  1.1  1.4  1.7  1.9  2.2 |
    #            | 0.5  -0.4  0.4  0.8  1.2  1.6  2.0  2.4  2.8 |
    #            | 0.5  -0.3  0.0  0.4  0.8  1.2  1.6  2.0  2.4 |
    #    A    =  | 0.6  -0.3  0.0  0.0  0.4  0.8  1.2  1.6  2.0 |
    #            | 0.7  -0.2  0.0  0.0  0.0  0.4  0.8  1.2  1.6 |
    #            | 0.8  -0.2  0.0  0.0  0.0  0.0  0.4  0.8  1.2 |
    #            | 0.8  -0.1  0.0  0.0  0.0  0.0  0.0  0.4  0.8 |
    #            | 0.9  -0.1  0.0  0.0  0.0  0.0  0.0  0.0  0.4 |

    #            | 1.0   2.0   3.0   4.0   5.0 |
    #            | 2.0   4.0   6.0   8.0  10.0 |
    #            | 3.0   6.0   9.0  12.0  15.0 |
    #            | 4.0   8.0  12.0  16.0  20.0 |
    #    BX   =  | 5.0  10.0  15.0  20.0  25.0 |
    #            | 6.0  12.0  18.0  24.0  30.0 |
    #            | 7.0  14.0  21.0  28.0  35.0 |
    #            | 8.0  16.0  24.0  32.0  40.0 |
    #            | 9.0  18.0  27.0  36.0  45.0 |

    #    INFO  =  0
    a = DoubleLapack[
            [ 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6 ],
            [ 1.2, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4 ],
            [ 1.4, 1.2, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2 ],
            [ 1.6, 1.4, 1.2, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 ],
            [ 1.8, 1.6, 1.4, 1.2, 1.0, 1.2, 1.4, 1.6, 1.8 ],
            [ 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 1.2, 1.4, 1.6 ],
            [ 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 1.2, 1.4 ],
            [ 2.4, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 1.2 ],
            [ 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0 ],
          ]
            
    bx = DoubleLapack[
              [ 93.0, 186.0, 279.0, 372.0, 465.0 ],
              [ 84.4, 168.8, 253.2, 337.6, 422.0 ],
              [ 76.6, 153.2, 229.8, 306.4, 383.0 ],
              [ 70.0, 140.0, 210.0, 280.0, 350.0 ],
              [ 65.0, 130.0, 195.0, 260.0, 325.0 ],
              [ 62.0, 124.0, 186.0, 248.0, 310.0 ],
              [ 61.4, 122.8, 184.2, 245.6, 307.0 ],
              [ 63.6, 127.2, 190.8, 254.4, 318.0 ],
              [ 69.0, 138.0, 207.0, 276.0, 345.0 ],
            ]
                  
    a_expected = DoubleLapack[
                [ 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0],
                [ 0.4, 0.3, 0.6, 0.8, 1.1, 1.4, 1.7, 1.9, 2.2 ],
                [ 0.5, -0.4, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4, 2.8],
                [ 0.5, -0.3, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4 ],
                [ 0.6, -0.3, 0.0, 0.0, 0.4, 0.8, 1.2, 1.6, 2.0],
                [ 0.7, -0.2, 0.0, 0.0, 0.0, 0.4, 0.8, 1.2, 1.6 ],
                [ 0.8, -0.2, 0.0, 0.0, 0.0, 0.0, 0.4, 0.8, 1.2 ],
                [ 0.8, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.8],
                [ 0.9, -0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.4 ],
              ]

    bx_expected =  DoubleLapack[
                  [ 1.0, 2.0, 3.0, 4.0, 5.0 ],
                  [ 2.0, 4.0, 6.0, 8.0, 10.0 ],
                  [ 3.0, 6.0, 9.0, 12.0, 15.0 ],
                  [ 4.0, 8.0, 12.0, 16.0, 20.0 ],
                  [ 5.0, 10.0, 15.0, 20.0, 25.0 ],
                  [ 6.0, 12.0, 18.0, 24.0, 30.0 ],
                  [ 7.0, 14.0, 21.0, 28.0, 35.0 ],
                  [ 8.0, 16.0, 24.0, 32.0, 40.0 ],
                  [ 9.0, 18.0, 27.0, 36.0, 45.0 ],
                ]

    ipiv = IntegerBlas.new(9)
    ipiv_expected = IntegerBlas[9, 9, 9, 9, 9, 9, 9, 9, 9]
    
    test_gesv(a, bx, ipiv, a_expected, bx_expected, ipiv_expected, 0.1, "IBM Example gesv-01")
    #Example 2
    #    This example shows how to solve the system AX = equal to B, where:

    #    Matrix A is the same used as input in Example 2 for ZGETRF.
    #    Matrix B is the same used as input in Example 2 for ZGETRS.
    #    Call Statement and Input:

    #                N  NRHS  A  LDA  IPVT  BX  LDB  INFO
    #                |   |    |   |    |    |    |    | 
    #    CALL ZGESV( 9 , 5 ,  A , 9 , IPVT, BX,  9 , INFO)

    #    A      = (same as input A in Example 2)
    #    IPVT   = (same as input IPVT in Example 2)
    #    BX     = (same as input BX in Example 2)
    #    Output:

    #            ┌                                                     ┐
    #            | (1.0,1.0)  (1.0,2.0)  (1.0,3.0) (1.0,4.0) (1.0,5.0) |
    #            | (2.0,1.0)  (2.0,2.0)  (2.0,3.0) (2.0,4.0) (2.0,5.0) |
    #            | (3.0,1.0)  (3.0,2.0)  (3.0,3.0) (3.0,4.0) (3.0,5.0) |
    #            | (4.0,1.0)  (4.0,2.0)  (4.0,3.0) (4.0,4.0) (4.0,5.0) |
    #    BX    = | (5.0,1.0)  (5.0,2.0)  (5.0,3.0) (5.0,4.0) (5.0,5.0) |
    #            | (6.0,1.0)  (6.0,2.0)  (6.0,3.0) (6.0,4.0) (6.0,5.0) |
    #            | (7.0,1.0)  (7.0,2.0)  (7.0,3.0) (7.0,4.0) (7.0,5.0) |
    #            | (8.0,1.0)  (8.0,2.0)  (8.0,3.0) (8.0,4.0) (8.0,5.0) |
    #            | (9.0,1.0)  (9.0,2.0)  (9.0,3.0) (9.0,4.0) (9.0,5.0) |
    #            └                                                     ┘

    #    INFO     =  0
    
  end
end

TestGesv.new