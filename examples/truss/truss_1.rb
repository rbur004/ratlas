require "../truss"

print "Tensions are positive, Compression negative\n"
print "0 point is floating end, Forces have opposite sign to members.\n"
print "Forces at supports can be added as Ay (floating end), Dx, Dy (fixed end). Set mount_forces to true\n\n"

print "Structural & Stress Analysis. 2nd Ed. Megson pg 109 Prob 4.1\n\n"
truss = Truss.new(
  [
      ['A', [0,0], [0,30]],
      ['B', [6,0], [0,0]],
      ['C', [12,0], [0,0]],
      ['D', [18, 0], [0,0]],
      ['E', [24, 0], [0,0]],
      ['F', [30, 0], [0,60]],
      ['G', [6, 8], [0,0]],
      ['H', [12, 8], [0,0]],
      ['J', [18, 8], [0,0]],
      ['K', [24, 8], [0,0]]
  ],
  [
    ['AG', ['A', 'G']],
    ['AB', ['A', 'B']],
    ['BG', ['B','G']],
    ['BC', ['B','C']],
    ['GC', ['G','C']],
    ['GH', ['G','H']],
    ['HC', ['H','C']],
    ['HJ', ['H','J']],
    ['CJ', ['C','J']],
    ['CD', ['C','D']],
    ['JD', ['J', 'D']],
    ['JK', ['J', 'K']],
    ['DK', ['D', 'K']],
    ['DE', ['D', 'E']],
    ['EK', ['E', 'K']],    
    ['FE', ['F', 'E']],    
    ['FK', ['F', 'K']],
    #supports
    ['By', ['B', 'Fy']],
    ['Bx', ['B', 'Fx']],
    ['Ey', ['E', 'Fy']]
  ],true
)
v = truss.solve
print "Truss = ", truss.joint_to_s, "\n"
#print "Solution Vector = \n", truss.solution_to_s, "\n"
#print "Test = " , truss.test_solution, "\n"
print "Forces\n", truss.forces_to_s, "\n"


=begin
print "\n\n3D test of truss code. Value from book example. Also wanting mount point forces\n"
truss = Truss.new(
             [['A', [-20.0, 8.66, -5.0], [0.0, 0.0, 0.0]],
              ['B', [-10.0, 6.93, -4.0], [0.0, 0.0, 0.0]],
              ['C', [  0.0, 5.2,  -3.0], [0.0, 6.0, 0.0]],
              ['D', [-20.0, 0.0,   0.0], [0.0, 0.0, 0.0]],
              ['E', [-10.0, 0.0,   0.0], [0.0, 0.0, 0.0]],
              ['F', [  0.0, 0.0,   0.0], [0.0, 10.0, 0.0]],
              ['G', [-20.0, 8.66,  5.0], [0.0, 0.0, 0.0]],
              ['H', [-10.0, 6.93,  4.0], [0.0, 0.0, 0.0]],
              ['J',  [ 0.0, 5.2,   3.0], [0.0, 3.0, 0.0]]
             ],
             [['FC', ['F', 'C']], 
              ['FJ', ['F', 'J']], 
              ['FB', ['F', 'B']],
              ['FE', ['F', 'E']], 
              ['FH', ['F', 'H']],
              ['CB', ['C', 'B']], 
              ['CJ', ['C', 'J']], 
              ['CH', ['C', 'H']], 
              ['JH', ['J', 'H']],
              ['EB', ['E', 'B']], 
              ['ED', ['E', 'D']], 
              ['EA', ['E', 'A']], 
              ['EG', ['E', 'G']], 
              ['EH', ['E', 'H']],
              ['BA', ['B', 'A']], 
              ['BG', ['B', 'G']], 
              ['BH', ['B', 'H']], 
              ['HG', ['H', 'G']],
              #forces at the support joints
              ['Ax', ['A', 'Fx']],
              ['Ay', ['A', 'Fy']],
              ['Az', ['A', 'Fz']],
              ['Gx', ['G', 'Fx']],
              ['Gy', ['G', 'Fy']],
              ['Gz', ['G', 'Fz']],
              ['Dx', ['D', 'Fx']],
              ['Dy', ['D', 'Fy']],
              ['Dz', ['D', 'Fz']]
             ], false
            )
v = truss.solve
#print "Matrix = \n", truss.matrix_to_s, "\n"
print "Force Vector = ", truss.force_vector_to_s, "\n"
print "Test = " , truss.test_solution, "\n"

print "\n\nTruss from text, with Bx, By and Cy already given\n"
truss2 = Truss.new(
          [ ['A', [-15.0, 8.66], [0.5, 0.0] ],
            ['B', [-20.0, 0.0], [-1.0, -2.2] ],
            ['C', [ 0.0, 0.0], [0.0, -0.8] ],
            ['D', [ -7.5, 4.33], [0.5, 0.0] ],
            ['E', [ -15.0, 0.0], [0.0, 2.0] ],
            ['F', [ -7.5, 0.0], [0.0, 1.0] ]
          ],
          [ ['CD', ['C', 'D']],
            ['CF', ['C', 'F']],
            ['FD', ['F', 'D']],
            ['FE', ['F', 'E']],
            ['BA', ['B', 'A']],
            ['BE', ['B', 'E']],
            ['DE', ['D', 'E']],
            ['DA', ['D', 'A']],
            ['AE', ['A', 'E']]
          ], false
        )
v = truss2.solve
print "Force Vector = ", truss2.force_vector_to_s, "\n"
print "Test = " , truss2.test_solution, "\n"

=end
=begin

print "\n\nTruss from text, with Bx, By and Cy NOT given\n"
truss2a = Truss.new(
          [ ['A', [-15.0, 8.66], [0.5, 0.0] ],
            ['B', [-20.0, 0.0], [0.0, 0.0] ],
            ['C', [ 0.0, 0.0], [0.0, 0.0] ],
            ['D', [ -7.5, 4.33], [0.5, 0.0] ],
            ['E', [ -15.0, 0.0], [0.0, 2.0] ],
            ['F', [ -7.5, 0.0], [0.0, 1.0] ]
          ],
          [ ['CD', ['C', 'D'] ],
            ['CF', ['C', 'F'] ],
            ['FD', ['F', 'D'] ],
            ['FE', ['F', 'E'] ],
            ['BA', ['B', 'A'] ],
            ['BE', ['B', 'E'] ],
            ['DE', ['D', 'E'] ],
            ['DA', ['D', 'A'] ],
            ['AE', ['A', 'E'] ],
            ['Bx', ['B', 'Fx'] ],
            ['By', ['B', 'Fy'] ],
            ['Cy', ['C', 'Fy'] ]
          ]
        )
v = truss2a.solve
print "Force Vector = ", truss2a.force_vector_to_s, "\n"
print "Test = " , truss2a.test_solution, "\n"
exit
=begin

class Numeric
  def degrees2radians
    #converts degrees to radians
    self*(Math::PI*2)/360
  end
end

def calc(truss_angle)
print "Calculate forces if truss @ angle #{truss_angle}\n"
y1 = 1.2*Math::sin((60+truss_angle).degrees2radians)
x1 = 1.2*Math::cos((60+truss_angle).degrees2radians)
y2 = 1.2*Math::sin(truss_angle.degrees2radians)
x2 = 1.2*Math::cos(truss_angle.degrees2radians)
#print "(#{x1}, #{y1}), (#{x2} , #{y2})\n"

truss3 = Truss.new(
                      #joint    Coordinates           x & y forces at joint
                    [ ['A', [ -(x1+x2*2), y1+y2*2],  [0.0, 600.0] ],
                      ['B', [ -x2*3, y2*3],          [0.0, 0.0] ],
                      ['C', [ 0.0, 0.0],             [0.0, 0.0] ],
                      ['D', [ -x1, y1],              [0.0, 600.0] ],
                      ['E', [ -x2, y2],              [0.0, 0.0] ],
                      ['F', [ -(x1+x2), y1+y2],      [0.0, 600.0] ],
                      ['G', [ -x2*2, y2*2],          [0.0, 0.0] ]
                    ],
                    # member-name   start & end joints
                    [ ['CD',        ['C', 'D']],
                      ['CE',        ['C', 'E']],
                      ['ED',        ['E', 'D']],
                      ['EF',        ['E', 'F']],
                      ['EG',        ['E', 'G']],
                      ['DF',        ['D', 'F']],
                      ['FA',        ['F', 'A']],
                      ['GF',        ['G', 'F']],
                      ['GA',        ['G', 'A']],
                      ['GB',        ['G', 'B']],
                      ['BA',        ['B', 'A']],
                      #Next three are the forces at the joints supporting the truss.
                      ['BX', ['B', 'Fx'] ], #dummy value 'Fx' means unknown force in x plane
                      ['BY', ['B', 'Fy'] ], #dummy value 'Fy' means unknown force in y plane
                      ['CY', ['C', 'Fy'] ]  #dummy value 'Fy' means unknown force in y plane
                    ]
                  )
v = truss3.solve
print "Force Vector = ", truss3.force_vector_to_s, "\n"
print "Test = " , truss3.test_solution, "\n"
end

"Made up 2D truss, to see the variation of forces at different angles\n"
calc(0)
calc(19)
calc(45)
=end


