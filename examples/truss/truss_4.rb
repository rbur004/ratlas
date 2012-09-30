require "../truss"

print "Tensions are positive, Compression negative\n"
print "0 point is floating end, Forces have opposite sign to members.\n"
print "Forces at supports can be added as Ay (floating end), Dx, Dy (fixed end). Set mount_forces to true\n\n"

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
            ['AE', ['A', 'E']],
            ['Bx', ['B', 'Fx'] ],
            ['By', ['B', 'Fy'] ],
            ['Cy', ['C', 'Fy'] ]
          ]
        )
v = truss2.solve
print "Truss from page 16 The analysis of Engineering Structures. Pippard & Barker Third edition\n"
print "Should give:\n Tensions CD -1.6, CF 1.39, FD 1, FE 1.39, BA -2.54, BE 0.27, DE -1.29, DA -0.89, AE 2.64\n"
print "Bx, By and Cy already given, but also asked for\n"
print "Hence the Bx, By, and Cy value will be 0 (or close to 0)\n"

print "Truss = ", truss2.joint_to_s, "\n"
#print "Solution Vector = \n", truss.solution_to_s, "\n"
#print "Test = " , truss.test_solution, "\n"
print "Forces\n", truss2.forces_to_s, "\n"
