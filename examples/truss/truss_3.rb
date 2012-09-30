require "../truss"

print "Tensions are positive, Compression negative\n"
print "0 point is floating end, Forces have opposite sign to members.\n"
print "Forces at supports can be added as Ay (floating end), Dx, Dy (fixed end). Set mount_forces to true\n\n"

print "Eng Mech Statics 4th ed. Meriam & Kraige problem 4/29 Pg 195 \n"
print "Problem expected BC Compression 4.12, BE Tension 0.901, EF Tension 3.38\n\n"
truss = Truss.new(
  [
      ['A', [0,0], [0,0]],
      ['B', [3,4], [0,0]],
      ['C', [9,4], [0,0]],
      ['D', [12, 0], [0,0]],
      ['E', [9, 0], [0,6]],
      ['F', [3, 0], [0,4]]
  ],
  [
    ['AB', ['A', 'B']],
    ['AF', ['A', 'F']],
    ['BF', ['B','F']],
    ['BC', ['B','C']],
    ['BE', ['B','E']],
    ['FE', ['F','E']],
    ['CE', ['C','E']],
    ['CD', ['C','D']],
    ['ED', ['E','D']],
    #forces at the support joints
    ['Ay', ['A', 'Fy']],
    ['Dx', ['D', 'Fx']],
    ['Dy', ['D', 'Fy']]
  ],true
)
v = truss.solve
print "Truss = ", truss.joint_to_s, "\n"
#print "Solution Vector = \n", truss.solution_to_s, "\n"
#print "Test = " , truss.test_solution, "\n"
print "Forces\n", truss.forces_to_s, "\n"
