require "../truss"

print "Tensions are positive, Compression negative\n"
print "0 point is floating end, Forces have opposite sign to members.\n"
print "Forces at supports can be added as Ay (floating end), Dx, Dy (fixed end). Set mount_forces to true\n\n"

print "Eng Mech Structures 4th Ed. Meriam & Kraige. pg 185 Prob 4/9\n\n"
height = Math::sqrt(8*8-4*4)
eweight = (400*9.81/2*3) #and D weight
bweight = (400*9.81/2*4) 
aweight = (400/2*9.81*2) #and C Weight
print "total at A = #{7*400/2*9.81}\n"
truss = Truss.new(
  [
      ['A', [0,0], [0,aweight]],
      ['B', [8,0], [0,bweight]],
      ['C', [16,0], [0,aweight]],
      ['D', [12, height], [0,eweight]],
      ['E', [4, height], [0,eweight]],
  ],
  [
    ['AB', ['A', 'B']],
    ['AE', ['A', 'E']],
    ['BC', ['B','C']],
    ['BD', ['B','D']],
    ['BE', ['B','E']],
    ['CD', ['C','D']],
    ['ED', ['E','D']],
    #supports
    ['Cy', ['C', 'Fy']],
    ['Cx', ['C', 'Fx']],
    ['Ay', ['A', 'Fy']]
  ],true
)
v = truss.solve
print "Truss = ", truss.joint_to_s, "\n"
#print "Solution Vector = \n", truss.solution_to_s, "\n"
#print "Test = " , truss.test_solution, "\n"
print "Forces\n", truss.forces_to_s, "\n"

