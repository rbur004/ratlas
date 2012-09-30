require "../truss"
def bd(k)
  #dummy call to keep this the same as the simultaneous.rb example
  k
end

print "Tensions are positive, Compression negative\n"
print "0 point is floating end, Forces have opposite sign to members.\n"
print "Forces at supports can be added as Ay (floating end), Dx, Dy (fixed end). Set mount_forces to true\n\n"


print "Truss from page 14,15 The analysis of Engineering Structures. Pippard & Barker Third edition\n"
print "Should give results in table on page 16\n"
print "3D test of truss code using method of tension coefficients.\n\n"
truss = Truss.new(
[['A', [bd(-20.0), bd(8.66), bd(-5.0)], [bd(0.0), bd(0.0), bd(0.0)]],
 ['B', [bd(-10.0), bd(6.93), bd(-4.0)], [bd(0.0), bd(0.0), bd(0.0)]],
 ['C', [bd(0.0),   bd(5.2),  bd(-3.0)], [bd(0.0), bd(6.0), bd(0.0)]],
 ['D', [bd(-20.0), bd(0.0),  bd(0.0)], [bd(0.0), bd(0.0), bd(0.0)]],
 ['E', [bd(-10.0), bd(0.0),  bd(0.0)], [bd(0.0), bd(0.0), bd(0.0)]],
 ['F', [bd(0.0),   bd(0.0),  bd(0.0)], [bd(0.0), bd(10.0), bd(0.0)]],
 ['G', [bd(-20.0), bd(8.66), bd(5.0)], [bd(0.0), bd(0.0), bd(0.0)]],
 ['H', [bd(-10.0), bd(6.93), bd(4.0)], [bd(0.0), bd(0.0), bd(0.0)]],
 ['J',  [bd(0.0),  bd(5.2),  bd(3.0)], [bd(0.0), bd(3.0), bd(0.0)]]
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
 ['HG', ['H', 'G']]
             ],
             false
            )
    
#print "Matrix = \n", truss.matrix_to_s, "\n"
v = truss.solve
print "Force Vector = \n", truss.joint_to_s, "\n"
#print "Solution Vector = \n", truss.solution_to_s, "\n"
print "Forces\n", truss.forces_to_s, "\n"
#print "Test = " , truss.test_solution, "\n"

