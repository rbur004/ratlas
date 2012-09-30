require "../truss"

print "Tensions are positive, Compression negative\n"
print "0 point is floating end, Forces have opposite sign to members.\n"
print "Forces at supports can be added as Ay (floating end), Dx, Dy (fixed end). Set mount_forces to true\n\n"


print "Rachel's truss 2\n\n"
joints = []
members = []
 (0..12).each do |i|
   joints[i] = [sprintf("L%02d", i+1), [i*2.4, 8.6-i*0.6], [0, 0.175]]
   joints[i+13] = [sprintf("U%02d", i+1), [i*2.4+0.2, 9.5-i*0.6], [0,-0.6]]
 end
 joints[5+13][2][1] = 0.825
 joints[6+13][2][1] = 1.825
 joints[7+13][2][1] = 0.825
 
 
#joints.each {|i| print "#{i[0]} #{i[1][0]}, #{i[1][1]} Force [#{i[2][0]}, #{i[2][1]}]\n"}
 
 m = -1
 (0..11).each do |i|
   j = i + 13
   members[m+=1] = [ "#{joints[i][0]}-#{joints[i+1][0]}", ["#{joints[i][0]}", "#{joints[i+1][0]}"] ]
   members[m+=1] = [ "#{joints[i][0]}-#{joints[j][0]}", ["#{joints[i][0]}", "#{joints[j][0]}"] ]
   members[m+=1] = [ "#{joints[j][0]}-#{joints[j+1][0]}", ["#{joints[j][0]}", "#{joints[j+1][0]}"] ]
   if(i%2 == 1)
     members[m+=1] = [ "#{joints[j][0]}-#{joints[i+1][0]}", ["#{joints[j][0]}", "#{joints[i+1][0]}"] ]
   else
     members[m+=1] = [ "#{joints[i][0]}-#{joints[j+1][0]}", ["#{joints[i][0]}", "#{joints[j+1][0]}"] ]
   end
 end
 members[m+=1] = [ "#{joints[12][0]}-#{joints[12+13][0]}", ["#{joints[12][0]}", "#{joints[12+13][0]}"]]
 
 members[m+=1] = ['L01y', ['L01', 'Fy']] #fixed at the top
 members[m+=1] = ['L13x', ['L13', 'Fx']] #roller at the bottom
 members[m+=1] = ['L13y', ['L13', 'Fy']]

#members.each { |m| print "#{m[0]} [ #{m[1][0]}, #{m[1][1]}]\n"}
 
truss = Truss.new( joints, members, true)

v = truss.solve
print "Truss = \n", truss.joint_to_s, "\n"
#print "Solution Vector = \n", truss.solution_to_s, "\n"
#print "Test = " , truss.test_solution, "\n"
print "Forces\n", truss.forces_to_s, "\n"


