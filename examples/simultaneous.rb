require 'mathn'
require 'bigdecimal/util'
require 'bigdecimal/math'
include BigMath


class Matrix
  def pp
    s = "[\n"
    self.row_vectors.each do |r|
      s += "  ["
      r.collect { |c| s += "#{sprintf("%2.2f", c)}, "}
      s += "]\n"
    end
    s += "]"
  end
end

class Array #pretty print add on for arrays
  def pp(format = nil)
    s =  "[ "
    if format
      s += sprintf(format, self[0]) if self[0]
    else
      s += "#{self[0]}" if self[0]
    end
    self[1..-1].each do |x| 
      if format
        s += ( ", " + sprintf(format, x) )
      else
        s += ", #{x}" 
      end
    end
    s += " ]"
  end
end

class Vector #add on retrieve by x,y,z
	def x
	  self[0]
	end
	def y
	  self[1]
  end
  def z
    self[2]
  end
  def length
    Math::sqrt(self.to_a.inject(0) { |i,j| i += j*j})
  end
end

    
class Joint 
  attr_accessor :name, :coord, :load
  def initialize(name, joint, load)    
    @name = name
    @coord = Vector[*joint]
    @load = Vector[*load] 
  end
  
  def distance(joint)
      if joint.class == Force
        joint.coord
      else
        joint.coord - self.coord
      end
  end
  
  def to_s
    "#{@name} => #{@coord} #{@load}"
  end
end

class Force < Joint
#dummy class to make unknown forces at a joint have the same structure as a joint definition
end

class Member
  attr_accessor :name, :length, :start_coord, :end_coord
  
  def initialize(name, start_coord, end_coord)
    @name = name
    @start_coord = start_coord
    @end_coord = end_coord
    @length = start_coord.distance(end_coord)
  end
  
  def linear_length
    @length.length
  end
  
end

class Truss
  attr_reader :joints,  :members, :dimensions, :solution, :force_vector
  
  def initialize( joints, members, mount_forces = true )
    @joints  = {}
    @members = []
    @dimensions = joints[0][1].length
    @mount_forces = mount_forces #these forces at the mount points need calculating?
    
    joints.each do |j| 
      @joints[ j[0] ] = Joint.new(*j)
      #print @joints[ j[0] ], "\n"
    end
    
    #print "     "
    members.each do |m|
      #print "    #{m[0]}, "
      case m[1][1]
      when 'Fx' : @members << Member.new( m[0], @joints[ m[1][0] ], Force.new(m[0],[bd(1.0), bd(0.0), bd(0.0)],nil) )
      when 'Fy' : @members << Member.new( m[0], @joints[ m[1][0] ], Force.new(m[0],[bd(0.0), bd(1.0), bd(0.0)],nil) )
      when 'Fz' : @members << Member.new( m[0], @joints[ m[1][0] ], Force.new(m[0],[bd(0.0), bd(0.0), bd(1.0)],nil) )
      else @members << Member.new( m[0], @joints[ m[1][0] ], @joints[ m[1][1] ] )
      end
    end
    #print "\n"
    
    create_matrix
    
  end
  
  def bd(f)
    #BigDecimal(f.to_s)
    f
  end
  
  def create_matrix
    pre_matrix = {} #we are creating a hash of arrays, each being an equation to solve
    pre_force_vector = {}
    
    if @mount_forces
      matrix_size = @joints.length * @dimensions 
      @joints.each do |j,k|
        pre_matrix[k.name + 'x'] = Array.new(matrix_size,0.0) 
        pre_matrix[k.name + 'y'] = Array.new(matrix_size,0.0)

        pre_force_vector[k.name + 'x'] = k.load.x
        pre_force_vector[k.name + 'y'] = k.load.y

        if @dimensions == 3
           pre_matrix[k.name + 'z'] = Array.new(matrix_size,0.0)
           pre_force_vector[k.name + 'z'] = k.load.z
        end
      end
    else
      matrix_size = @members.length 
      @members.each do |m|
        pre_matrix[m.start_coord.name + 'x'] = Array.new(matrix_size,0.0) 
        pre_matrix[m.start_coord.name + 'y'] = Array.new(matrix_size,0.0)

        pre_force_vector[m.start_coord.name + 'x'] = m.start_coord.load.x
        pre_force_vector[m.start_coord.name + 'y'] = m.start_coord.load.y

        if @dimensions == 3
           pre_matrix[m.start_coord.name + 'z'] = Array.new(matrix_size,0.0)
           pre_force_vector[m.start_coord.name + 'z'] = m.start_coord.load.z
        end
      end
    end
    
    @members.each_with_index do |m,j| 
        pre_matrix[m.start_coord.name + 'x'][j] = m.length.x
        pre_matrix[m.start_coord.name + 'y'][j] = m.length.y
        pre_matrix[m.start_coord.name + 'z'][j] = m.length.z if @dimensions == 3
        
         
        if pre_matrix[m.end_coord.name + 'x']  #ignore the redundand possibilities
          pre_matrix[m.end_coord.name + 'x'][j] = -m.length.x
          pre_matrix[m.end_coord.name + 'y'][j] = -m.length.y
          pre_matrix[m.end_coord.name + 'z'][j] = -m.length.z if @dimensions == 3
        end      
    end
    
    m = []
    v = []
    @row_names = []
    pre_matrix.sort.each do |k,row|
      #print "#{k}, #{row.pp("%6.3f")} = #{pre_force_vector[k]}\n"
      m << row
      v << pre_force_vector[k]
      @row_names << k
    end

    if @dimensions == 2 && @mount_forces == false
      (0...m.length - 1).each do |i|
        @matrix =  Matrix.rows(m[0...i] + m[i+1..-1])
        @force_vector = Vector.elements(v[0...i] + v[i+1..-1])
        break if @matrix.determinant != 0
      end
    else
      @matrix =  Matrix.rows(m)
      @force_vector = Vector.elements(v)
    end
    print "determinant is 0\n" if @matrix.determinant == 0
    #print "#{@matrix.row_size} x #{@matrix.column_size}\n"
  end
  
  
  def solve
    print @matrix.determinant, "\n", @matrix.pp, "\n"
    @solution = @matrix.inverse * @force_vector
=begin
    #alternate method of solving simultaneous equations
    #attemping to find why we are getting errors in matrix inversion.
    #print @matrix.pp, "\n"
    d = @matrix.determinant
    l = @force_vector.size
    
    m_array = []
    @matrix.column_vectors.each do |r|
      col = []
      r.collect { |c| col << c }
      m_array << col
    end
    
    f_array = []
    @force_vector.collect { |c| f_array << c }
      
    solution = []
    (0...l).each do |i|
      tmp_matrix = Matrix.columns(m_array[0...i] + [f_array] + m_array[i+1..-1] )
      #print tmp_matrix.pp, "\n"
      d2 = tmp_matrix.determinant
      print d2, "/", d, "= #{d2/d}\n"
      solution << d2/d
    end
    @solution = Vector.elements(solution)
=end
  end
  
  def solution_to_s
    print "["
    result.collect { |c| print "#{sprintf("%6.4f", c)}, "}
    print "]\n"
  end
  
  def matrix_to_s
    s = "[\n"
    @matrix.row_vectors.each do |r|
      s += "  ["
      r.collect { |c| s += "#{sprintf("%6.4f", c)}, "}
      s += "]\n"
    end
    s += "]"
  end

  def force_vector_to_s
    s = "\n "
    @row_names.each do |n|
      s += "    #{n}, "
    end
    s += "\n["
    @force_vector.collect { |c| s += "#{sprintf("%6.4f", c)}, "}
    s += "]"
  end

  def solution_to_s
    s = "\n["
    result.collect { |c| s += "#{sprintf("%6.4f", c)}, "}
    s += "]"
  end
  
  def test_solution(sol = @solution)
    s = @matrix.determinant.to_s
    result = @matrix * sol
    s += "\n["
    result.collect { |c| s += "#{sprintf("%6.4f", c)}, "}
    s += "]\n"
  end
  
  def forces_to_s
    s = ""
    @members.each_with_index do|m,i|
      s += "Member #{m.name}, length #{m.length.length}, force = #{m.length.length*@solution[i]} \n"
    end
    return s
  end
  
end

def bd(f)
  #BigDecimal(f.to_s)
  f
end

print "Statics Book problem 4/29\n"
truss = Truss.new(
  [
      ['A', [-12,0], [0,0]],
      ['B', [-9,4], [0,0]],
      ['C', [-3,4], [0,0]],
      ['D', [0, 0], [0, 0]],
      ['E', [-3, 0], [0,-6]],
      ['F', [-9, 0], [0,-4]],
  ],
  [
    ['DC', ['D','C']],
    ['DE', ['D','E']],
    ['CE', ['C','E']],
    ['CB', ['C','B']],
    ['EB', ['E','B']],
    ['EF', ['E','F']],
    ['BF', ['B','F']],
    ['FA', ['F', 'A']],
    ['BA', ['B', 'A']]
  ],
  false
)
v = truss.solve
print "Solution\n", truss.forces_to_s

#print "Force Vector = \n", truss.force_vector_to_s, "\n"
#print "Solution Vector = \n", truss.solution, "\n"
#print "Forces\n", truss.print_forces, "\n"
#print "Test = " , truss.test_solution, "\n"

=begin
print "3D test of truss code. Value from book example 1.\n"
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
    
v = truss.solve
#print "Matrix = \n", truss.matrix_to_s, "\n"
print "Force Vector = ", truss.force_vector_to_s, "\n"
print "Test = " , truss.test_solution, "\n"

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

print "Truss from text, with Bx, By and Cy already given, but also asked for\n"
print "Hence the Bx, By, and Cy value will be 0 (or close to 0)\n"
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
print "Force Vector = ", truss2.force_vector_to_s, "\n"
print "Test = " , truss2.test_solution, "\n"
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
v = truss2a.solve
print "Force Vector = ", truss2a.force_vector_to_s, "\n"
print "Test = " , truss2a.test_solution, "\n"
exit
#=begin

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


