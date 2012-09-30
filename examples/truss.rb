require 'ratlas'
include RAtlas
include Math

class Lapack
  def pp
    s = "[\n"
    self.each_row do |r,i|
      s += "  ["
      r.each { |c| s += "#{sprintf("%2.2f", c)}, "}
      s[-1] = "]"
      s += "\n"
    end
    s += "]"
  end
end

class Blas
  def x
    self[0]
  end
  def y
    self[1]
  end
  def z
    self[2]
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
  
  def to_bigdecimal
    self.collect! { |v|  BigDecimal(v.to_s) }
  end
end
  
class Joint 
  attr_accessor :name, :coord, :load
  def initialize(name, joint, load)
    @name = name
    @coord = DoubleLapack[*joint]
    @load = DoubleLapack[*load] 
  end
  
  def distance(joint)
      if joint.class == Force
        joint.coord
      else
        joint.coord - self.coord
      end
  end
    
  def to_s
    #print "name *** = #{@name}\n"
    axis = ["x","y","z"]
    s = "#{@name}"
    (0...@coord.nrows).each do |i|
      s += " #{axis[i]} = #{@coord[i]}"
    end
    s += " Load Force"
    (0...@load.nrows).each do |i|
      s += " #{axis[i]} = #{@load[i]}"
    end
    return s
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
  
end

class Truss
  attr_reader :joints,  :members, :dimensions, :solution, :force_vector
  
  def initialize( joints, members, mount_forces = true )
    @joints  = {}
    @members = []
    @dimensions = joints[0][1].length
    @mount_forces = mount_forces #these forces at the mount points need calculating?
    
    #puts members.pp
    
    joints.each do |j| 
      @joints[ j[0] ] = Joint.new(*j)
      #print @joints[ j[0] ], "\n"
    end
    
    #print "     "
    members.each do |m|
      #print "    #{m[0]}, #{m[1][0]}, #{m[1][1]}"
      if(@dimension == 2)
        case m[1][1]
        when 'Fx' : @members << Member.new( m[0], @joints[ m[1][0] ], Force.new(m[0],[1.0, 0.0], nil) )
        when 'Fy' : @members << Member.new( m[0], @joints[ m[1][0] ], Force.new(m[0],[0.0, 1.0], nil) )
        else @members << Member.new( m[0], @joints[ m[1][0] ], @joints[ m[1][1] ] )
        end
      else
        case m[1][1]
        when 'Fx' : @members << Member.new( m[0], @joints[ m[1][0] ], Force.new(m[0],[1.0, 0.0, 0.0], nil) )
        when 'Fy' : @members << Member.new( m[0], @joints[ m[1][0] ], Force.new(m[0],[0.0, 1.0, 0.0], nil) )
        when 'Fz' : @members << Member.new( m[0], @joints[ m[1][0] ], Force.new(m[0],[0.0, 0.0, 1.0], nil) )
        else @members << Member.new( m[0], @joints[ m[1][0] ], @joints[ m[1][1] ] )
        end
      end
    end
    #print "\n"
    
    create_matrix
    
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
        @matrix =  DoubleLapack[*(m[0...i] + m[i+1..-1])]
        @force_vector = DoubleLapack[*(v[0...i] + v[i+1..-1])]
        break if @matrix.determinant != 0
      end
    else
      @matrix =  DoubleLapack[ *m]
      @force_vector = DoubleLapack[ *v]
    end
    #print @matrix
    #print "determinant is 0\n" if @matrix.determinant == 0
    #print "#{@matrix.row_size} x #{@matrix.column_size}\n"
  end
  
  
  def solve
    #Use rlapack routines (i.e. the ruby clpack interfaces to the f77 lapack routines provided by atlas)
    p = IntegerBlas.new(@force_vector.nrows)
    #puts @matrix.pp
    #puts @force_vector.pp
    @solution = @matrix.xgesv(@force_vector, p)
=begin
    #Basic determinant method
    print @matrix.determinant, "\n", @matrix.pp, "\n"
    @solution = @matrix.inverse * @force_vector
=
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
    @solution.pp
  end
  
  def matrix_to_s
    @matrix.pp
  end

  def force_vector_to_s
    @force_vector.pp
  end
  
  def test_solution(sol = @solution)
    s = @matrix.determinant.to_s
    result = @matrix * sol
    s += result.to_s
  end
  
  def forces_to_s
    s = ""
    @members.each_with_index do|m,i|
      s += "Member #{m.name}, Tension #{@solution[i]}, length #{m.length.length}, force = #{m.length.length*@solution[i]} \n"
    end
    return s
  end
  
  def joint_to_s    
    s = ""
    @joints.each { |j,k| s += k.to_s + "\n"}
    return s
  end
end

