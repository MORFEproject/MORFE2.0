"""
Overview: collections of shape functions and their derivatives
"""


function select_quadrature_rule(etype::Symbol)
  #
  if (etype==:TE4)
    return get_quadrature_points(Val{:GLTET15})
  elseif (etype==:T10)
    return get_quadrature_points(Val{:GLTET15})
  elseif (etype==:HE8)
    return get_quadrature_points(Val{:GLHEX8})
  elseif (etype==:H20)
    return get_quadrature_points(Val{:GLHEX27})
  elseif (etype==:H27)
    return get_quadrature_points(Val{:GLHEX27})
  elseif (etype==:PE6)
    return get_quadrature_points(Val{:GLWED21})
  elseif (etype==:P15)
    return get_quadrature_points(Val{:GLWED21})
  elseif (etype==:P18)
    return get_quadrature_points(Val{:GLWED21})
  end
  #
  return nothing
end



"""
> N!(N::Array{Float64},gp::Tuple,::Type{Val{:T10}})
It fills N with the shape functions of a T10 element evaluated at gauss point gp.
- N : shape functions Array
- gp : gauss point coordinates
- :Type{Val{:T10}} : element type specifier
"""
function N!(N::Array{Float64},gp::Tuple,::Type{Val{:T10}})
  #
  gp4 = 1.0-gp[1]-gp[2]-gp[3]
  #

  N[1] = gp[1]*(2.0*gp[1]-1)
  N[2] = gp[2]*(2.0*gp[2]-1)
  N[3] = gp[3]*(2.0*gp[3]-1)
  N[4] = gp4*(2.0*gp4-1)
  #
  N[5]  =  4.0*gp[1]*gp[2]
  N[6]  =  4.0*gp[2]*gp[3]
  N[7]  =  4.0*gp[1]*gp[3]
  N[8]  =  4.0*gp[1]*gp4
  N[9]  =  4.0*gp[2]*gp4
  N[10] =  4.0*gp[3]*gp4
  #
  return nothing
end


"""
> ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:T10}})
It fills ∂N∂a with the derivatives of the shape functions of a T4 element evaluated at gauss point gp.
- ∂N∂a : shape functions derivatives Array
- gp : gauss point coordinates
- :Type{Val{:T10}} : element type specifier
"""
function ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:T10}})
  #
  qa1 = 4.0*gp[1]
  qa2 = 4.0*gp[2]
  qa3 = 4.0*gp[3]
  qa4 = 4.0*(1-gp[1]-gp[2]-gp[3])
  #
  ∂N∂a[1,1] = -1.0+qa1
  ∂N∂a[1,2] = 0.0
  ∂N∂a[1,3] = 0.0
  #
  ∂N∂a[2,1] = 0.0
  ∂N∂a[2,2] = -1.0+qa2
  ∂N∂a[2,3] = 0.0
  #
  ∂N∂a[3,1] = 0.0
  ∂N∂a[3,2] = 0.0
  ∂N∂a[3,3] = -1.0+qa3
  #
  ∂N∂a[4,1] = 1.0-qa4
  ∂N∂a[4,2] = 1.0-qa4
  ∂N∂a[4,3] = 1.0-qa4
  #
  ∂N∂a[5,1] = qa2
  ∂N∂a[5,2] = qa1
  ∂N∂a[5,3] = 0.0
  #
  ∂N∂a[6,1] = 0.0
  ∂N∂a[6,2] = qa3
  ∂N∂a[6,3] = qa2
  #
  ∂N∂a[7,1] = qa3
  ∂N∂a[7,2] = 0.0
  ∂N∂a[7,3] = qa1
  #
  ∂N∂a[8,1] = qa4-qa1
  ∂N∂a[8,2] = -qa1
  ∂N∂a[8,3] = -qa1
  #
  ∂N∂a[9,1] = -qa2
  ∂N∂a[9,2] = qa4-qa2
  ∂N∂a[9,3] = -qa2
  #
  ∂N∂a[10,1] = -qa3
  ∂N∂a[10,2] = -qa3
  ∂N∂a[10,3] = qa4-qa3
  #
  return nothing
  #
end




"""
> N!(N::Array{Float64},gp::Tuple,::Type{Val{:HE8}})
It fills N with the shape functions of a H8 element evaluated at gauss point gp.
- N : shape functions Array
- gp : gauss point coordinates
- :Type{Val{:HE8}} : element type specifier
"""
function N!(N::Array{Float64},gp::Tuple,::Type{Val{:HE8}})
  #
  x,y,z = gp
  # -1, +1, -1
  N[1] = 0.125*(1 - x)*(1 + y)*(1 - z)
  # +1, +1, -1
  N[2] = 0.125*(1 + x)*(1 + y)*(1 - z)
  # +1, -1, -1
  N[3] = 0.125*(1 + x)*(1 - y)*(1 - z)
  # -1, -1, -1
  N[4] = 0.125*(1 - x)*(1 - y)*(1 - z)
  # -1, +1, +1
  N[5] = 0.125*(1 - x)*(1 + y)*(1 + z)
  # +1, +1, +1
  N[6] = 0.125*(1 + x)*(1 + y)*(1 + z)
  # +1, -1, +1
  N[7] = 0.125*(1 + x)*(1 - y)*(1 + z)
  # -1, -1, +1
  N[8] = 0.125*(1 - x)*(1 - y)*(1 + z)
  #
  return nothing
end



"""
> ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:HE8}})
It fills ∂N∂a with the derivatives of the shape functions of a H8 element evaluated at gauss point gp.
- ∂N∂a : shape functions derivatives Array
- gp : gauss point coordinates
- :Type{Val{:HE8}} : element type specifier
"""
function ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:HE8}})
  #
  x,y,z = gp
  #
  ∂N∂a[1,1] = -0.125*(1 + y)*(1 - z)
  ∂N∂a[1,2] = +0.125*(1 - x)*(1 - z)
  ∂N∂a[1,3] = -0.125*(1 - x)*(1 + y)
  #
  ∂N∂a[2,1] = +0.125*(1 + y)*(1 - z)
  ∂N∂a[2,2] = +0.125*(1 + x)*(1 - z)
  ∂N∂a[2,3] = -0.125*(1 + x)*(1 + y)
  #
  ∂N∂a[3,1] = +0.125*(1 - y)*(1 - z)
  ∂N∂a[3,2] = -0.125*(1 + x)*(1 - z)
  ∂N∂a[3,3] = -0.125*(1 + x)*(1 - y)
  #
  ∂N∂a[4,1] = -0.125*(1 - y)*(1 - z)
  ∂N∂a[4,2] = -0.125*(1 - x)*(1 - z)
  ∂N∂a[4,3] = -0.125*(1 - x)*(1 - y)
  #
  ∂N∂a[5,1] = -0.125*(1 + y)*(1 + z)
  ∂N∂a[5,2] = +0.125*(1 - x)*(1 + z)
  ∂N∂a[5,3] = +0.125*(1 - x)*(1 + y)
  #
  ∂N∂a[6,1] = +0.125*(1 + y)*(1 + z)
  ∂N∂a[6,2] = +0.125*(1 + x)*(1 + z)
  ∂N∂a[6,3] = +0.125*(1 + x)*(1 + y)
  #
  ∂N∂a[7,1] = +0.125*(1 - y)*(1 + z)
  ∂N∂a[7,2] = -0.125*(1 + x)*(1 + z)
  ∂N∂a[7,3] = +0.125*(1 + x)*(1 - y)
  #
  ∂N∂a[8,1] = -0.125*(1 - y)*(1 + z)
  ∂N∂a[8,2] = -0.125*(1 - x)*(1 + z)
  ∂N∂a[8,3] = +0.125*(1 - x)*(1 - y)
  #
  return nothing
  #
end





"""
> N!(N::Array{Float64},gp::Tuple,::Type{Val{:TE4}})
It fills N with the shape functions of a T4 element evaluated at gauss point gp.
- N : shape functions Array
- gp : gauss point coordinates
- :Type{Val{:TE4}} : element type specifier
"""
function N!(N::Array{Float64},gp::Tuple,::Type{Val{:TE4}})
  #
  N[1] = gp[1]
  N[2] = gp[2]
  N[3] = gp[3]
  N[4] = 1-gp[1]-gp[2]-gp[3]
  #
  return nothing
end





"""
> ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:TE4}})
It fills ∂N∂a with the derivatives of the shape functions of a T4 element evaluated at gauss point gp.
- ∂N∂a : shape functions derivatives Array
- gp : gauss point coordinates
- :Type{Val{:TE4}} : element type specifier
"""
function ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:TE4}})
  #
  ∂N∂a[1,1] = 1.0
  ∂N∂a[1,2] = 0.0
  ∂N∂a[1,3] = 0.0
  #
  ∂N∂a[2,1] = 0.0
  ∂N∂a[2,2] = 1.0
  ∂N∂a[2,3] = 0.0
  #
  ∂N∂a[3,1] = 0.0
  ∂N∂a[3,2] = 0.0
  ∂N∂a[3,3] = 1.0
  #
  ∂N∂a[4,1] = -1.0
  ∂N∂a[4,2] = -1.0
  ∂N∂a[4,3] = -1.0
  #
  return nothing
  #
end



"""
> N!(N::Array{Float64},gp::Tuple,::Type{Val{:PE6}})
It fills N with the shape functions of a PE6 element evaluated at gauss point gp.
- N : shape functions Array
- gp : gauss point coordinates
- :Type{Val{:PE6}} : element type specifier
"""
function N!(N::Array{Float64},gp::Tuple,::Type{Val{:PE6}})
  #
  x,y,z = gp
  #
  N[1] = x * 0.5 * (1.0 - z)
  N[2] = y * 0.5 * (1.0 - z)
  N[3] = (1.0-x-y) * 0.5 * (1.0 - z)
  N[4] = x * 0.5 * (1.0 + z)
  N[5] = y * 0.5 * (1.0 + z)
  N[6] = (1.0-x-y) * 0.5 * (1.0 + z)
  #
  return nothing
end



"""
> ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:PE6}})
It fills ∂N∂a with the derivatives of the shape functions of a T4 element evaluated at gauss point gp.
- ∂N∂a : shape functions derivatives Array
- gp : gauss point coordinates
- :Type{Val{:PE6}} : element type specifier
"""
function ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:PE6}})
  #
  x,y,z = gp
  #
  ∂N∂a[1,1] = 0.5*(1 - z)
  ∂N∂a[1,2] = 0.0
  ∂N∂a[1,3] = -0.5*x
  #
  ∂N∂a[2,1] = 0.0
  ∂N∂a[2,2] = 0.5*(1 - z)
  ∂N∂a[2,3] = -0.5*y
  #
  ∂N∂a[3,1] = -0.5*(1 - z)
  ∂N∂a[3,2] = -0.5*(1 - z)
  ∂N∂a[3,3] = -0.5*(1 - x - y)
  # 
  ∂N∂a[4,1] = +0.5*(1 + z)
  ∂N∂a[4,2] = +0.0
  ∂N∂a[4,3] = +0.5*x
  #
  ∂N∂a[5,1] = +0.0
  ∂N∂a[5,2] = +0.5*(1 + z)
  ∂N∂a[5,3] = +0.5*y
  #
  ∂N∂a[6,1] = -0.5*(1 + z)
  ∂N∂a[6,2] = -0.5*(1 + z)
  ∂N∂a[6,3] = +0.5*(1 - x - y)
  #
  return nothing
  #
end




"""
> N!(N::Array{Float64},gp::Tuple,::Type{Val{:P15}})
It fills N with the shape functions of a PE6 element evaluated at gauss point gp.
- N : shape functions Array
- gp : gauss point coordinates
- :Type{Val{:P15}} : element type specifier
"""
function N!(N::Array{Float64},gp::Tuple,::Type{Val{:P15}})
    #
    @inbounds c1 = 1.0 - gp[1] - gp[2]
    @inbounds c2 = 1.0 - gp[3]
    @inbounds c3 = 1.0 + gp[3]
    #
    @inbounds da1 = 2.0*gp[1]
    @inbounds da2 = 2.0*gp[2]
    #
    @inbounds N[1] = -c1*c2*(da1+da2+gp[3])/2.0
    @inbounds N[2] = -gp[1]*c2*(2.0-da1+gp[3])/2.0
    @inbounds N[3] = -gp[2]*c2*(2.0-da2+gp[3])/2.0
    @inbounds N[4] = -c1*c3*(da1+da2-gp[3])/2.0
    @inbounds N[5] = -gp[1]*c3*(2.0-da1-gp[3])/2.0
    @inbounds N[6] = -gp[2]*c3*(2.0-da2-gp[3])/2.0
    @inbounds N[7] = da1*c1*c2
    @inbounds N[8] = da1*gp[2]*c2
    @inbounds N[9] = da2*c1*c2
    @inbounds N[10] = da1*c1*c3
    @inbounds N[11] = da1*gp[2]*c3
    @inbounds N[12] = da2*c1*c3
    @inbounds N[13] = c1*c2*c3
    @inbounds N[14] = gp[1]*c2*c3
    @inbounds N[15] = gp[2]*c2*c3
    return nothing
end



"""
> ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:P15}})
It fills ∂N∂a with the derivatives of the shape functions of a T4 element evaluated at gauss point gp.
- ∂N∂a : shape functions derivatives Array
- gp : gauss point coordinates
- :Type{Val{:P15}} : element type specifier
"""
function ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:P15}})
  @inbounds  c1 = 1.0 - gp[1] - gp[2]
  @inbounds  c2 = 1.0 - gp[3]
  @inbounds  c3 = 1.0 + gp[3]
  #
  @inbounds  da1 = 2.0*gp[1]
  @inbounds  da2 = 2.0*gp[2]
  @inbounds  da3 = 2.0*gp[3]
  #
  @inbounds  qa1 = 4.0*gp[1]
  @inbounds  qa2 = 4.0*gp[2]
  # 
  @inbounds  ∂N∂a[1,1] = c2*(qa1+qa2+gp[3]-2.0)/2.0
  @inbounds  ∂N∂a[1,2] = c2*(qa1+qa2+gp[3]-2.0)/2.0
  @inbounds  ∂N∂a[1,3] = c1*(da1+da2+da3-1.0)/2.0
  # 
  @inbounds ∂N∂a[2,1] = c2*(qa1-gp[3]-2.0)/2.0
  @inbounds ∂N∂a[2,2] = 0.0
  @inbounds ∂N∂a[2,3] = gp[1]*(-da1+da3+1.0)/2.0
  # 
  @inbounds ∂N∂a[3,1] = 0.0
  @inbounds ∂N∂a[3,2] = c2*(qa2-gp[3]-2.0)/2.0
  @inbounds ∂N∂a[3,3] = gp[2]*(-da2+da3+1.0)/2.0
  # 
  @inbounds ∂N∂a[4,1] = c3*(qa1+qa2-gp[3]-2.0)/2.0
  @inbounds ∂N∂a[4,2] = c3*(qa1+qa2-gp[3]-2.0)/2.0
  @inbounds ∂N∂a[4,3] = -c1*(da1+da2-da3-1.0)/2.0
  # 
  @inbounds ∂N∂a[5,1] = c3*(qa1+gp[3]-2.0)/2.0
  @inbounds ∂N∂a[5,2] = 0.0
  @inbounds ∂N∂a[5,3] = gp[1]*(da1+da3-1.0)/2.0
  # 
  @inbounds ∂N∂a[6,1] = 0.0
  @inbounds ∂N∂a[6,2] = c3*(qa2+gp[3]-2.0)/2.0
  @inbounds ∂N∂a[6,3] = gp[2]*(da2+da3-1.0)/2.0
  # 
  @inbounds ∂N∂a[7,1] = -c2*(da1+gp[2]-1.0)*2.0
  @inbounds ∂N∂a[7,2] = -gp[1]*c2*2.0
  @inbounds ∂N∂a[7,3] = -gp[1]*c1*2.0
  # 
  @inbounds ∂N∂a[8,1] = gp[2]*c2*2.0
  @inbounds ∂N∂a[8,2] = gp[1]*c2*2.0
  @inbounds ∂N∂a[8,3] = -gp[1]*gp[2]*2.0
  # 
  @inbounds ∂N∂a[9,1] = -gp[2]*c2*2.0
  @inbounds ∂N∂a[9,2] = -c2*(gp[1]+da2-1.0)*2.0
  @inbounds ∂N∂a[9,3] = -gp[2]*c1*2.0
  # 
  @inbounds ∂N∂a[10,1] = -c3*(da1+gp[2]-1.0)*2.0
  @inbounds ∂N∂a[10,2] = -gp[1]*c3*2.0
  @inbounds ∂N∂a[10,3] = gp[1]*c1*2.0
  # 
  @inbounds ∂N∂a[11,1] = gp[2]*c3*2.0
  @inbounds ∂N∂a[11,2] = gp[1]*c3*2.0
  @inbounds ∂N∂a[11,3] = gp[1]*gp[2]*2.0
  # 
  @inbounds ∂N∂a[12,1] = -gp[2]*c3*2.0
  @inbounds ∂N∂a[12,2] = -c3*(gp[1]+da2-1.0)*2.0
  @inbounds ∂N∂a[12,3] = gp[2]*c1*2.0
  # 
  @inbounds ∂N∂a[13,1] = -c2*c3
  @inbounds ∂N∂a[13,2] = -c2*c3
  @inbounds ∂N∂a[13,3] = -da3*c1
  # 
  @inbounds ∂N∂a[14,1] = c2*c3
  @inbounds ∂N∂a[14,2] = 0.0
  @inbounds ∂N∂a[14,3] = -da1*gp[3]
  # 
  @inbounds ∂N∂a[15,1] = 0.0
  @inbounds ∂N∂a[15,2] = c2*c3
  @inbounds ∂N∂a[15,3] = -da2*gp[3]
  #
  return nothing
  #
end


"""
> N!(N::Array{Float64},gp::Tuple,::Type{Val{:H27}})
It fills N with the shape functions of a PE6 element evaluated at gauss point gp.
- N : shape functions Array
- gp : gauss point coordinates
- :Type{Val{:H27}} : element type specifier
"""
function N!(N::Array{Float64},gp::Tuple,::Type{Val{:H27}})
  #
  x1 = gp[1]*(gp[1]-1.0)/2.0
  x2 = (1-gp[1])*(1+gp[1])
  x3 = gp[1]*(gp[1]+1.0)/2.0
  #
  y1 = gp[2]*(gp[2]-1.0)/2.0
  y2 = (1-gp[2])*(1+gp[2])
  y3 = gp[2]*(gp[2]+1.0)/2.0
  #
  z1 = gp[3]*(gp[3]-1.0)/2.0
  z2 = (1-gp[3])*(1+gp[3])
  z3 = gp[3]*(gp[3]+1.0)/2.0
  #
  N[1] = x1*y1*z1 
  N[2] = x3*y1*z1
  N[3] = x1*y3*z1
  N[4] = x3*y3*z1
  #
  N[5] = x1*y1*z3
  N[6] = x3*y1*z3
  N[7] = x1*y3*z3
  N[8] = x3*y3*z3
  #
  N[9] = x2*y1*z1
  N[10] = x1*y2*z1
  N[11] = x2*y2*z1
  N[12] = x3*y2*z1
  #
  N[13] = x2*y3*z1
  N[14] = x1*y1*z2
  N[15] = x2*y1*z2
  N[16] = x3*y1*z2
  #
  N[17] = x1*y2*z2
  N[18] = x2*y2*z2
  N[19] = x3*y2*z2
  N[20] = x1*y3*z2
  #
  N[21] = x2*y3*z2
  N[22] = x3*y3*z2
  N[23] = x2*y1*z3
  N[24] = x1*y2*z3
  N[25] = x2*y2*z3
  N[26] = x3*y2*z3
  N[27] = x2*y3*z3
  #
    return nothing
end




"""
> ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:H27}})
It fills ∂N∂a with the derivatives of the shape functions of a T4 element evaluated at gauss point gp.
- ∂N∂a : shape functions derivatives Array
- gp : gauss point coordinates
- :Type{Val{:H27}} : element type specifier
"""
function ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:H27}})
  #
  x1 = gp[1]*(gp[1]-1.0)/2.0
  x2 = (1-gp[1])*(1+gp[1])
  x3 = gp[1]*(gp[1]+1.0)/2.0
  #
  y1 = gp[2]*(gp[2]-1.0)/2.0
  y2 = (1-gp[2])*(1+gp[2])
  y3 = gp[2]*(gp[2]+1.0)/2.0
  #
  z1 = gp[3]*(gp[3]-1.0)/2.0
  z2 = (1-gp[3])*(1+gp[3])
  z3 = gp[3]*(gp[3]+1.0)/2.0
  #
  dx1 = (2.0*gp[1]-1.0)/2.0
  dx2 = -2.0*gp[1]
  dx3 = (2.0*gp[1]+1.0)/2.0
  #
  dy1 = (2.0*gp[2]-1.0)/2.0
  dy2 = -2.0*gp[2]
  dy3 = (2.0*gp[2]+1.0)/2.0
  #
  dz1 = (2.0*gp[3]-1.0)/2.0
  dz2 = -2.0*gp[3]
  dz3 = (2.0*gp[3]+1.0)/2.0
  # 1
  ∂N∂a[1,1] = dx1*y1*z1 
  ∂N∂a[1,2] = x1*dy1*z1 
  ∂N∂a[1,3] = x1*y1*dz1 
  # 2
  ∂N∂a[2,1] = dx3*y1*z1
  ∂N∂a[2,2] = x3*dy1*z1
  ∂N∂a[2,3] = x3*y1*dz1
  # 3
  ∂N∂a[3,1] = dx1*y3*z1
  ∂N∂a[3,2] = x1*dy3*z1
  ∂N∂a[3,3] = x1*y3*dz1
  # 4
  ∂N∂a[4,1] = dx3*y3*z1
  ∂N∂a[4,2] = x3*dy3*z1
  ∂N∂a[4,3] = x3*y3*dz1
  # 5
  ∂N∂a[5,1] = dx1*y1*z3
  ∂N∂a[5,2] = x1*dy1*z3
  ∂N∂a[5,3] = x1*y1*dz3
  # 6
  ∂N∂a[6,1] = dx3*y1*z3
  ∂N∂a[6,2] = x3*dy1*z3
  ∂N∂a[6,3] = x3*y1*dz3
  # 7
  ∂N∂a[7,1] = dx1*y3*z3
  ∂N∂a[7,2] = x1*dy3*z3
  ∂N∂a[7,3] = x1*y3*dz3
  # 8
  ∂N∂a[8,1] = dx3*y3*z3
  ∂N∂a[8,2] = x3*dy3*z3
  ∂N∂a[8,3] = x3*y3*dz3
  # 9
  ∂N∂a[9,1] = dx2*y1*z1
  ∂N∂a[9,2] = x2*dy1*z1
  ∂N∂a[9,3] = x2*y1*dz1
  # 10
  ∂N∂a[10,1] = dx1*y2*z1
  ∂N∂a[10,2] = x1*dy2*z1
  ∂N∂a[10,3] = x1*y2*dz1
  # 11
  ∂N∂a[11,1] = dx2*y2*z1
  ∂N∂a[11,2] = x2*dy2*z1
  ∂N∂a[11,3] = x2*y2*dz1
  # 12
  ∂N∂a[12,1] = dx3*y2*z1
  ∂N∂a[12,2] = x3*dy2*z1
  ∂N∂a[12,3] = x3*y2*dz1
  # 13
  ∂N∂a[13,1] = dx2*y3*z1
  ∂N∂a[13,2] = x2*dy3*z1
  ∂N∂a[13,3] = x2*y3*dz1
  # 14
  ∂N∂a[14,1] = dx1*y1*z2
  ∂N∂a[14,2] = x1*dy1*z2
  ∂N∂a[14,3] = x1*y1*dz2
  # 15
  ∂N∂a[15,1] = dx2*y1*z2
  ∂N∂a[15,2] = x2*dy1*z2
  ∂N∂a[15,3] = x2*y1*dz2
  # 16
  ∂N∂a[16,1] = dx3*y1*z2
  ∂N∂a[16,2] = x3*dy1*z2
  ∂N∂a[16,3] = x3*y1*dz2
  # 17
  ∂N∂a[17,1] = dx1*y2*z2
  ∂N∂a[17,2] = x1*dy2*z2
  ∂N∂a[17,3] = x1*y2*dz2
  # 18
  ∂N∂a[18,1] = dx2*y2*z2
  ∂N∂a[18,2] = x2*dy2*z2
  ∂N∂a[18,3] = x2*y2*dz2
  # 19
  ∂N∂a[19,1] = dx3*y2*z2
  ∂N∂a[19,2] = x3*dy2*z2
  ∂N∂a[19,3] = x3*y2*dz2
  # 20
  ∂N∂a[20,1] = dx1*y3*z2
  ∂N∂a[20,2] = x1*dy3*z2
  ∂N∂a[20,3] = x1*y3*dz2
  # 21
  ∂N∂a[21,1] = dx2*y3*z2
  ∂N∂a[21,2] = x2*dy3*z2
  ∂N∂a[21,3] = x2*y3*dz2
  # 22
  ∂N∂a[22,1] = dx3*y3*z2
  ∂N∂a[22,2] = x3*dy3*z2
  ∂N∂a[22,3] = x3*y3*dz2
  # 23
  ∂N∂a[23,1] = dx2*y1*z3
  ∂N∂a[23,2] = x2*dy1*z3
  ∂N∂a[23,3] = x2*y1*dz3
  # 24
  ∂N∂a[24,1] = dx1*y2*z3
  ∂N∂a[24,2] = x1*dy2*z3
  ∂N∂a[24,3] = x1*y2*dz3
  # 25
  ∂N∂a[25,1] = dx2*y2*z3
  ∂N∂a[25,2] = x2*dy2*z3
  ∂N∂a[25,3] = x2*y2*dz3
  # 26
  ∂N∂a[26,1] = dx3*y2*z3
  ∂N∂a[26,2] = x3*dy2*z3
  ∂N∂a[26,3] = x3*y2*dz3
  # 27
  ∂N∂a[27,1] = dx2*y3*z3
  ∂N∂a[27,2] = x2*dy3*z3
  ∂N∂a[27,3] = x2*y3*dz3
  #
  return nothing
    #
end



"""
> N!(N::Array{Float64},gp::Tuple,::Type{Val{:H20}})
It fills N with the shape functions of a PE6 element evaluated at gauss point gp.
- N : shape functions Array
- gp : gauss point coordinates
- :Type{Val{:H20}} : element type specifier
"""
function N!(N::Array{Float64},gp::Tuple,::Type{Val{:H20}})
  #
  x,y,z = gp
  # [-1,-1,-1] [+1,-1,-1] [+1,+1,-1] [-1,+1,-1]
  N[1]  = (1.0-x) * (1.0-y) * (1.0-z) * (-x-y-z-2.0) / 8.0
  N[4]  = (1.0+x) * (1.0-y) * (1.0-z) * (+x-y-z-2.0) / 8.0
  N[3]  = (1.0+x) * (1.0+y) * (1.0-z) * (+x+y-z-2.0) / 8.0
  N[2]  = (1.0-x) * (1.0+y) * (1.0-z) * (-x+y-z-2.0) / 8.0
  # [-1,-1,+1] [+1,-1,+1] [+1,+1,+1] [-1,+1,+1]
  N[5]  = (1.0-x) * (1.0-y) * (1.0+z) * (-x-y+z-2.0) / 8.0
  N[8]  = (1.0+x) * (1.0-y) * (1.0+z) * (+x-y+z-2.0) / 8.0
  N[7]  = (1.0+x) * (1.0+y) * (1.0+z) * (+x+y+z-2.0) / 8.0
  N[6]  = (1.0-x) * (1.0+y) * (1.0+z) * (-x+y+z-2.0) / 8.0
  # [0,-1,-1] [0,+1,-1] [0,+1,+1] [0,-1,+1]
  N[12]  = (1-x^2) * (1.0-y) * (1.0-z) * 0.25
  N[10] = (1-x^2) * (1.0+y) * (1.0-z) * 0.25
  N[14] = (1-x^2) * (1.0+y) * (1.0+z) * 0.25
  N[16] = (1-x^2) * (1.0-y) * (1.0+z) * 0.25
  # [-1,0,-1] [+1,0,-1] [+1,0,+1] [-1,0,+1]
  N[9] = (1-y^2) * (1.0-x) * (1.0-z) * 0.25
  N[11] = (1-y^2) * (1.0+x) * (1.0-z) * 0.25
  N[15] = (1-y^2) * (1.0+x) * (1.0+z) * 0.25
  N[13] = (1-y^2) * (1.0-x) * (1.0+z) * 0.25
  # [-1,-1,0] [+1,-1,0] [+1,+1,0] [-1,+1,0]
  N[17] = (1-z^2) * (1.0-y) * (1.0-x) * 0.25
  N[20] = (1-z^2) * (1.0+y) * (1.0-x) * 0.25
  N[19] = (1-z^2) * (1.0+y) * (1.0+x) * 0.25
  N[18] = (1-z^2) * (1.0-y) * (1.0+x) * 0.25
  #
  return nothing
end




"""
> ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:H20}})
It fills ∂N∂a with the derivatives of the shape functions of a T4 element evaluated at gauss point gp.
- ∂N∂a : shape functions derivatives Array
- gp : gauss point coordinates
- :Type{Val{:H20}} : element type specifier
"""
function ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:H20}})
  #
  x,y,z = gp
  #
  ∂N∂a[1,1] = 1/8*(-1 + y)*(-1 + z)*(1 + 2*x + y + z)
  ∂N∂a[1,2] = 1/8*(-1 + x)*(-1 + z)*(1 + x + 2*y + z)
  ∂N∂a[1,3] = 1/8*(-1 + x)*(-1 + y)*(1 + x + y + 2*z)
  # 2
  ∂N∂a[4,1] = -(1/8)*(-1 + y)*(-1 + z)*(1 - 2*x + y + z)
  ∂N∂a[4,2] = 1/8*(1 + x)*(-1 + x - 2*y - z)*(-1 + z)
  ∂N∂a[4,3] = 1/8*(1 + x)*(-1 + y)*(-1 + x - y - 2*z)
  # 3
  ∂N∂a[3,1] = -(1/8)*(1 + y)*(-1 + 2*x + y - z)*(-1 + z)
  ∂N∂a[3,2] = -(1/8)*(1 + x)*(-1 + x + 2*y - z)*(-1 + z)
  ∂N∂a[3,3] = -(1/8)*(1 + x)*(1 + y)*(-1 + x + y - 2*z)
  # 4
  ∂N∂a[2,1] = 1/8*(1 + y)*(-1 - 2*x + y - z)*(-1 + z)
  ∂N∂a[2,2] = -(1/8)*(-1 + x)*(-1 + z)*(1 + x - 2*y + z)
  ∂N∂a[2,3] = -(1/8)*(-1 + x)*(1 + y)*(1 + x - y + 2*z)
  # 5
  ∂N∂a[5,1] = -(1/8)*(-1 + y)*(1 + 2*x + y - z)*(1 + z)
  ∂N∂a[5,2] = -(1/8)*(-1 + x)*(1 + x + 2*y - z)*(1 + z)
  ∂N∂a[5,3] = -(1/8)*(-1 + x)*(-1 + y)*(1 + x + y - 2*z)
  # 6
  ∂N∂a[8,1] = 1/8*(-1 + y)*(1 - 2*x + y - z)*(1 + z)
  ∂N∂a[8,2] = -(1/8)*(1 + x)*(1 + z)*(-1 + x - 2*y + z)
  ∂N∂a[8,3] = -(1/8)*(1 + x)*(-1 + y)*(-1 + x - y + 2*z)
  # 7
  ∂N∂a[7,1] = 1/8*(1 + y)*(1 + z)*(-1 + 2*x + y + z)
  ∂N∂a[7,2] = 1/8*(1 + x)*(1 + z)*(-1 + x + 2*y + z)
  ∂N∂a[7,3] = 1/8*(1 + x)*(1 + y)*(-1 + x + y + 2*z)
  # 8
  ∂N∂a[6,1] = -(1/8)*(1 + y)*(1 + z)*(-1 - 2*x + y + z)
  ∂N∂a[6,2] = 1/8*(-1 + x)*(1 + x - 2*y - z)*(1 + z)
  ∂N∂a[6,3] = 1/8*(-1 + x)*(1 + y)*(1 + x - y - 2*z)
  # 9
  ∂N∂a[12,1] = -(1/2)*x*(-1 + y)*(-1 + z)
  ∂N∂a[12,2] = -(1/4)*(-1 + x^2)*(-1 + z)
  ∂N∂a[12,3] = -(1/4)*(-1 + x^2)*(-1 + y)
  # 10
  ∂N∂a[10,1] = 1/2*x*(1 + y)*(-1 + z)
  ∂N∂a[10,2] = 1/4*(-1 + x^2)*(-1 + z)
  ∂N∂a[10,3] = 1/4*(-1 + x^2)*(1 + y)
  # 11
  ∂N∂a[14,1] = -(1/2)*x*(1 + y)*(1 + z)
  ∂N∂a[14,2] = -(1/4)*(-1 + x^2)*(1 + z)
  ∂N∂a[14,3] = -(1/4)*(-1 + x^2)*(1 + y)
  # 12
  ∂N∂a[16,1] = 1/2*x*(-1 + y)*(1 + z)
  ∂N∂a[16,2] = 1/4*(-1 + x^2)*(1 + z)
  ∂N∂a[16,3] = 1/4*(-1 + x^2)*(-1 + y)
  # 13
  ∂N∂a[9,1] = -(1/4)*(-1 + y^2)*(-1 + z)
  ∂N∂a[9,2] = -(1/2)*(-1 + x)*y*(-1 + z)
  ∂N∂a[9,3] = -(1/4)*(-1 + x)*(-1 + y^2)
  # 14
  ∂N∂a[11,1] = 1/4*(-1 + y^2)*(-1 + z)
  ∂N∂a[11,2] = 1/2*(1 + x)*y*(-1 + z)
  ∂N∂a[11,3] = 1/4*(1 + x)*(-1 + y^2)
  # 15
  ∂N∂a[15,1] = -(1/4)*(-1 + y^2)*(1 + z)
  ∂N∂a[15,2] = -(1/2)*(1 + x)*y*(1 + z)
  ∂N∂a[15,3] = -(1/4)*(1 + x)*(-1 + y^2)
  # 16
  ∂N∂a[13,1] = 1/4*(-1 + y^2)*(1 + z)
  ∂N∂a[13,2] = 1/2*(-1 + x)*y*(1 + z)
  ∂N∂a[13,3] = 1/4*(-1 + x)*(-1 + y^2)
  # 17
  ∂N∂a[17,1] = -(1/4)*(-1 + y)*(-1 + z^2)
  ∂N∂a[17,2] = -(1/4)*(-1 + x)*(-1 + z^2)
  ∂N∂a[17,3] = -(1/2)*(-1 + x)*(-1 + y)*z
  # 18
  ∂N∂a[20,1] = 1/4*(-1 + y)*(-1 + z^2)
  ∂N∂a[20,2] = 1/4*(1 + x)*(-1 + z^2)
  ∂N∂a[20,3] = 1/2*(1 + x)*(-1 + y)*z
  # 19
  ∂N∂a[19,1] = -(1/4)*(1 + y)*(-1 + z^2)
  ∂N∂a[19,2] = -(1/4)*(1 + x)*(-1 + z^2)
  ∂N∂a[19,3] = -(1/2)*(1 + x)*(1 + y)*z
  # 20
  ∂N∂a[18,1] = 1/4*(1 + y)*(-1 + z^2)
  ∂N∂a[18,2] = 1/4*(-1 + x)*(-1 + z^2)
  ∂N∂a[18,3] = 1/2*(-1 + x)*(1 + y)*z
  #
  return nothing
    #
end





"""
> N!(N::Array{Float64},gp::Tuple,::Type{Val{:P18}})
It fills N with the shape functions of a PE6 element evaluated at gauss point gp.
- N : shape functions Array
- gp : gauss point coordinates
- :Type{Val{:P18}} : element type specifier
"""
function N!(N::Array{Float64},gp::Tuple,::Type{Val{:P18}})
  #
  x,y,z = gp
  a1 = x
  a2 = y
  a3 = 1-x-y
  # 
  N[1] = (2.0*a1-1.0)*a1 * z*(z-1)/2.0
  N[2] = (2.0*a2-1.0)*a2 * z*(z-1)/2.0
  N[3] = (2.0*a3-1.0)*a3 * z*(z-1)/2.0
  #
  N[4] = 4.0*a1*a2       * z*(z-1)/2.0
  N[5] = 4.0*a3*a2       * z*(z-1)/2.0
  N[6] = 4.0*a1*a3       * z*(z-1)/2.0
  #
  N[7] = (2.0*a1-1.0)*a1 * z*(z+1)/2.0
  N[8] = (2.0*a2-1.0)*a2 * z*(z+1)/2.0
  N[9] = (2.0*a3-1.0)*a3 * z*(z+1)/2.0
  #
  N[10] = 4.0*a1*a2       * z*(z+1)/2.0
  N[11] = 4.0*a3*a2       * z*(z+1)/2.0
  N[12] = 4.0*a1*a3       * z*(z+1)/2.0
  #
  N[13] = (2.0*a1-1.0)*a1 * (1+z)*(1-z)
  N[14] = (2.0*a2-1.0)*a2 * (1+z)*(1-z)
  N[15] = (2.0*a3-1.0)*a3 * (1+z)*(1-z)
  #
  N[16] = 4.0*a1*a2       * (1+z)*(1-z)
  N[17] = 4.0*a3*a2       * (1+z)*(1-z)
  N[18] = 4.0*a1*a3       * (1+z)*(1-z)
  #
  return nothing
end



"""
> ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:P18}})
It fills ∂N∂a with the derivatives of the shape functions of a T4 element evaluated at gauss point gp.
- ∂N∂a : shape functions derivatives Array
- gp : gauss point coordinates
- :Type{Val{:P18}} : element type specifier
"""
function ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:P18}})
  #
  x,y,z = gp
  #
  ∂N∂a[1,1] = 1/2.0*(-1.0 + 4.0* x)*(-1.0 + z)*z
  ∂N∂a[1,2] = 0.0
  ∂N∂a[1,3] = 1/2*x*(-1 + 2*x)*(-1 + 2*z)
  #
  ∂N∂a[2,1] = 0.0
  ∂N∂a[2,2] = 1/2*(-1 + 4*y)*(-1 + z)*z
  ∂N∂a[2,3] = 1/2*y*(-1 + 2*y)*(-1 + 2*z)
  #
  ∂N∂a[3,1] = 1/2*(-3 + 4*x + 4*y)*(-1 + z)*z
  ∂N∂a[3,2] = 1/2*(-3 + 4*x + 4*y)*(-1 + z)*z
  ∂N∂a[3,3] = 1/2*(1 + 2*x^2 - 3*y + 2*y^2 + x*(-3 + 4*y))*(-1 + 2*z)
  #
  ∂N∂a[4,1] = 2*y*(-1 + z)*z
  ∂N∂a[4,2] = 2*x*(-1 + z)*z
  ∂N∂a[4,3] = 2*x*y*(-1 + 2*z)
  #
  ∂N∂a[5,1] = -2*y*(-1 + z)*z
  ∂N∂a[5,2] = -2*(-1 + x + 2*y)*(-1 + z)*z
  ∂N∂a[5,3] = -2*y*(-1 + x + y)*(-1 + 2*z)
  #
  ∂N∂a[6,1] = -2*(-1 + 2*x + y)*(-1 + z)*z
  ∂N∂a[6,2] = -2*x*(-1 + z)*z
  ∂N∂a[6,3] = -2*x*(-1 + x + y)*(-1 + 2*z)
  #
  ∂N∂a[7,1] = 1/2*(-1 + 4*x)*z*(1 + z)
  ∂N∂a[7,2] = 0.0
  ∂N∂a[7,3] = 1/2*x*(-1 + 2*x)*(1 + 2*z)
  #
  ∂N∂a[8,1] = 0.0
  ∂N∂a[8,2] = 1/2*(-1 + 4*y)*z*(1 + z)
  ∂N∂a[8,3] = 1/2*y*(-1 + 2*y)*(1 + 2*z)
  #
  ∂N∂a[9,1] = 1/2*(-3 + 4*x + 4*y)*z*(1 + z)
  ∂N∂a[9,2] = 1/2*(-3 + 4*x + 4*y)*z*(1 + z)
  ∂N∂a[9,3] = 1/2*(1 + 2*x^2 - 3*y + 2*y^2 + x*(-3 + 4*y))*(1 + 2*z)
  #
  ∂N∂a[10,1] = 2*y*z*(1 + z)
  ∂N∂a[10,2] = 2*x*z*(1 + z)
  ∂N∂a[10,3] = 2*x*y*(1 + 2*z)
  #
  ∂N∂a[11,1] = -2*y*z*(1 + z)
  ∂N∂a[11,2] = -2*(-1 + x + 2*y)*z*(1 + z)
  ∂N∂a[11,3] = -2*y*(-1 + x + y)*(1 + 2*z)
  #
  ∂N∂a[12,1] = -2*(-1 + 2*x + y)*z*(1 + z)
  ∂N∂a[12,2] = -2*x*z*(1 + z)
  ∂N∂a[12,3] = -2*x*(-1 + x + y)*(1 + 2*z)
  #
  ∂N∂a[13,1] = -((-1 + 4*x)*(-1 + z^2))
  ∂N∂a[13,2] = 0.0
  ∂N∂a[13,3] = 2* (1 - 2*x)*x*z
  #
  ∂N∂a[14,1] = 0.0
  ∂N∂a[14,2] = -((-1 + 4*y)*(-1 + z^2))
  ∂N∂a[14,3] = 2*(1 - 2*y)*y*z
  #
  ∂N∂a[15,1] = -((-3 + 4*x + 4*y)*(-1 + z^2))
  ∂N∂a[15,2] = -((-3 + 4*x + 4*y)*(-1 + z^2))
  ∂N∂a[15,3] = -2*(1 + 2*x^2 - 3*y + 2*y^2 + x*(-3 + 4*y))*z
  #
  ∂N∂a[16,1] = -4*y*(-1 + z)*(1 + z)
  ∂N∂a[16,2] = -4*x*(-1 + z)*(1 + z)
  ∂N∂a[16,3] = -8*x*y*z
  #
  ∂N∂a[17,1] = 4*y*(-1 + z)*(1 + z)
  ∂N∂a[17,2] = 4*(-1 + x + 2*y)*(-1 + z^2)
  ∂N∂a[17,3] = 8*y*(-1 + x + y)*z
  #
  ∂N∂a[18,1] = 4*(-1 + 2*x + y)*(-1 + z^2)
  ∂N∂a[18,2] = 4*x*(-1 + z)*(1 + z)
  ∂N∂a[18,3] = 8*x*(-1 + x + y)*z
  #
  return nothing
    #
end






"""
> N!(N::Array{Float64},gp::Tuple,::Type{Val{:H27}})
It fills N with the shape functions of a PE6 element evaluated at gauss point gp.
- N : shape functions Array
- gp : gauss point coordinates
- :Type{Val{:H27}} : element type specifier
"""
function N27!(N::Array{Float64},gp::Tuple)
  #
  x1 = gp[1]*(gp[1]-1.0)/2.0
  x2 = (1-gp[1])*(1+gp[1])
  x3 = gp[1]*(gp[1]+1.0)/2.0
  #
  y1 = gp[2]*(gp[2]-1.0)/2.0
  y2 = (1-gp[2])*(1+gp[2])
  y3 = gp[2]*(gp[2]+1.0)/2.0
  #
  z1 = gp[3]*(gp[3]-1.0)/2.0
  z2 = (1-gp[3])*(1+gp[3])
  z3 = gp[3]*(gp[3]+1.0)/2.0
  #
  N[1] = x1*y1*z1 
  N[2] = x3*y1*z1
  N[3] = x1*y3*z1
  N[4] = x3*y3*z1
  #
  N[5] = x1*y1*z3
  N[6] = x3*y1*z3
  N[7] = x1*y3*z3
  N[8] = x3*y3*z3
  #
  N[9] = x2*y1*z1
  N[10] = x1*y2*z1
  N[11] = x2*y2*z1
  N[12] = x3*y2*z1
  #
  N[13] = x2*y3*z1
  N[14] = x1*y1*z2
  N[15] = x2*y1*z2
  N[16] = x3*y1*z2
  #
  N[17] = x1*y2*z2
  N[18] = x2*y2*z2
  N[19] = x3*y2*z2
  N[20] = x1*y3*z2
  #
  N[21] = x2*y3*z2
  N[22] = x3*y3*z2
  N[23] = x2*y1*z3
  N[24] = x1*y2*z3
  N[25] = x2*y2*z3
  N[26] = x3*y2*z3
  N[27] = x2*y3*z3
  #
    return nothing
end




"""
> ∂N∂a!(∂N∂a::Array{Float64},gp::Tuple,::Type{Val{:H27}})
It fills ∂N∂a with the derivatives of the shape functions of a T4 element evaluated at gauss point gp.
- ∂N∂a : shape functions derivatives Array
- gp : gauss point coordinates
- :Type{Val{:H27}} : element type specifier
"""
function ∂N∂a27!(∂N∂a::Array{Float64},gp::Tuple)
  #
  x1 = gp[1]*(gp[1]-1.0)/2.0
  x2 = (1-gp[1])*(1+gp[1])
  x3 = gp[1]*(gp[1]+1.0)/2.0
  #
  y1 = gp[2]*(gp[2]-1.0)/2.0
  y2 = (1-gp[2])*(1+gp[2])
  y3 = gp[2]*(gp[2]+1.0)/2.0
  #
  z1 = gp[3]*(gp[3]-1.0)/2.0
  z2 = (1-gp[3])*(1+gp[3])
  z3 = gp[3]*(gp[3]+1.0)/2.0
  #
  dx1 = (2.0*gp[1]-1.0)/2.0
  dx2 = -2.0*gp[1]
  dx3 = (2.0*gp[1]+1.0)/2.0
  #
  dy1 = (2.0*gp[2]-1.0)/2.0
  dy2 = -2.0*gp[2]
  dy3 = (2.0*gp[2]+1.0)/2.0
  #
  dz1 = (2.0*gp[3]-1.0)/2.0
  dz2 = -2.0*gp[3]
  dz3 = (2.0*gp[3]+1.0)/2.0
  # 1
  ∂N∂a[1,1] = dx1*y1*z1 
  ∂N∂a[1,2] = x1*dy1*z1 
  ∂N∂a[1,3] = x1*y1*dz1 
  # 2
  ∂N∂a[2,1] = dx3*y1*z1
  ∂N∂a[2,2] = x3*dy1*z1
  ∂N∂a[2,3] = x3*y1*dz1
  # 3
  ∂N∂a[3,1] = dx1*y3*z1
  ∂N∂a[3,2] = x1*dy3*z1
  ∂N∂a[3,3] = x1*y3*dz1
  # 4
  ∂N∂a[4,1] = dx3*y3*z1
  ∂N∂a[4,2] = x3*dy3*z1
  ∂N∂a[4,3] = x3*y3*dz1
  # 5
  ∂N∂a[5,1] = dx1*y1*z3
  ∂N∂a[5,2] = x1*dy1*z3
  ∂N∂a[5,3] = x1*y1*dz3
  # 6
  ∂N∂a[6,1] = dx3*y1*z3
  ∂N∂a[6,2] = x3*dy1*z3
  ∂N∂a[6,3] = x3*y1*dz3
  # 7
  ∂N∂a[7,1] = dx1*y3*z3
  ∂N∂a[7,2] = x1*dy3*z3
  ∂N∂a[7,3] = x1*y3*dz3
  # 8
  ∂N∂a[8,1] = dx3*y3*z3
  ∂N∂a[8,2] = x3*dy3*z3
  ∂N∂a[8,3] = x3*y3*dz3
  # 9
  ∂N∂a[9,1] = dx2*y1*z1
  ∂N∂a[9,2] = x2*dy1*z1
  ∂N∂a[9,3] = x2*y1*dz1
  # 10
  ∂N∂a[10,1] = dx1*y2*z1
  ∂N∂a[10,2] = x1*dy2*z1
  ∂N∂a[10,3] = x1*y2*dz1
  # 11
  ∂N∂a[11,1] = dx2*y2*z1
  ∂N∂a[11,2] = x2*dy2*z1
  ∂N∂a[11,3] = x2*y2*dz1
  # 12
  ∂N∂a[12,1] = dx3*y2*z1
  ∂N∂a[12,2] = x3*dy2*z1
  ∂N∂a[12,3] = x3*y2*dz1
  # 13
  ∂N∂a[13,1] = dx2*y3*z1
  ∂N∂a[13,2] = x2*dy3*z1
  ∂N∂a[13,3] = x2*y3*dz1
  # 14
  ∂N∂a[14,1] = dx1*y1*z2
  ∂N∂a[14,2] = x1*dy1*z2
  ∂N∂a[14,3] = x1*y1*dz2
  # 15
  ∂N∂a[15,1] = dx2*y1*z2
  ∂N∂a[15,2] = x2*dy1*z2
  ∂N∂a[15,3] = x2*y1*dz2
  # 16
  ∂N∂a[16,1] = dx3*y1*z2
  ∂N∂a[16,2] = x3*dy1*z2
  ∂N∂a[16,3] = x3*y1*dz2
  # 17
  ∂N∂a[17,1] = dx1*y2*z2
  ∂N∂a[17,2] = x1*dy2*z2
  ∂N∂a[17,3] = x1*y2*dz2
  # 18
  ∂N∂a[18,1] = dx2*y2*z2
  ∂N∂a[18,2] = x2*dy2*z2
  ∂N∂a[18,3] = x2*y2*dz2
  # 19
  ∂N∂a[19,1] = dx3*y2*z2
  ∂N∂a[19,2] = x3*dy2*z2
  ∂N∂a[19,3] = x3*y2*dz2
  # 20
  ∂N∂a[20,1] = dx1*y3*z2
  ∂N∂a[20,2] = x1*dy3*z2
  ∂N∂a[20,3] = x1*y3*dz2
  # 21
  ∂N∂a[21,1] = dx2*y3*z2
  ∂N∂a[21,2] = x2*dy3*z2
  ∂N∂a[21,3] = x2*y3*dz2
  # 22
  ∂N∂a[22,1] = dx3*y3*z2
  ∂N∂a[22,2] = x3*dy3*z2
  ∂N∂a[22,3] = x3*y3*dz2
  # 23
  ∂N∂a[23,1] = dx2*y1*z3
  ∂N∂a[23,2] = x2*dy1*z3
  ∂N∂a[23,3] = x2*y1*dz3
  # 24
  ∂N∂a[24,1] = dx1*y2*z3
  ∂N∂a[24,2] = x1*dy2*z3
  ∂N∂a[24,3] = x1*y2*dz3
  # 25
  ∂N∂a[25,1] = dx2*y2*z3
  ∂N∂a[25,2] = x2*dy2*z3
  ∂N∂a[25,3] = x2*y2*dz3
  # 26
  ∂N∂a[26,1] = dx3*y2*z3
  ∂N∂a[26,2] = x3*dy2*z3
  ∂N∂a[26,3] = x3*y2*dz3
  # 27
  ∂N∂a[27,1] = dx2*y3*z3
  ∂N∂a[27,2] = x2*dy3*z3
  ∂N∂a[27,3] = x2*y3*dz3
  #
  return nothing
    #
end