"""
Overview: functions used to manage dofs handling
"""


"""
> mutable struct Field
- neq : number of equations
- nen : total number of entries
- cpn : nodal uknowns per node
- dof : dofs ordering
- val : field values
"""
mutable struct Field
  neq::Int64
  nen::Int64
  cpn::Int64
  dof::Vector{Int64}
  val::Vector{Float64}
end


"""
> Field(mesh::Grid, cpn::Int64)
Field constructor
- mesh : Grid 
- cpn : number of nodal unknowns per node
"""
function Field(mesh::Grid, cpn::Int64)
  #
  nen = mesh.nn*cpn
  dof = zeros(Int64,nen)
  val = zeros(Float64,nen)
  # filter nodes without domain
  for i = 1:mesh.nΩ
    iΩ = mesh.Ω[i]
    nn = size(iΩ.e2n)[1]
    for j = 1:nn
      for k = 1:cpn
        dof[k+(iΩ.e2n[j]-1)*cpn] = 1
      end
    end
  end
  # fix Dirichlet boundary conditions
  for i = 1:mesh.nΓ
    iΓ = mesh.Γ[i]
    nn = size(iΓ.e2n)[1]
    bcd = iΓ.bcdofs
    bcv = iΓ.bcvals
    for j = 1:cpn
      if (bcd[j]>0)
        for k = 1:nn
          dof[j+(iΓ.e2n[k]-1)*cpn] = -1
          val[j+(iΓ.e2n[k]-1)*cpn] = bcv[j]
        end
      end
    end
  end
  # count number of equations
  neq = 0
  for i = 1:nen
    if (dof[i]>0)
      neq += 1
      dof[i] = neq
    end
  end
  #
  return Field(neq,nen,cpn,dof,val)
  #
end


"""
> dofs!(ϕ::Field, nn::Int64, nodes::Vector{Int64}, dofs::Vector{nt64})
It copies the dofs associated to nodes list nodes into dofs.
- ϕ : field of interest 
- nn : number of queries
- nodes : list of nodes
- dofs : number of degrees of freedom
"""
function dofs!(ϕ::Field, nn::Int64, nodes::Vector{Int64},
               dofs::Vector{Int64})
  #
  for i = 1:nn
    for j = 1:ϕ.cpn
      @inbounds dofs[j+(i-1)*ϕ.cpn] = ϕ.dof[j+(nodes[i]-1)*ϕ.cpn]
    end
  end
  #
  return nothing
end

"""
> dofs!(ϕ::Field, nn::Int64, nodes::SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, dofs::Vector{nt64})
It copies the dofs associated to nodes list nodes into dofs.
This version accepts SubArrays.
- ϕ : field of interest 
- nn : number of queries
- nodes : list of nodes
- dofs : number of degrees of freedom
"""
function dofs!(ϕ::Field, nn::Int64, nodes::SubArray{Int64, 1, Vector{Int64}, Tuple{UnitRange{Int64}}, true},
               dofs::Vector{Int64})
  #
  for i = 1:nn
    for j = 1:ϕ.cpn
      @inbounds dofs[j+(i-1)*ϕ.cpn] = ϕ.dof[j+(nodes[i]-1)*ϕ.cpn]
    end
  end
  #
  return nothing
end


"""
> dofs_vals!(ϕ::Field, nn::Int64, nodes::Vector{Int64}, dofs::Vector{Int64}, vals::Vector{Float64})
It copies the dofs and values associated to nodes list nodes into dofs and vals.
- ϕ : field of interest 
- nn : number of queries
- nodes : list of nodes
- dofs : number of degrees of freedom
- vals : values of the field
"""
function dofs_vals!(ϕ::Field, nn::Int64, nodes::Vector{Int64}, 
                    dofs::Vector{Int64}, vals::Vector{Float64})
  #
  for i = 1:nn
    for j = 1:ϕ.cpn
      @inbounds dofs[j+(i-1)*ϕ.cpn] = ϕ.dof[j+(nodes[i]-1)*ϕ.cpn]
      @inbounds vals[j+(i-1)*ϕ.cpn] = ϕ.val[j+(nodes[i]-1)*ϕ.cpn]
    end
  end
  #
  return nothing
end


"""
> update_field!(ϕ::Field,U::Vector{Float64})
it overwrites free dofs values using the solution U
- ϕ : field to update
- U : solution vector
"""
function update_field!(ϕ::Field,U::Vector{Float64})
  #
  for i = 1:ϕ.nen
    if (ϕ.dof[i]>0)
      @inbounds ϕ.val[i] = U[ϕ.dof[i]]
    end
  end
  #
  return nothing
end


"""
> increment_field!(ϕ::Field,U::Vector{Float64})
it increments field values. ϕ.val += U
- ϕ : field to update
- U : solution vector
"""
function increment_field!(ϕ::Field,U::Vector{Float64})
  #
  for i = 1:ϕ.nen
    if (ϕ.dof[i]>0)
      @inbounds ϕ.val[i] += U[ϕ.dof[i]]
    end
  end
  #
  return nothing
end