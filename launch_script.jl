#push!(LOAD_PATH,joinpath(pwd(),"FEM"))
#using FEM
using SparseArrays
using ExtendableSparse
using LinearAlgebra
using FEMQuad
using Arpack
using MAT
using Dates
#using Pardiso
#using KrylovKit

include("./source/defs.jl")
include("./source/param_struct.jl")
include("./source/shape_functions.jl")
include("./source/materials.jl")
include("./source/mesh.jl")
include("./source/field.jl")
include("./source/assembler.jl")
include("./source/assembler_dummy.jl")
include("./source/elemental.jl")
include("./source/utils.jl")
include("./source/dpim.jl")
include("./source/dpim_routines.jl")
include("./source/realification.jl")
include("./source/export_solution.jl")

##################################################
#  PREPROCESSOR
##################################################

#mesh_file = "cantilever_hexa_mini.mphtxt"
#mesh_file = "cantilever_hexa.mphtxt"
#mesh_file = "cantilever_hexa_uno.mphtxt"
#info_file = "cantilever.jl"

#mesh_file = "blade_5.mphtxt"
#info_file = "blade.jl"
#=
mesh_file = "arch_I.mphtxt"
info_file = "arch_I.jl"


mesh_file = "beam.mphtxt"
info_file = "beam.jl"

mesh_file = "cantilever_prism.mphtxt"
info_file = "cantilever_prism.jl"
=#
#=
mesh_file = "arch_hot.mphtxt"
info_file = "arch_hot2.jl"
=#




mesh_file = "beam.mphtxt"
info_file = "beam_damp.jl"


mesh_file = "beam.mphtxt"
info_file = "beam_superharm.jl"


mesh_file = "arch_2_force.mphtxt"
info_file = "arch_2_force.jl"


mesh_file = "arch_superh.mphtxt"
info_file = "arch_superh.jl"

#mesh_file = "beam.mphtxt"
#info_file = "beam_superharm.jl"




# read the data file 
include("./input/"*info_file)
# read the mesh file and initialise the grid data structure
println("Reading mesh")
mesh = read_mesh(mesh_file,domains_list,materials_list,materials_dict,
                 boundaries_list,constrained_dof,bc_vals)

# initialise a dummy field to store dofs ordering and static solutions
U = Field(mesh,dim)

info.nm=length(info.Φ)   # master modes
info.nz=2*info.nm         
info.nzforce=2  # imposes only two nonautonomous
if info.Ffreq==0 
  info.nzforce=0 
end    
info.nrom=info.nz+info.nzforce
info.nK=U.neq   # dim of FEM problem
info.nA=2*info.nK  # dim of first order sys
info.nMat=info.nA+info.nz  # dim of system to be solved


##################################################
#  CORE ANALYSIS
##################################################

@time Cp=dpim(mesh,U,info) # computes parametrization
@time realification!(Cp,info)  # performs realification

##################################################
#  OUTPUT ON FILE
##################################################

# I write daown in a compact way the output

# Alfa vector preparation

@time howmany=count_terms_dyn(Cp,info)
@time mappings,mappings_vel,mappings_modal,mappings_modal_vel,Avector,fdyn=store_dyn_and_map(Cp,info,howmany)

#export_maps(mappings,Avector,info,howmany)

file = matopen(info.output_dir*"/param.mat", "w")
write(file, "ndof",info.nA)
write(file, "max_order",info.max_order)
write(file, "mappings",mappings)
write(file, "mappings_modal",mappings_modal)
write(file, "mappings_vel",mappings)
write(file, "mappings_modal_vel",mappings_modal_vel)
write(file, "Avector",Avector)
write(file, "fdyn",fdyn)
close(file)

println("done!")




##################################################
#  OUTPUT ON FILE
##################################################

# writes complex parametrization on file
ofile = open(info.output_dir*"/outCFULL.txt","w")

write(ofile,"\nParametrization f\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:Cp[p].nc # for every alpha vector
    write(ofile,string(Cp[p].Avector[i])*"\n")
    for j in 1:info.nz
      write(ofile,string(Cp[p].f[j,i])*"\n")
    end  
  end
end

write(ofile,"\nParametrization W\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:Cp[p].nc # for every alpha vector
    write(ofile,"\n"*string(Cp[p].Avector[i])*"\n")
    for j in 1:info.nA
      write(ofile,string(Cp[p].W[j,i])*"\n")
    end  
  end
end
close(ofile)               

# writes real parametrization on file
ofile = open(info.output_dir*"/outRFULL.txt","w")
write(ofile,"\nParametrization fr\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:Cp[p].nc # for every alpha vector
    write(ofile,string(Cp[p].Avector[i])*"\n")
    for j in 1:info.nm
      write(ofile,string(2*real(Cp[p].fr[j,i]))*"\n")
    end  
    for j in 1:info.nm
      write(ofile,string(-2*imag(Cp[p].fr[j,i]))*"\n")
    end  
  end
end

write(ofile,"\nParametrization Wr\n")
for p in 1:info.max_order   # for every order
  write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:Cp[p].nc # for every alpha vector
    write(ofile,"\n"*string(Cp[p].Avector[i])*"\n")
    for j in 1:info.nA
      write(ofile,string(real(Cp[p].Wr[j,i]))*"\n")
    end  
  end
end
close(ofile)               

write_rdyn(info,Cp)

println("done!")
#==#
#==#