

"""
> make_output_dir(info)
It creates the folder that stores the output of the analysis
"""
function compose_name_output_dir(model_name::String,info::Infostruct)

  dirname= "./output/"*model_name
  dirname=dirname*"_o_"*string(info.max_order)
  dirname=dirname*"_eps_"*string(info.max_orderNA)
  dirname=dirname*"_Fmode_"*replace(string(info.Fmodes)[2:end-1], ", "=>"_")
  dirname=dirname*"_"*replace(string(now()), ":"=>"_")[1:end-4]

  return dirname
end



"""
> make_output_dir(info)
It creates the folder that stores the output of the analysis
"""
function make_output_dir(info::Infostruct)

  try
    mkdir(info.output_dir)
  catch
#    println("existing "*dd)
  end

  return nothing
end


"""
> structure_odir(odir::String, id::Int64, tag = 'a')
It creates the sub-folders that stores the output of the analysis
- odir : global output folder path
- id : number of computed step. Autonomous is zero, non-autonomous in increasing order
- tag : 'a' = autonomous, 'na' = non-autonomous
"""
function structure_odir(odir::String, id::Int64, tag = 'a')
  #
  odir_C = odir*"/comb_"*tag*"_"*string(id)
  odir_W = odir*"/mapping_"*tag*"_"*string(id)
  odir_f = odir*"/rdyn_"*tag*"_"*string(id)
  odir_M = odir*"/manifold_"*tag*"_"*string(id)
  try
    mkdir(odir_W)
  catch
#    println("existing "*odir_W)
  end
  try
    mkdir(odir_C)
  catch
#    println("existing "*odir_C)
  end
  try
    mkdir(odir_f)
  catch
#    println("existing "*odir_f)
  end
  try
    mkdir(odir_M)
  catch
#    println("existing "*odir_M)
  end

return odir_C*"/", odir_W*"/", odir_f*"/", odir_M*"/"
end




"""
> export_eig(mesh_file::String, U::field_type, λ::Vector{ComplexF64},ϕ::Matrix{ComplexF64},out_dir::String)
It exports the eigenfunctions in vtk format and eigenfrequencies in .txt file
- mesh_file : name of the mesh file
- U : dummy field variable
- λ : eigenvalues
- ϕ : eigenfunctions
- out_dir : output directory
"""
function export_eig(mesh::Grid, U::Field, 
                    λ::Vector{Float64},ϕ::Matrix{Float64},
                    info::Infostruct, tag::Int64)
  odir = info.output_dir *"/eig"
  #
  try
    mkdir(odir)
  catch
#    println("existing "*odir)
  end
  #
  oeigv = odir*"/eigenfrequencies"*string(tag)*".txt"
  oeigf = odir*"/eigenfunctions"*string(tag)*".vtk"
  #
  io = open(oeigv,"w")
  #
  for i = 1:size(λ)[1]
    write(io,string(sqrt(λ[i])))
    write(io,"\n")
  end
  close(io)
  #
  io = open(oeigf,"w")
  write(io,"# vtk DataFile Version 3.0\n")
  write(io,"eigenfunctions\n")
  write(io,"ASCII\n")
  write(io,"\n")
  write(io,"DATASET UNSTRUCTURED_GRID\n")
  write(io,"POINTS "*string(mesh.nn)*" float\n")
  #
  for n = 1:mesh.nn
    for i = 1:dim
      write(io,string(mesh.n2c[i+(n-1)*dim])*" ")
    end
    write(io,"\n")
  end
  write(io,"\n")
  # compute total cell number and entries
  tcn = 0
  tce = 0
  for i = 1:mesh.nΩ
    for j = 1:mesh.Ω[i].Sen
      nn = mesh.Ω[i].Senn[j]
      tcn += mesh.Ω[i].ne[j]
      tce += mesh.Ω[i].ne[j]*nn
    end
  end
  #
  write(io,"CELLS "*string(tcn)*" "*string(tce+tcn)*"\n")
  #
  for d = 1:mesh.nΩ
    iΩ = mesh.Ω[d]
    for set = 1:iΩ.Sen
      etype = iΩ.Set[set]
      nn = iΩ.Senn[set]
      skip = iΩ.eskip[set]
      #
      print_cell_nodes(io,iΩ,set,nn,skip,Val{etype})
      #
    end
  end
  #
  write(io,"\n")
  write(io,"CELL_TYPES "*string(tcn)*"\n")
  #
  for d = 1:mesh.nΩ
    iΩ = mesh.Ω[d]
    for set = 1:iΩ.Sen
      etype = iΩ.Set[set]
      print_cell_type(io,iΩ,set,Val{etype})

    end
  end
  #
  write(io,"\n")
  write(io,"POINT_DATA "*string(mesh.nn)*"\n")
  #
  for eigm = 1:size(λ)[1]
    write(io,"VECTORS eigenmode"*string(eigm)*" float\n")
    for i = 1:mesh.nn
      for j = 1:dim
        if (U.dof[j+(i-1)*dim]>0)
          write(io,string(ϕ[U.dof[j+(i-1)*dim],eigm])*" ")
        else
          write(io,string(U.val[j+(i-1)*dim])*" ")
        end
      end
      write(io,'\n')
    end
    write(io,'\n')
  end
  #
  close(io)
  #
  return nothing
end







"""
> print_cell_type(io,iΩ,set,::Type{Val{:TE4}})
It writes to the output file the correct cell type assuming .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:TE4}} : four-nodes linear tetrahedron identifier
"""
function print_cell_type(io,iΩ,set,::Type{Val{:TE4}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"10\n")
  end
  #
  return nothing
  #
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:T10}})
It writes to the output file the correct cell type assuming .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:T10}} : ten-nodes quadratic tetrahedron identifier
"""
function print_cell_type(io,iΩ,set,::Type{Val{:T10}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"24\n")
  end
  #
  return nothing
  #
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:HE8}})
It writes to the output file the correct cell type assuming .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:HE8}} : eight-nodes lienar hexahedron identifier
"""
function print_cell_type(io,iΩ,set,::Type{Val{:HE8}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"12\n")
  end
  #
  return nothing
  #
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:PE6}})
It writes to the output file the correct cell type assuming .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:PE6}} : six-nodes lienar prism identifier
"""
function print_cell_type(io,iΩ,set,::Type{Val{:PE6}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"13\n")
  end
  #
  return nothing
  #
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:H27}})
It writes to the output file the correct cell type assuming .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:H27}} : twentyseven-nodes quadratic hexahedron identifier
"""
function print_cell_type(io,iΩ,set,::Type{Val{:H27}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"29\n")
  end
  #
  return nothing
  #
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:H20}})
It writes to the output file the correct cell type assuming .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:H20}} : twenty-nodes quadratic hexahedron identifier
"""
function print_cell_type(io,iΩ,set,::Type{Val{:H20}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"25\n")
  end
  #
  return nothing
  #
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:P15}})
It writes to the output file the correct cell type assuming .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:P15}} : fifteen-nodes quadratic prism identifier
"""
function print_cell_type(io,iΩ,set,::Type{Val{:P15}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"26\n")
  end
  #
  return nothing
  #
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:P18}})
It writes to the output file the correct cell type assuming .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:P18}} : fifteen-nodes quadratic prism identifier
"""
function print_cell_type(io,iΩ,set,::Type{Val{:P18}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"32\n")
  end
  #
  return nothing
  #
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:TE4}})
It writes to the output file the nodes of an element .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:TE4}} : four-nodes linear tetrahedron identifier
"""
function print_cell_nodes(io,iΩ,set,nn,skip,::Type{Val{:TE4}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"4 ")
    write(io,string(iΩ.e2n[skip+1+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+2+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+3+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+4+(e-1)*nn]-1)*" ")
    write(io,"\n")
  end
  #
  return nothing
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:T10}})
It writes to the output file the nodes of an element .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:T10}} : ten-nodes quadratic tetrahedron identifier
"""
function print_cell_nodes(io,iΩ,set,nn,skip,::Type{Val{:T10}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"10 ")
    write(io,string(iΩ.e2n[skip+1+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+2+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+3+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+4+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+5+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+6+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+7+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+8+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+9+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+10+(e-1)*nn]-1)*" ")
    write(io,"\n")
  end
  #
  return nothing
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:HE8}})
It writes to the output file the nodes of an element .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:HE8}} : eight-nodes linear hexahedron identifier
"""
function print_cell_nodes(io,iΩ,set,nn,skip,::Type{Val{:HE8}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"8 ")
    write(io,string(iΩ.e2n[skip+1+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+2+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+3+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+4+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+5+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+6+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+7+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+8+(e-1)*nn]-1)*" ")
    write(io,"\n")
  end
  #
  return nothing
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:H20}})
It writes to the output file the nodes of an element .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:H20}} : twenty-nodes quadratic hexahedron identifier
"""
function print_cell_nodes(io,iΩ,set,nn,skip,::Type{Val{:H20}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"20 ")
    write(io,string(iΩ.e2n[skip+1+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+2+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+3+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+4+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+5+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+6+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+7+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+8+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+9+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+10+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+11+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+12+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+13+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+14+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+15+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+16+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+17+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+18+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+19+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+20+(e-1)*nn]-1)*" ")
    write(io,"\n")
  end
  #
  return nothing
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:H27}})
It writes to the output file the nodes of an element .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:H27}} : twentyseven-nodes quadratic hexahedron identifier
"""
function print_cell_nodes(io,iΩ,set,nn,skip,::Type{Val{:H27}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"27 ")
    write(io,string(iΩ.e2n[skip+1+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+2+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+4+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+3+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+5+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+6+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+8+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+7+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+9+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+12+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+13+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+10+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+23+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+26+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+27+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+24+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+14+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+16+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+22+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+20+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+17+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+19+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+15+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+21+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+11+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+25+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+18+(e-1)*nn]-1)*" ")
    write(io,"\n")
  end
  #
  return nothing
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:PE6}})
It writes to the output file the nodes of an element .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:PE6}} : six-nodes linear prism identifier
"""
function print_cell_nodes(io,iΩ,set,nn,skip,::Type{Val{:PE6}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"6 ")
    write(io,string(iΩ.e2n[skip+1+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+2+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+3+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+4+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+5+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+6+(e-1)*nn]-1)*" ")
    write(io,"\n")
  end
  #
  return nothing
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:P15}})
It writes to the output file the nodes of an element .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:P15}} : fifteen-nodes quadratic prism identifier
"""
function print_cell_nodes(io,iΩ,set,nn,skip,::Type{Val{:P15}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"15 ")
    write(io,string(iΩ.e2n[skip+1+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+2+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+3+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+4+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+5+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+6+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+7+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+8+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+9+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+10+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+11+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+12+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+13+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+14+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+15+(e-1)*nn]-1)*" ")
    write(io,"\n")
  end
  #
  return nothing
end


"""
> print_cell_type(io,iΩ,set,::Type{Val{:P18}})
It writes to the output file the nodes of an element .vtk format.
- io : file handler
- iΩ : reference domain
- set : reference element number
- ::Type{Val{:P18}} : eighteen-nodes quadratic prism identifier
"""
function print_cell_nodes(io,iΩ,set,nn,skip,::Type{Val{:P18}})
  #
  for e = 1:iΩ.ne[set]
    write(io,"18 ")
    write(io,string(iΩ.e2n[skip+1+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+2+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+3+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+4+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+5+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+6+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+7+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+8+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+9+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+10+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+11+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+12+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+13+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+14+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+15+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+16+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+17+(e-1)*nn]-1)*" ")
    write(io,string(iΩ.e2n[skip+18+(e-1)*nn]-1)*" ")
    write(io,"\n")
  end
  #
  return nothing
end