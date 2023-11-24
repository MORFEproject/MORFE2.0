function count_terms_dyn(Cp::Vector{Parametrisation},info::Infostruct)
    # conto quanti sono
    howmany=0
    for p in 1:info.max_order   # for every order
        for i in 1:Cp[p].nc # for every alpha vector
          #corresp=Cp[p].corresp[i] # that is resonant or under the given order NA
          #  if corresp!=0
              howmany=howmany+1
          #  end
        end
    end
    return howmany
end

function store_dyn_and_map(Cp::Vector{Parametrisation},info::Infostruct,howmany::Integer)

mappings=zeros(howmany,mesh.nn,3)
mappings_vel=zeros(howmany,mesh.nn,3)
mappings_modal=zeros(howmany,info.neig)
mappings_modal_vel=zeros(howmany,info.neig)
Avector=zeros(howmany,Cp[1].nc)
fdyn=zeros(howmany,info.nm*2)

U = Field(mesh,dim)
println("Assemblying M K")
colptr,rowval = assembler_dummy_MK(mesh,U)
val=zeros(Float64,length(rowval))
K=SparseMatrixCSC(info.nK,info.nK,colptr,rowval,val)
M = deepcopy(K)
assembler_MK!(mesh,U,K,M)

neig = maximum([info.Φ; info.neig])   
println("Computing eigenvalues")
λ, ϕ = eigs(K,M,nev=neig,which=:SM)
λ = real(λ)
ϕ = real(ϕ)
  #set the mode directions in a predictable way
  for ieig=1:neig
    if sum(ϕ[:,ieig])<0.0
       ϕ[:,ieig]=-ϕ[:,ieig]
    end
  end

index=1
for p in 1:info.max_order   # for every order
  for i in 1:Cp[p].nc # for every alpha vector
    #corresp=Cp[p].corresp[i] # that is resonant or under the given order
    #  if corresp!=0
        Avector[index,:]=Cp[p].Avector[i]
        index=index+1
    #  end    
  end
end

index=1
for p in 1:info.max_order   # for every order
  #write(ofile,"\nOrder "*string(p)*"\n")
  for i in 1:Cp[p].nc # for every alpha vector
    #corresp=Cp[p].corresp[i] # that is resonant or under the given order
    #if corresp!=0
      for j in 1:info.nm
        fdyn[index,j]=2*real(Cp[p].fr[j,i])
      end  
      for j in 1:info.nm
        fdyn[index,j+info.nm]=-2*imag(Cp[p].fr[j,i])
      end
      index=index+1
    #end    
  end
end

# mi segno le mappe
# mi sergno il numero di dofs
index=1
  for p in 1:info.max_order   # for every order
    for i in 1:Cp[p].nc # for every alpha vector
      #corresp=Cp[p].corresp[i] # that is resonant or under the given order
      #if corresp!=0
        for inode=1:mesh.nn # loop nei nodi
            for idof=1:3 # lop nei dof
              dof=U.dof[(inode-1)*3+idof]# estraggo il dof
              if dof>0 # check se dof è positivo
                mappings[index,inode,idof]=real(Cp[p].Wr[dof,i])
                mappings_vel[index,inode,idof]=real(Cp[p].Wr[info.nK+dof,i])
              end
            end
        end
        for imode=1:info.neig # loop nei modi
          #println(size(ϕ[:,imode]))
          #println(size(M))
          #println(size(real.(Cp[p].Wr[:,i])))
          mappings_modal[index,imode]=ϕ[:,imode]'*M*real.(Cp[p].Wr[1:info.nK,i])
          mappings_modal_vel[index,imode]=ϕ[:,imode]'*M*real.(Cp[p].Wr[info.nK+1:2*info.nK,i])
        end
        index=index+1
      #end
    end      
  end
  
  return  mappings,mappings_vel,mappings_modal,mappings_modal_vel,Avector,fdyn


end


function export_maps(mappings::Array{Float64},Avector::Matrix{Float64},info::Infostruct,howmany::Integer)

    ofile = info.output_dir*"/mappings.vtk"
    #
    # scrivo la parte iniziale
    io = open(ofile,"w")
    write(io,"# vtk DataFile Version 3.0\n")
    write(io,"eigenfunctions\n")
    write(io,"ASCII\n")
    write(io,"\n")
    write(io,"DATASET UNSTRUCTURED_GRID\n")
    write(io,"POINTS "*string(mesh.nn)*" float\n")
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
  
    # parto con le mappe calibrate sugli alpha vector
    for imaps=1:howmany
      write(io,"VECTORS mapping_"*replace(string(Int.(Avector[imaps,:])),", "=>"_")[2:end-1]*" float\n")
      for i = 1:mesh.nn
        for j = 1:dim
          write(io,string(mappings[imaps,i,j])*" ")
        end
        write(io,'\n')
      end 
    end
    close(io)
    return nothing
  end
  