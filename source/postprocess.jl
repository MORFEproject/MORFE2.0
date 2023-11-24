"""
Overview: post process subroutines
"""



"""
> save_matcont_rdyn_automatic_veps(rdyn,ndofs,Ω_list,outdir)
It saves the autonomous reduced dynamics in a format suitable for matcont.
- rdyn : reduced dynamics in matcont format style
- ndofs : number of degrees of freedom
- outdir : output directory
"""
function save_matcont_rdyn_automatic_veps(rdyn,ndofs,Ω_list,outdir)
  #
  odir = outdir*"/matcont_automatic/"
  #
  mkdir(odir)
  #
  var_file = open(odir*"MORFEsystem.m","w")
  control_parameters = "mu"
  for i = 1:size(Ω_list)[1]
    control_parameters *= ",beta"*string(i)
  end
  for i = 1:size(Ω_list)[1]
    control_parameters *= ",w"*string(i)
  end

  core_file = "function out = DNF_example
  out{1} = @init;
  out{2} = @fun_eval;
  out{3} = [];
  out{4} = [];
  out{5} = [];
  out{6} = [];
  out{7} = [];
  out{8} = [];
  out{9} = [];
  out{10}= @userf1;
  end
  % --------------------------------------------------------------------------
  function [tspan,y0,options] = init
  end
  % --------------------------------------------------------------------------
  function jac = jacobian(t,x,$control_parameters)
  end
  % --------------------------------------------------------------------------
  function jacp = jacobianp(t,x,$control_parameters)
  end
  % --------------------------------------------------------------------------
  function hess = hessians(t,kmrgd,$control_parameters)
  end
  % --------------------------------------------------------------------------
  function hessp = hessiansp(t,kmrgd,$control_parameters)
  end
  %---------------------------------------------------------------------------
  function tens3  = der3(t,kmrgd,$control_parameters)
  end
  %---------------------------------------------------------------------------
  function tens4  = der4(t,kmrgd,$control_parameters)
  end
  %---------------------------------------------------------------------------
  function tens5  = der5(t,kmrgd,$control_parameters)
  end
  %
  % --------------------------------------------------------------------------
  function dydt = fun_eval(t,x,$control_parameters)\n"
  #
  write(var_file,core_file)

  for i = 1:ndofs+size(Ω_list)[1]*2
    write(var_file,"z"*string(i)*"="*"x("*string(i)*");\n")
  end

  write(var_file,"dydt=[\n")

  for i = 1:ndofs+size(Ω_list)[1]*2
    write(var_file,rdyn[i][6:end]*";\n")
  end

  write(var_file,"];\n")
  write(var_file,"end")

  #
  close(var_file)
  #
  return nothing
  #
end





"""
> save_matcont_rdyn_automatic(rdyn,ndofs,outdir)
It saves the autonomous reduced dynamics in a format suitable for matcont.
- rdyn : reduced dynamics in matcont format style
- ndofs : number of degrees of freedom
- outdir : output directory
"""
function save_matcont_rdyn_automatic(rdyn,ndofs,outdir)
  #
  odir = outdir*"/matcont_automatic/"
  #
  mkdir(odir)
  #
  var_file = open(odir*"MORFEsystem.m","w")
  #
  write(var_file,
  "function out = DNF_example
  out{1} = @init;
  out{2} = @fun_eval;
  out{3} = [];
  out{4} = [];
  out{5} = [];
  out{6} = [];
  out{7} = [];
  out{8} = [];
  out{9} = [];
  out{10}= @userf1;
  end
  % --------------------------------------------------------------------------
  function [tspan,y0,options] = init
  end
  % --------------------------------------------------------------------------
  function jac = jacobian(t,x,mu)
  end
  % --------------------------------------------------------------------------
  function jacp = jacobianp(t,x,mu)
  end
  % --------------------------------------------------------------------------
  function hess = hessians(t,kmrgd,mu)
  end
  % --------------------------------------------------------------------------
  function hessp = hessiansp(t,kmrgd,mu)
  end
  %---------------------------------------------------------------------------
  function tens3  = der3(t,kmrgd,mu)
  end
  %---------------------------------------------------------------------------
  function tens4  = der4(t,kmrgd,mu)
  end
  %---------------------------------------------------------------------------
  function tens5  = der5(t,kmrgd,mu)
  end
  %
  % --------------------------------------------------------------------------
  function dydt = fun_eval(t,x,mu)\n")

  for i = 1:ndofs
    write(var_file,"z"*string(i)*"="*"x("*string(i)*");\n")
  end

  write(var_file,"dydt=[\n")

  for i = 1:ndofs
    write(var_file,rdyn[i][6:end]*";\n")
  end

  write(var_file,"];\n")
  write(var_file,"end")

  #
  close(var_file)
  #
  return nothing
  #
end




"""
> save_matcont_rdyn(rdyn,ndofs,outdir)
It saves the non-autonomous reduced dynamics in a format suitable for matcont GUI usage.
- rdyn : reduced dynamics in matcont format style
- ndofs : number of degrees of freedom
- outdir : output directory
"""
function save_matcont_rdyn_nonautonomous(rdyn,ndofs,outdir,Ω_list)
  #
  odir = outdir*"/matcont/"
  #
  mkdir(odir)
  #
  var_file = open(odir*"variables.txt","w")
  par_file = open(odir*"parameters.txt","w")
  eqn_file = open(odir*"equations.txt","w")
  #
  var = ["" for i = 1:ndofs+size(Ω_list)[1]*2]
  for i = 1:ndofs+size(Ω_list)[1]*2
    var[i] = "z"*string(i)
    write(eqn_file,rdyn[i]*";\n")
  end
  #
  vars = var[1]
  for i = 2:ndofs+size(Ω_list)[1]*2
    vars *= ","*var[i]
  end
  write(var_file,vars)
  #
  vars = "mu"
  for i = 1:size(Ω_list)[1]
    vars *= ",w"*string(i)
    vars *= ",beta"*string(i)
  end
  #
  write(par_file,vars)
  #
  close(var_file)
  close(par_file)
  close(eqn_file)
  #
  return nothing
  #
end


"""
> save_matcont_rdyn(rdyn,ndofs,outdir)
It saves the autonomous reduced dynamics in a format suitable for matcont GUI usage.
- rdyn : reduced dynamics in matcont format style
- ndofs : number of degrees of freedom
- outdir : output directory
"""
function save_matcont_rdyn(rdyn,ndofs,outdir)
  #
  odir = outdir*"/matcont/"
  #
  mkdir(odir)
  #
  var_file = open(odir*"variables.txt","w")
  par_file = open(odir*"parameters.txt","w")
  eqn_file = open(odir*"equations.txt","w")
  #
  var = ["" for i = 1:ndofs]
  for i = 1:ndofs
    var[i] = "z"*string(i)
    write(eqn_file,rdyn[i]*";\n")
  end
  vars = var[1]
  for i = 2:ndofs
    vars *= ","*var[i]
  end
  write(var_file,vars)
  write(par_file,"mu")
  #
  close(var_file)
  close(par_file)
  close(eqn_file)
  #
  return nothing
  #
end


"""
> fill_rdyn_veps!(rdyn,ndofs,fr,Cp,p,ith_Ω)
- rdyn : reduced dynamics
- ndofs : number of degrees of freedom
- fr : realified reduced dynamics
- Cp : parametrisation of order p
- p : order of the asymptotic development
- ith_Ω : integer represeting which excitation frequency you are exciting
"""
function fill_rdyn_veps!(rdyn,ndofs,fr,Cp,p,ith_Ω)
  #
  coeff = zeros(Int64,ndofs)
  var1 = "z"*string(ndofs+2*ith_Ω-1)
  var2 = "z"*string(ndofs+2*ith_Ω)
  β = "beta"*string(ith_Ω)
  #
  for c = 1:Cp.nc
    #
    fill!(coeff,0)
    #
    for d = 1:p
      coeff[Cp.comb[d,c]] += 1
    end
    #
    monomial = ""
    for d = 1:ndofs
      #
      if (coeff[d]!=0)
        monomial *= "*z"*string(d)*"^"*string(coeff[d])
      end
      #
    end

    for d = 1:ndofs
      #
      if (abs(imag(fr[d,c]))>1e-100)
        rdyn[d] *= " + "*string(-2.0*imag(fr[d,c]))*monomial*"*"*β*"*"*var1
      end
      #
      if (abs(real(fr[d,c]))>1e-100)
        rdyn[d] *= " + "*string(+2.0*real(fr[d,c]))*monomial*"*"*β*"*"*var2
      end
      #
    end
    #
  end
  #
  return nothing
end


"""
> append_rdyn_frequency!(rdyn,ndofs,ith_Ω)
It appends auxiliary variables required to recast the non-autonomous reduced dynamics as autonomous
- rdyn : reduced dynamics
- ndofs : number of degrees of freedom
- ith_Ω : integer represeting which excitation frequency you are exciting
"""
function append_rdyn_frequency!(rdyn,ndofs,ith_Ω)
  #
  freq = "w"*string(ith_Ω)
  var1 = "z"*string(ndofs+2*ith_Ω-1)
  var2 = "z"*string(ndofs+2*ith_Ω)
  #
  eq1 = var1*"' = "*var1*" + "*freq*"*"*var2*" - "*var1*"*("*var1*"^2+"*var2*"^2"*")"
  eq2 = var2*"' = "*var2*" - "*freq*"*"*var1*" - "*var2*"*("*var1*"^2+"*var2*"^2"*")"
  #
  push!(rdyn,eq1)
  push!(rdyn,eq2)
  #
  return nothing
end


"""
> init_rdyn(ndofs)
It initializes the real-valued reduced dynamics in a format compatible with MATCONT
- ndofs : number of degrees of freedom
"""
function init_rdyn(ndofs)
  #
  rdyn = ["" for i ∈ 1:ndofs]
  for i = 1:ndofs
    rdyn[i] = "z"*string(i)*"' = "
  end
  # add unfolding parameter to first dofs identity-tangent to a modal displacement
  nm = Int(ndofs/2)
  for i = 1:nm
    rdyn[nm+i] *= "+mu*z"*string(nm+i)
  end
  #
  return rdyn
  #
end


"""
> init_rdyn(ndofs)
It initializes the real-valued reduced dynamics in a format compatible with MATCONT
- ndofs : number of degrees of freedom
"""
function fill_rdyn!(rdyn,ndofs,fr,Cp,p)
  #
  coeff = zeros(Int64,ndofs)
  #
  for c = 1:Cp.nc
    fill!(coeff,0)
    #
    for d = 1:p
      coeff[Cp.comb[d,c]] += 1
    end
    #
    monomial = ""
    for d = 1:ndofs
      if (coeff[d]!=0)
        monomial *= "*z"*string(d)*"^"*string(coeff[d])
      end
    end
    #
    for d = 1:ndofs
      if (abs(real(fr[d,c]))>1e-100)
        rdyn[d] *= " + "*string(real(fr[d,c]))*monomial
      end
    end
    #
  end
  #
  return nothing
  #
end


"""
> mass_normalization!(ϕ,M,neig)
Modes are normalised such that ϕᵢMϕⱼ=δᵢⱼ
- ϕ : eigenmodes
- M : mass matrix
- neig : number of computed eigenvalues
"""
function mass_normalization!(ϕ,M,neig)
  #
  for i = 1:neig
    c = transpose(ϕ[:,i])*M*ϕ[:,i]
    for j = 1:M.data.m
      ϕ[j,i] /= sqrt(c)
    end
  end
  #
  return nothing
end


"""
> make_output_dir(mesh_file::String)
It creates the folder that stores the output of the analysis
- mesh_file : name of the mesh file
"""
function make_output_dir(mesh_file::String)
  #
  dd = string(now())
  dd = mesh_file*"_"*dd
  dd = replace(dd,":"=>"_")
  dd = replace(dd,"."=>"_")
  dd = "./"*dd
  mkdir(dd)
  #
  return dd
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
  mkdir(odir_C)
  mkdir(odir_W)
  mkdir(odir_f)
  mkdir(odir_M)  
  #
  return odir_C*"/", odir_W*"/", odir_f*"/", odir_M*"/"
end


"""
> export_data!(Cp,Wr,fr,p,ndofs,U,M,Φ,neig,odir_C,odir_W,odir_f,odir_M)
It saves the output of a given order of the parametrisation procedure. Only realified quantities are saved.
- Cp : order p of a given parametrisation
- Wr : realified mapping of order p
- fr : realified reduced dynamics of order p
- p : order of the asymptotic expansion
- ndofs : number of degrees of freedom
- U : displacement field
- M : mass matrix
- Φ : eigenmodes list
- neig : number of computed eigenvalues
- odir_C : path to the folder that stores monomials exponents
- odir_W : path to the folder that stores the realified mappings
- odir_f : path to the folder that stores the realified reduced dynamics
- odir_M : path to the folder that stores the mass-wiedghted projected mappings onto the computed modes
"""
function export_data!(Cp,Wr,fr,p,ndofs,U,M,Φ,neig,
                      odir_C,odir_W,odir_f,odir_M)
  #
  Ψ = zeros(Float64,U.nen)
  Υ = zeros(Float64,U.nen)
  rdyn = zeros(Float64,ndofs)
  #
  for i = 1:Cp.nc
    # compute monomial exponents
    coeff = zeros(Int64,ndofs)
    for j = 1:p
      coeff[Cp.comb[j,i]] += 1
    end
    #
    fill_mapping!(Ψ,Υ,Wr,i,U)
    rdyn[:] = real(fr[:,i])
    #
    export_C(odir_C,coeff,ndofs)
    export_W(Ψ,Υ,coeff,odir_W)
    export_f(rdyn,ndofs,odir_f)
    export_M(Wr,i,odir_M,M,Φ,neig)
    #
  end
  #
  return nothing
end


"""
> export_data_veps!(Cp,Wr,fr,p,ndofs,U,M,Φ,neig,odir_C_c, odir_W_c, odir_f_c, odir_M_c, odir_C_s, odir_W_s, odir_f_s, odir_M_s)
It saves the output of a given order of the non-autonomous parametrisation procedure. Only realified quantities are saved.
- Cp : order p of a given parametrisation
- Wr : realified mapping of order p
- fr : realified reduced dynamics of order p
- p : order of the asymptotic expansion
- ndofs : number of degrees of freedom
- U : displacement field
- M : mass matrix
- Φ : eigenmodes list
- neig : number of computed eigenvalues
- odir_C_c : path to the folder that stores monomials exponents, cosine harmonic
- odir_W_c : path to the folder that stores the realified mappings, cosine harmonic
- odir_f_c : path to the folder that stores the realified reduced dynamics, cosine harmonic
- odir_M_c : path to the folder that stores the mass-wiedghted projected mappings onto the computed modes, cosine harmonic
- odir_C_s : path to the folder that stores monomials exponents, sine harmonic
- odir_W_s : path to the folder that stores the realified mappings, sine harmonic
- odir_f_s : path to the folder that stores the realified reduced dynamics, sine harmonic
- odir_M_s : path to the folder that stores the mass-wiedghted projected mappings onto the computed modes, sine harmonic
"""
function export_data_veps!(Cp,Wr,fr,p,ndofs,U,M,Φ,neig,
                           odir_C_c, odir_W_c, odir_f_c, odir_M_c,
                           odir_C_s, odir_W_s, odir_f_s, odir_M_s)
  #
  Ψ = zeros(Float64,U.nen)
  Υ = zeros(Float64,U.nen)
  rdyn = zeros(Float64,ndofs)
  #
  for i = 1:Cp.nc
    # compute monomial exponents
    coeff = zeros(Int64,ndofs)
    for j = 1:p
      coeff[Cp.comb[j,i]] += 1
    end
    #
    #
    fill_mapping_c!(Ψ,Υ,Wr,i,U)
    fill_rdyn_c!(rdyn,fr,i,ndofs)
    #
    export_C(odir_C_c,coeff,ndofs)
    export_W(Ψ,Υ,coeff,odir_W_c)
    export_f(rdyn,ndofs,odir_f_c)
    export_M(Wr,i,odir_M_c,M,Φ,neig,+2.0)
    #
    fill_mapping_s!(Ψ,Υ,Wr,i,U)
    fill_rdyn_s!(rdyn,fr,i,ndofs)
    #
    export_C(odir_C_s,coeff,ndofs)
    export_W(Ψ,Υ,coeff,odir_W_s)
    export_f(rdyn,ndofs,odir_f_s)
    export_M(Wr,i,odir_M_s,M,Φ,neig,-2.0)
    #
  end
  #
  return nothing
end


"""
> export_M(Wr,comb,odir,M,Φ,neig)
it saves the mass-wieighted projections of the mappings onto the computed modal basis.
- Wr : realified mapping of order p
- comb : saved index combination
- odir : output directory
- M : mass-matrix
- Φ : eigenmodes
- neig : number of computed eigenmodes
"""
function export_M(Wr,comb,odir,M,Φ,neig,scale=1.0)
  #
  io1 = open(odir*"/upsilon.txt","a")
  io2 = open(odir*"/psi.txt","a")
  #
  var_Υ = M*Wr[1:M.data.m,comb]
  var_Ψ = M*Wr[M.data.m+1:2*M.data.m,comb]
  #
  for i = 1:neig
    #
    Υ = real(transpose(Φ[:,i])*var_Υ)*scale
    Ψ = real(transpose(Φ[:,i])*var_Ψ)*scale
    #
    write(io1,string(real(Υ))*" ")
    write(io2,string(real(Ψ))*" ")
    #
  end
  #
  write(io1,"\n")
  write(io2,"\n")
  #
  close(io1)
  close(io2)
  #
  return nothing
end


"""
> export_C(odir,coeff,ndofs)
it saves the monomials exponents.
- odir : output directory
- coeff : exponents of the monomial
- ndofs : number of degrees of freedom of the reduced model
"""
function export_C(odir,coeff,ndofs)
  #
  io = open(odir*"/monomials.txt","a")
  #
  for i = 1:ndofs
    write(io,string(coeff[i])*" ")
  end
  #
  write(io,"\n")
  close(io)
  #
  return nothing
end


"""
> export_f(fr,ndofs,odir)
it saves the realified reduced dynamics.
- fr : realified reduced dynamics vector
- ndofs : number of degrees of freedom of the reduced model
- odir : output directory
"""
function export_f(fr,ndofs,odir)
  #
  io = open(odir*"/coefficients.txt","a")
  #
  for i = 1:ndofs
    write(io,string(real(fr[i]))*" ")
  end
  #
  write(io,"\n")
  close(io)
  #
  return nothing
end


"""
> export_W(Ψ,Υ,coeff,U,odir_W)
it saves the realified mappings
- Ψ : realified displacement mapping
- Υ : realified velocity mapping
- odir : output directory
"""
function export_W(Ψ,Υ,coeff,odir_W)
  #
  fname_p = "map_"
  #
  tag = ""
  #
  for i = 1:size(coeff)[1]
    tag *= string(coeff[i])*"_"
  end
  fname_p = odir_W*fname_p*tag*".jld2"
  #
  jldsave(fname_p; Ψ, Υ)
  #
  return nothing
  #
end


"""
> fill_mapping!(Ψ,Υ,Wr,e,U)
it copies the mappings but taking into account the global nodes ordering
- Ψ : realified displacement mapping
- Υ : realified velocity mapping
- Wr : realified mapping
- e : combination index
- U : displacement field
"""
function fill_mapping!(Ψ,Υ,Wr,e,U)
  #
  neq = U.neq
  #
  for i = 1:U.nen
    if (U.dof[i]>0)
      Υ[i] = real(Wr[U.dof[i],e])
      Ψ[i] = real(Wr[U.dof[i]+neq,e])
    end
  end
  #
  return nothing
  #
end


"""
> fill_mapping!(Ψ,Υ,Wr,e,U)
it copies the mappings but taking into account the global nodes ordering. Cosine variant.
- Ψ : realified displacement mapping
- Υ : realified velocity mapping
- Wr : realified mapping
- e : combination index
- U : displacement field
"""
function fill_mapping_c!(Ψ,Υ,Wr,e,U)
  #
  neq = U.neq
  #
  for i = 1:U.nen
    if (U.dof[i]>0)
      Υ[i] = 2.0*real(Wr[U.dof[i],e])
      Ψ[i] = 2.0*real(Wr[U.dof[i]+neq,e])
    end
  end
  #
  return nothing
  #
end


"""
> fill_mapping!(Ψ,Υ,Wr,e,U)
it copies the mappings but taking into account the global nodes ordering. Since variant.
- Ψ : realified displacement mapping
- Υ : realified velocity mapping
- Wr : realified mapping
- e : combination index
- U : displacement field
"""
function fill_mapping_s!(Ψ,Υ,Wr,e,U)
  #
  neq = U.neq
  #
  for i = 1:U.nen
    if (U.dof[i]>0)
      Υ[i] = -2.0*imag(Wr[U.dof[i],e])
      Ψ[i] = -2.0*imag(Wr[U.dof[i]+neq,e])
    end
  end
  #
  return nothing
  #
end


"""
> fill_rdyn_c!(rdyn,fr,i,ndofs)
Harmonics realification step. Cosine variant.
- rdyn : realified reduced dynamics
- fr : reduced dynamics
- i : index combination
- ndofs : number of degrees of freedom
"""
function fill_rdyn_c!(rdyn,fr,i,ndofs)
  #
  for d = 1:ndofs
    rdyn[d] = 2.0*real(fr[d,i])
  end
  #
  return nothing
end


"""
> fill_rdyn_s!(rdyn,fr,i,ndofs)
Harmonics realification step. Sine variant.
- rdyn : realified reduced dynamics
- fr : reduced dynamics
- i : index combination
- ndofs : number of degrees of freedom
"""
function fill_rdyn_s!(rdyn,fr,i,ndofs)
  #
  for d = 1:ndofs
    rdyn[d] = -2.0*imag(fr[d,i])
  end
  #
  return nothing
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
                    out_dir::String, tag::Int64)
  odir = out_dir*"/eig"
  #
  try
    mkdir(odir)
  catch
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