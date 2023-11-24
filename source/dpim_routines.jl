
# directly builds half matrix 
function homological_HALF!(Sol::Vector{ComplexF64},Rhs::Vector{ComplexF64},Mat::SparseMatrixCSC{ComplexF64},
                           Λ::Vector{ComplexF64},Cp::Vector{Parametrisation},p::Int64,Apos::Int64,
                           M::SparseMatrixCSC{Float64},K::SparseMatrixCSC{Float64},C::SparseMatrixCSC{Float64},
                           AY::Matrix{ComplexF64},info::Infostruct)

  
  σ=dot(Cp[p].Avector[Apos],Λ[1:info.nrom])
 
  resonant_modes = zeros(Bool,info.nrom)
  if info.style=='c'
  fill!(resonant_modes,false)
    for i = 1:info.nz
      λ = Λ[i]  
      if (abs(imag(σ-λ))/abs(imag(λ))<=ϵ_tol)
        resonant_modes[i] = true
      end
    end
  elseif(info.style=='g')
    fill!(resonant_modes,true)
  end
  println(Cp[p].Avector[Apos],"  ",resonant_modes)
  
  fill!(Rhs,0.0)
  Rhs[1:info.nK]=Cp[p].rhs[info.nK+1:info.nA,Apos]-M*Cp[p].Wf[info.nK+1:info.nA,Apos]-(σ*M+C)*Cp[p].Wf[1:info.nK,Apos]
  fill!(Mat.nzval,0.0)
  Mat[1:info.nK,1:info.nK] = (σ^2+info.α*σ)*M + (info.β*σ+1)*K
  for d = 1:info.nz
    Mat[info.nK+d,info.nK+d] = 1.0      
    if resonant_modes[d]  # if resonant
      λ = Λ[d]
      Mat[info.nK+d,1:info.nK] = (σ-conj(λ))*AY[1:info.nK,d] 
      Mat[1:info.nK,info.nK+d] = (σ-conj(λ))*AY[1:info.nK,d] 
      if d<=info.nm && resonant_modes[d+info.nm]  
        Mat[info.nK+d,info.nK+d+info.nm] = 1.0 
      end  
      if d>info.nm && resonant_modes[d-info.nm]  
        Mat[info.nK+d,info.nK+d-info.nm] = 1.0 
      end  
      Rhs[info.nK+d]=-transpose(AY[1:info.nK,d])*Cp[p].Wf[1:info.nK,Apos]
    end  
  end  

#  println("Solving system2")
#  ps = MKLPardisoSolver()
#  solve!(ps,Sol,Mat,Rhs)
  Sol.=Mat\Rhs
#  println("end soluzione")
   
  Cp[p].W[1:info.nK,Apos]=Sol[1:info.nK]
  Cp[p].W[info.nK+1:info.nA,Apos]=σ*Sol[1:info.nK]+Cp[p].Wf[1:info.nK,Apos]
  Cp[p].f[1:info.nz,Apos]=Sol[info.nK+1:info.nK+info.nz]
  for d = 1:info.nz
    if resonant_modes[d]  # if resonant
      Cp[p].W[info.nK+1:info.nA,Apos]+=Cp[p].f[d,Apos]*Cp[1].W[1:info.nK,d]
    end  
  end

end


function homological_FULL!(Sol::Vector{ComplexF64},Rhs::Vector{ComplexF64},Mat::SparseMatrixCSC{ComplexF64},
                           Λ::Vector{ComplexF64},Cp::Parametrisation,Apos::Int64,
                           M::SparseMatrixCSC{Float64},K::SparseMatrixCSC{Float64},C::SparseMatrixCSC{Float64},
                           AY::Matrix{ComplexF64},XTA::Matrix{ComplexF64},info::Infostruct)

                            
  σ=dot(Cp.Avector[Apos],Λ[1:info.nrom])
  resonant_modes = zeros(Bool,info.nrom)
  if info.style=='c'
    fill!(resonant_modes,false)
      for i = 1:info.nz
        λ = Λ[i]  
        if (abs(imag(σ-λ))/abs(imag(λ))<=ϵ_tol)
          resonant_modes[i] = true
        end
      end
    elseif(info.style=='g')
      fill!(resonant_modes,true)
    end
  println(Cp.Avector[Apos],"  ",resonant_modes)

  fill!(Rhs,0.0)
  Rhs[1:info.nA]=Cp.rhs[:,Apos]
  Rhs[1:info.nK]-=M*Cp.Wf[1:info.nK,Apos]
  Rhs[info.nK+1:info.nA]-=M*Cp.Wf[info.nK+1:info.nA,Apos]

  fill!(Mat.nzval,0.0)
  Mat[1:info.nK,1:info.nK] = σ*M
  Mat[info.nK+1:info.nA,info.nK+1:info.nA] = σ*M+C
  Mat[1:info.nK,info.nK+1:info.nA] = -M
  Mat[info.nK+1:info.nA,1:info.nK] = K
  for j = 1:info.nz
    if resonant_modes[j]  # if resonant
      Mat[1:info.nA,info.nA+j]=AY[:,j]
      Mat[info.nA+j,1:info.nA]=XTA[j,:]
    else
      Mat[info.nA+j,info.nA+j]=1   
    end  
  end  

#  println("Solving system2")
#  ps = MKLPardisoSolver()
#  solve!(ps,Sol,Mat,Rhs)
  Sol.=Mat\Rhs
#  println("end soluzione")
  
  Cp.W[:,Apos]=Sol[1:info.nA]
  Cp.f[1:info.nz,Apos]=Sol[info.nA+1:info.nA+info.nz]

end


function fillrhsG!(mesh::Grid,U::Field,Cp::Vector{Parametrisation},p::Int64)

  for p1 in 1:p-1, p2 in 1:p-1
    if (p1+p2)==p  
      for i in 1:Cp[p1].nc, j in 1:Cp[p2].nc 
        Avector=Cp[p1].Avector[i]+Cp[p2].Avector[j]   
        pos=findfirst(x->x==Avector,Cp[p].Avector)
        if Cp[p].corresp[pos]>0
          Ψ₁ = Cp[p1].W[:,i]
          Ψ₂ = Cp[p2].W[:,j]
          assembly_G!(Cp[p],pos,Ψ₁,Ψ₂,mesh,U)
        end 
      end
    end
  end 

end


function fillrhsH!(mesh::Grid,U::Field,Cp::Vector{Parametrisation},p::Int64)

  for p1 in 1:p-2, p2 in 1:p-2, p3 in 1:p-2
    if (p1+p2+p3)==p  
      for i in 1:Cp[p1].nc, j in 1:Cp[p2].nc, k in 1:Cp[p3].nc  
        Avector=Cp[p1].Avector[i]+Cp[p2].Avector[j]+Cp[p3].Avector[k]   
        pos=findfirst(x->x==Avector,Cp[p].Avector)
        if Cp[p].corresp[pos]>0
          Ψ₁ = Cp[p1].W[:,i]
          Ψ₂ = Cp[p2].W[:,j]
          Ψ₃ = Cp[p3].W[:,k]
          assembly_H!(Cp[p],pos,Ψ₁,Ψ₂,Ψ₃,mesh,U)
        end
      end
    end
  end 

end


function fillWf!(Cp::Vector{Parametrisation},p::Int64,info::Infostruct)

  for p1 in 2:p-1, p2 in 2:p-1
    if (p1+p2)==p+1  
      for i in 1:Cp[p1].nc
        A1=Cp[p1].Avector[i][:]
        for j in 1:Cp[p2].nc 
          A2=Cp[p2].Avector[j][:]      
          for s in 1:info.nrom
            if A1[s]>0
              Avector=A1+A2
              Avector[s]-=1
              pos=findfirst(x->x==Avector,Cp[p].Avector)
              Cp[p].Wf[1:info.nK,pos]+=A1[s]*Cp[p1].W[1:info.nK,i]*Cp[p2].f[s,j]
              Cp[p].Wf[info.nK+1:info.nA,pos]+=A1[s]*Cp[p1].W[info.nK+1:info.nA,i]*Cp[p2].f[s,j]
            end  
          end    
        end 
      end
    end
  end

end 
    

function fillWfnonaut!(Cp::Vector{Parametrisation},p::Int64,Avector::Vector{Int64},ind_rhs::Int64,info::Infostruct)

  for r in info.nz+1:info.nrom
    if Avector[r]>0
      for s in 1:info.nz
        fs_r = Cp[1].f[s,r]
        if abs(fs_r)>10^(-8)    # only fills if the reduced dyn is nzero
          Av_W=Avector[:]
          Av_W[s]+=1
          Av_W[r]-=1
          ind_W=findfirst(x->x==Av_W,Cp[p].Avector)
          Cp[p].Wf[1:info.nK,ind_rhs]+=Cp[p].W[1:info.nK,ind_W]*fs_r*Av_W[s]
          Cp[p].Wf[info.nK+1:info.nA,ind_rhs]+=Cp[p].W[info.nK+1:info.nA,ind_W]*fs_r*Av_W[s]
        end
      end
    end  
  end

end


#=
# non optimal creation of half matrix
function homological_HALF_expl!(Sol::Vector{ComplexF64},Rhs::Vector{ComplexF64},Mat::SparseMatrixCSC{ComplexF64},
                                Λ::Vector{ComplexF64},Cp::Vector{Parametrisation},p::Int64,Apos::Int64,
                                M::SparseMatrixCSC{Float64},K::SparseMatrixCSC{Float64},C::SparseMatrixCSC{Float64},
                                AY::Matrix{ComplexF64},XTA::Matrix{ComplexF64},info::Infostruct)

   
  σ=dot(Cp[p].Avector[Apos],Λ[1:info.nrom])
  resonant_modes = zeros(Bool,info.nrom)
  fill!(resonant_modes,false)
  for i = 1:info.nz
    λ = Λ[i]
    if (abs(σ-λ)/abs(λ)<=ϵ_tol)
      resonant_modes[i] = true
    end
  end
  println(Cp[p].Avector[Apos],"  ",resonant_modes)

  fill!(Rhs,0.0)
  fill!(Mat.nzval,0.0)
  
  Mat[1:info.nK,1:info.nK] = σ^2*M+σ*C+K
  Rhs[1:info.nK]=Cp[p].rhs[info.nK+1:info.nA,Apos]-M*Cp[p].Wf[info.nK+1:info.nA,Apos]-(σ*M+C)*Cp[p].Wf[1:info.nK,Apos]
#  Rhs[1:info.nK]-=(σ*M+C)*Cp[p].Wf[1:info.nK,Apos]

  for j = 1:info.nz
    if resonant_modes[j]  # if resonant
      Mat[info.nK+j,1:info.nK] = XTA[j,1:info.nK]+σ*XTA[j,info.nK+1:info.nA]   
      Mat[1:info.nK,info.nK+j] = AY[info.nK+1:info.nA,j]+(σ*M+C)*Cp[1].W[1:info.nK,j]
      Rhs[info.nK+j] -= transpose(XTA[j,info.nK+1:info.nA])*Cp[p].Wf[1:info.nK,Apos]
      for k = 1:info.nz
        if resonant_modes[k]  # if resonant
          Mat[info.nK+k,info.nK+j] = transpose(XTA[k,info.nK+1:info.nA])*Cp[1].W[1:info.nK,j]
        end 
      end 
    else
      Mat[info.nK+j,info.nK+j]=1   
    end
  end  

  println("Solving system2")
  ps = MKLPardisoSolver()
  solve!(ps,Sol,Mat,Rhs)
  #  Sol.=Mat\Rhs
  println("end soluzione")

  Cp[p].W[1:info.nK,Apos]=Sol[1:info.nK]
  Cp[p].W[info.nK+1:info.nA,Apos]=σ*Sol[1:info.nK]+Cp[p].Wf[1:info.nK,Apos]
  Cp[p].f[1:info.nz,Apos]=Sol[info.nK+1:info.nK+info.nz]
  for d = 1:info.nz
    if resonant_modes[d]  # if resonant
      Cp[p].W[info.nK+1:info.nA,Apos]+=Cp[p].f[d,Apos]*Cp[1].W[1:info.nK,d]
    end  
  end

end

=#