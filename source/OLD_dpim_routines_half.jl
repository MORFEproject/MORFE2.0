
function homological_HALF!(Sol::Vector{ComplexF64},Rhs::Vector{ComplexF64},Mat::SparseMatrixCSC{ComplexF64},
                           Λ::Vector{ComplexF64},Cp::Vector{Parametrisation},p::Int64,Apos::Int64,
                           M::SparseMatrixCSC{Float64},K::SparseMatrixCSC{Float64},
                           AY::Matrix{ComplexF64},info::Infostruct)

  
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
#  Rhs[1:info.nK]=Cp[p].rhs[info.nK+1:info.nA,Apos]-((σ+info.α)*M+info.β*K)*Cp[p].Wf[1:info.nK,Apos]
  Rhs[1:info.nK]=Cp[p].rhs[info.nK+1:info.nA,Apos]-M*Cp[p].Wf[info.nK+1:info.nA,Apos]
                 -((σ+info.α)*M+info.β*K)*Cp[p].Wf[1:info.nK,Apos]
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
#      Rhs[info.nK+j] -= transpose(XTA[j,info.nK+1:info.nA])*Cp[p].Wf[1:info.nK,Apos]
    end  
  end  

  println("Solving system2")
  ps = MKLPardisoSolver()
  solve!(ps,Sol,Mat,Rhs)
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





function homological_HALFFULL!(Sol::Vector{ComplexF64},Rhs::Vector{ComplexF64},Mat::SparseMatrixCSC{ComplexF64},
  Λ::Vector{ComplexF64},Cp::Vector{Parametrisation},p::Int64,Apos::Int64,
  M::SparseMatrixCSC{Float64},K::SparseMatrixCSC{Float64},
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
Rhs[1:info.nA]=Cp[p].rhs[:,Apos]
Rhs[1:info.nK]-=M*Cp[p].Wf[1:info.nK,Apos]
Rhs[info.nK+1:info.nA]-=M*Cp[p].Wf[info.nK+1:info.nA,Apos]

fill!(Mat.nzval,0.0)
Mat[1:info.nK,1:info.nK] = σ*M
Mat[info.nK+1:info.nA,info.nK+1:info.nA] = σ*M+(info.α*M+info.β*K)
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

Rhs_HALF = Array{ComplexF64}(undef,info.nK+info.nz)
Sol_HALF = Array{ComplexF64}(undef,info.nK+info.nz)
Mat_HALF=spzeros(ComplexF64,info.nK+info.nz,info.nK+info.nz)

Mat_HALF[1:info.nK+info.nz,1:info.nK] = Mat[info.nK+1:info.nMat,1:info.nK]+σ*Mat[info.nK+1:info.nMat,info.nK+1:info.nA]
Mat_HALF[1:info.nK+info.nz,info.nK+1:info.nK+info.nz]=Mat[info.nK+1:info.nMat,info.nA+1:info.nMat]
for j = 1:info.nz
if resonant_modes[j]  # if resonant
Mat_HALF[1:info.nK+info.nz,info.nK+j]=Mat[info.nK+1:info.nMat,info.nA+j]+Mat[info.nK+1:info.nMat,info.nK+1:info.nA]*Cp[1].W[1:info.nK,j]
end  
end  
Rhs_HALF[1:info.nK+info.nz]=Rhs[info.nK+1:info.nMat]-Mat[info.nK+1:info.nMat,info.nK+1:info.nA]*Cp[p].Wf[1:info.nK,Apos]

println("Solving system2")
ps = MKLPardisoSolver()
solve!(ps,Sol_HALF,Mat_HALF,Rhs_HALF)
#  Sol.=Mat\Rhs
println("end soluzione")

Cp[p].W[1:info.nK,Apos]=Sol_HALF[1:info.nK]
Cp[p].W[info.nK+1:info.nA,Apos]=σ*Sol_HALF[1:info.nK]+Cp[p].Wf[1:info.nK,Apos]
Cp[p].f[1:info.nz,Apos]=Sol_HALF[info.nK+1:info.nK+info.nz]
for d = 1:info.nz
if resonant_modes[d]  # if resonant
Cp[p].W[info.nK+1:info.nA,Apos]+=Cp[p].f[d,Apos]*Cp[1].W[1:info.nK,d]
end  
end


end

