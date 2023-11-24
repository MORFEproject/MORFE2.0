
struct material 
  #
  ρ::Float64
  Dᵢⱼₖₗ::Matrix{Float64}
  #
end


function MORFE_newmaterial(name::String, ρ::Float64, E::Float64, ν::Float64)
  #
  Dᵢⱼₖₗ = zeros(Float64,(dim*(dim-1),dim*(dim-1)))
  #
  λ = E*ν/((1.0+ν)*(1.0-2.0*ν))
  μ = E/(1+ν)
  #
  Dᵢⱼₖₗ[1,1] = λ+μ
  Dᵢⱼₖₗ[1,2] = λ
  Dᵢⱼₖₗ[2,1] = λ
  Dᵢⱼₖₗ[2,2] = λ+μ
  Dᵢⱼₖₗ[1,3] = λ
  Dᵢⱼₖₗ[3,1] = λ
  Dᵢⱼₖₗ[2,3] = λ
  Dᵢⱼₖₗ[3,2] = λ
  Dᵢⱼₖₗ[3,3] = λ+μ
  Dᵢⱼₖₗ[4,4] = μ/2.0
  Dᵢⱼₖₗ[5,5] = μ/2.0
  Dᵢⱼₖₗ[6,6] = μ/2.0
  #
  mat = material(ρ, Dᵢⱼₖₗ)
  return mat
end


