"""
Overview: elemental subroutines
"""


# (Mₑ::Matrix{Float64},Kₑ::Matrix{Float64},
#                        X::Array{Float64},ρ::Float64,Dᵢⱼₖₗ::Matrix{Float64},
#                        N::Vector{Float64},∂N∂a::Matrix{Float64},∂N∂x::Matrix{Float64},
#                        sym∇::Matrix{Float64},Jac::Matrix{Float64},
#                        Jac⁻¹::Matrix{Float64},nn::Int64,etype::Symbol,qr)



function integrate_H!(Fₑ::Vector{ComplexF64},X::Vector{Float64},
                      Dᵢⱼₖₗ::Matrix{Float64},
                      Ψ₁ₑ,Ψ₂ₑ,Ψ₃ₑ,
                      N::Vector{Float64},
                      ∂N∂a::Matrix{Float64},∂N∂x::Matrix{Float64},
                      Jac::Matrix{Float64},Jac⁻¹::Matrix{Float64},
                      ∇U₁::Matrix{ComplexF64},
                      ∇U₂::Matrix{ComplexF64},
                      ∇U₃::Matrix{ComplexF64},
                      e₁₂::Matrix{ComplexF64},
                      e₂₃::Matrix{ComplexF64},
                      e₁₃::Matrix{ComplexF64},
                      eᵛ₁₂::Vector{ComplexF64},eᵛ₂₃::Vector{ComplexF64},eᵛ₁₃::Vector{ComplexF64},
                      σᵛ₁₂::Vector{ComplexF64},σᵛ₂₃::Vector{ComplexF64},σᵛ₁₃::Vector{ComplexF64},
                      sym∇ⁿˡ₁::Matrix{ComplexF64},
                      sym∇ⁿˡ₂::Matrix{ComplexF64},
                      sym∇ⁿˡ₃::Matrix{ComplexF64},
                      nn::Int64,etype::Symbol,qr)
  #
  for i = 1:nn*dim
    @inbounds Fₑ[i] = 0.0
  end
  #
  for (w,gp) in qr
    #
    N!(N,gp,Val{etype})
    ∂N∂a!(∂N∂a,gp,Val{etype})
    #
    Jac_det = metric!(X,∂N∂a,∂N∂x,nn,Jac,Jac⁻¹)
    #
    fill!(∇U₁,0.0)
    fill!(∇U₂,0.0)
    fill!(∇U₃,0.0)
    #
    @inbounds for j = 1:dim
      for i = 1:dim
        for k = 1:nn
          ∇U₁[i,j] += ∂N∂x[k,j]*Ψ₁ₑ[i+(k-1)*dim]
          ∇U₂[i,j] += ∂N∂x[k,j]*Ψ₂ₑ[i+(k-1)*dim]
          ∇U₃[i,j] += ∂N∂x[k,j]*Ψ₃ₑ[i+(k-1)*dim]
        end
      end
    end
    #
    fill!(e₁₂,0.0)
    fill!(e₂₃,0.0)
    fill!(e₁₃,0.0)
    @inbounds for i = 1:dim
      for j = i:dim 
        for k = 1:dim
          e₁₂[i,j] += 0.5*(∇U₁[k,i]*∇U₂[k,j] + ∇U₂[k,i]*∇U₁[k,j])
          e₂₃[i,j] += 0.5*(∇U₂[k,i]*∇U₃[k,j] + ∇U₃[k,i]*∇U₂[k,j])
          e₁₃[i,j] += 0.5*(∇U₁[k,i]*∇U₃[k,j] + ∇U₃[k,i]*∇U₁[k,j])
        end
      end
    end
    #
    tv_sym!(e₁₂,eᵛ₁₂)
    tv_sym!(e₂₃,eᵛ₂₃)
    tv_sym!(e₁₃,eᵛ₁₃)
    #
    fill!(σᵛ₁₂,0.0)
    fill!(σᵛ₂₃,0.0)
    fill!(σᵛ₁₃,0.0)
    for i = 1:dim*(dim-1)
      for j = 1:dim*(dim-1)
        σᵛ₁₂[i] += Dᵢⱼₖₗ[j,i]*eᵛ₁₂[j]
        σᵛ₂₃[i] += Dᵢⱼₖₗ[j,i]*eᵛ₂₃[j]
        σᵛ₁₃[i] += Dᵢⱼₖₗ[j,i]*eᵛ₁₃[j]
      end
    end
    #
    @inbounds for i = 1:nn
      for j = 1:dim
        sym∇ⁿˡ₁[j+(i-1)*dim,1] = (∂N∂x[i,1]*∇U₁[j,1])
        sym∇ⁿˡ₂[j+(i-1)*dim,1] = (∂N∂x[i,1]*∇U₂[j,1])
        sym∇ⁿˡ₃[j+(i-1)*dim,1] = (∂N∂x[i,1]*∇U₃[j,1])
        sym∇ⁿˡ₁[j+(i-1)*dim,2] = (∂N∂x[i,2]*∇U₁[j,2])
        sym∇ⁿˡ₂[j+(i-1)*dim,2] = (∂N∂x[i,2]*∇U₂[j,2])
        sym∇ⁿˡ₃[j+(i-1)*dim,2] = (∂N∂x[i,2]*∇U₃[j,2])
        sym∇ⁿˡ₁[j+(i-1)*dim,3] = (∂N∂x[i,3]*∇U₁[j,3])
        sym∇ⁿˡ₂[j+(i-1)*dim,3] = (∂N∂x[i,3]*∇U₂[j,3])
        sym∇ⁿˡ₃[j+(i-1)*dim,3] = (∂N∂x[i,3]*∇U₃[j,3])
        sym∇ⁿˡ₁[j+(i-1)*dim,4] = (∂N∂x[i,1]*∇U₁[j,2] + ∂N∂x[i,2]*∇U₁[j,1])
        sym∇ⁿˡ₂[j+(i-1)*dim,4] = (∂N∂x[i,1]*∇U₂[j,2] + ∂N∂x[i,2]*∇U₂[j,1])
        sym∇ⁿˡ₃[j+(i-1)*dim,4] = (∂N∂x[i,1]*∇U₃[j,2] + ∂N∂x[i,2]*∇U₃[j,1])
        sym∇ⁿˡ₁[j+(i-1)*dim,5] = (∂N∂x[i,2]*∇U₁[j,3] + ∂N∂x[i,3]*∇U₁[j,2])
        sym∇ⁿˡ₂[j+(i-1)*dim,5] = (∂N∂x[i,2]*∇U₂[j,3] + ∂N∂x[i,3]*∇U₂[j,2])
        sym∇ⁿˡ₃[j+(i-1)*dim,5] = (∂N∂x[i,2]*∇U₃[j,3] + ∂N∂x[i,3]*∇U₃[j,2])
        sym∇ⁿˡ₁[j+(i-1)*dim,6] = (∂N∂x[i,1]*∇U₁[j,3] + ∂N∂x[i,3]*∇U₁[j,1])
        sym∇ⁿˡ₂[j+(i-1)*dim,6] = (∂N∂x[i,1]*∇U₂[j,3] + ∂N∂x[i,3]*∇U₂[j,1])
        sym∇ⁿˡ₃[j+(i-1)*dim,6] = (∂N∂x[i,1]*∇U₃[j,3] + ∂N∂x[i,3]*∇U₃[j,1])
      end
    end
    #
    @inbounds for j = 1:dim*(dim-1)
      for i = 1:nn*dim
        Fₑ[i] += (1.0/6.0) * ( sym∇ⁿˡ₁[i,j]*σᵛ₂₃[j] + sym∇ⁿˡ₂[i,j]*σᵛ₁₃[j] + sym∇ⁿˡ₃[i,j]*σᵛ₁₂[j] )*w*Jac_det
      end
    end
    #
  end
  #
  return nothing
  #
end


# (Mₑ::Matrix{Float64},Kₑ::Matrix{Float64},
#                        X::Array{Float64},ρ::Float64,Dᵢⱼₖₗ::Matrix{Float64},
#                        N::Vector{Float64},∂N∂a::Matrix{Float64},∂N∂x::Matrix{Float64},
#                        sym∇::Matrix{Float64},Jac::Matrix{Float64},
#                        Jac⁻¹::Matrix{Float64},nn::Int64,etype::Symbol,qr)


function integrate_G!(Fₑ::Vector{ComplexF64},X::Vector{Float64},
                      Dᵢⱼₖₗ::Matrix{Float64},
                      Ψ₁ₑ,Ψ₂ₑ,
                      N::Vector{Float64},∂N∂a::Matrix{Float64},
                      ∂N∂x::Matrix{Float64},Jac::Matrix{Float64},Jac⁻¹::Matrix{Float64},
                      ∇U₁::Matrix{ComplexF64},∇U₂::Matrix{ComplexF64},
                      sym∇U₁::Matrix{ComplexF64},sym∇U₂::Matrix{ComplexF64},e₁₂::Matrix{ComplexF64},
                      εᵛ₁::Vector{ComplexF64},εᵛ₂::Vector{ComplexF64},eᵛ₁₂::Vector{ComplexF64},
                      σᵛ₁::Vector{ComplexF64},σᵛ₂::Vector{ComplexF64},σᵛ₁₂::Vector{ComplexF64},
                      sym∇::Matrix{ComplexF64},sym∇ⁿˡ₁::Matrix{ComplexF64},sym∇ⁿˡ₂::Matrix{ComplexF64},
                      nn::Int64,etype::Symbol,qr)
  #
  for i = 1:nn*dim
    @inbounds Fₑ[i] = 0.0
  end
  #
  for (w,gp) in qr
    #
    N!(N,gp,Val{etype})
    ∂N∂a!(∂N∂a,gp,Val{etype})
    #
    Jac_det = metric!(X,∂N∂a,∂N∂x,nn,Jac,Jac⁻¹)
    #
    fill!(∇U₁,0.0)
    fill!(∇U₂,0.0)
    #
    @inbounds for j = 1:dim
      for i = 1:dim
        for k = 1:nn
          ∇U₁[i,j] += ∂N∂x[k,j]*Ψ₁ₑ[i+(k-1)*dim]
          ∇U₂[i,j] += ∂N∂x[k,j]*Ψ₂ₑ[i+(k-1)*dim]
        end
      end
    end
    #
    fill!(e₁₂,0.0)
    @inbounds for i = 1:dim
      for j = i:dim 
        sym∇U₁[i,j] = 0.5*(∇U₁[i,j]+∇U₁[j,i])
        sym∇U₂[i,j] = 0.5*(∇U₂[i,j]+∇U₂[j,i])
        for k = 1:dim
          e₁₂[i,j] += 0.5*(∇U₁[k,i]*∇U₂[k,j] + ∇U₂[k,i]*∇U₁[k,j])
        end
      end
    end
    #
    tv_sym!(sym∇U₁,εᵛ₁)
    tv_sym!(sym∇U₂,εᵛ₂)
    tv_sym!(e₁₂,eᵛ₁₂)
    #
    fill!(σᵛ₁,0.0)
    fill!(σᵛ₂,0.0)
    fill!(σᵛ₁₂,0.0)
    #
    @inbounds for i = 1:dim*(dim-1)
      for j = 1:dim*(dim-1)
        σᵛ₁[i]  += Dᵢⱼₖₗ[j,i]*εᵛ₁[j]
        σᵛ₂[i]  += Dᵢⱼₖₗ[j,i]*εᵛ₂[j]
        σᵛ₁₂[i] += Dᵢⱼₖₗ[j,i]*eᵛ₁₂[j]
      end
    end
    #
    @inbounds for i = 1:nn
      sym∇[1+(i-1)*dim,1] = ∂N∂x[i,1]
      sym∇[2+(i-1)*dim,2] = ∂N∂x[i,2]
      sym∇[3+(i-1)*dim,3] = ∂N∂x[i,3]
      sym∇[1+(i-1)*dim,4] = ∂N∂x[i,2]
      sym∇[2+(i-1)*dim,4] = ∂N∂x[i,1]
      sym∇[2+(i-1)*dim,5] = ∂N∂x[i,3]
      sym∇[3+(i-1)*dim,5] = ∂N∂x[i,2]
      sym∇[1+(i-1)*dim,6] = ∂N∂x[i,3]
      sym∇[3+(i-1)*dim,6] = ∂N∂x[i,1]
    end
    #
    #fill!(sym∇ⁿˡ₁,0.0)
    #fill!(sym∇ⁿˡ₂,0.0)
    #
    @inbounds for i = 1:nn
      for j = 1:dim
        sym∇ⁿˡ₁[j+(i-1)*dim,1] = (∂N∂x[i,1]*∇U₁[j,1])
        sym∇ⁿˡ₂[j+(i-1)*dim,1] = (∂N∂x[i,1]*∇U₂[j,1])

        sym∇ⁿˡ₁[j+(i-1)*dim,2] = (∂N∂x[i,2]*∇U₁[j,2])
        sym∇ⁿˡ₂[j+(i-1)*dim,2] = (∂N∂x[i,2]*∇U₂[j,2])

        sym∇ⁿˡ₁[j+(i-1)*dim,3] = (∂N∂x[i,3]*∇U₁[j,3])
        sym∇ⁿˡ₂[j+(i-1)*dim,3] = (∂N∂x[i,3]*∇U₂[j,3])

        sym∇ⁿˡ₁[j+(i-1)*dim,4] = (∂N∂x[i,1]*∇U₁[j,2] + ∂N∂x[i,2]*∇U₁[j,1])
        sym∇ⁿˡ₂[j+(i-1)*dim,4] = (∂N∂x[i,1]*∇U₂[j,2] + ∂N∂x[i,2]*∇U₂[j,1])

        sym∇ⁿˡ₁[j+(i-1)*dim,5] = (∂N∂x[i,2]*∇U₁[j,3] + ∂N∂x[i,3]*∇U₁[j,2])
        sym∇ⁿˡ₂[j+(i-1)*dim,5] = (∂N∂x[i,2]*∇U₂[j,3] + ∂N∂x[i,3]*∇U₂[j,2])

        sym∇ⁿˡ₁[j+(i-1)*dim,6] = (∂N∂x[i,1]*∇U₁[j,3] + ∂N∂x[i,3]*∇U₁[j,1])
        sym∇ⁿˡ₂[j+(i-1)*dim,6] = (∂N∂x[i,1]*∇U₂[j,3] + ∂N∂x[i,3]*∇U₂[j,1])
      end
    end
    #
    @inbounds for j = 1:dim*(dim-1)
      for i = 1:nn*dim
        Fₑ[i] += 0.5 * ( sym∇[i,j]*σᵛ₁₂[j] + sym∇ⁿˡ₁[i,j]*σᵛ₂[j] + sym∇ⁿˡ₂[i,j]*σᵛ₁[j] )*w*Jac_det
      end
    end
    #
  end
  #
  return nothing
  #
end



"""
> integrate_MK!(Mₑ::Matrix{Float64},Kₑ::Matrix{Float64},X::Array{Float64},
ρ::Float64,Dᵢⱼₖₗ::Matrix{Float64},
N::Array{Float64},∂N∂a::Matrix{Float64},∂N∂x::Matrix{Float64},
sym∇::Matrix{Float64},Jac::Matrix{Float64},
Jac⁻¹::Matrix{Float64},nn,etype::Symbol)
It integrates elemental and stiffness matrices
- Mₑ : elemental mass matrix
- Kₑ : elemental stiffness matrix
- X  : nodal coordinates
- ρ : density
- Dᵢⱼₖₗ : elasticity tensor
- N : shape functions
- ∂N∂a : shape functions master derivatives
- ∂N∂x : shape functions sparial derivatives
- sym∇ : 0.5*(∇(⋅)+ ∇ᵀ(⋅))
- Jac : ∂x∂a
- Jac⁻¹: ∂a∂x
- nn : nodes number
- etype : element type
"""
function integrate_MK!(Mₑ::Matrix{Float64},Kₑ::Matrix{Float64},
                       X::Array{Float64},ρ::Float64,Dᵢⱼₖₗ::Matrix{Float64},
                       N::Vector{Float64},∂N∂a::Matrix{Float64},∂N∂x::Matrix{Float64},
                       sym∇::Matrix{Float64},Jac::Matrix{Float64},
                       Jac⁻¹::Matrix{Float64},nn::Int64,etype::Symbol,qr)
  #
  @inbounds for j = 1:nn*dim
    for i = 1:nn*dim
      Mₑ[i,j] = 0.0
      Kₑ[i,j] = 0.0
    end
  end
  #
  for (w,gp) in qr
    #
    N!(N,gp,Val{etype})
    ∂N∂a!(∂N∂a,gp,Val{etype})
    #
    Jac_det = metric!(X,∂N∂a,∂N∂x,nn,Jac,Jac⁻¹)
    #
    @inbounds for i = 1:nn
      sym∇[1+(i-1)*dim,1] = ∂N∂x[i,1]
      sym∇[2+(i-1)*dim,2] = ∂N∂x[i,2]
      sym∇[3+(i-1)*dim,3] = ∂N∂x[i,3]
      sym∇[1+(i-1)*dim,4] = ∂N∂x[i,2]
      sym∇[2+(i-1)*dim,4] = ∂N∂x[i,1]
      sym∇[2+(i-1)*dim,5] = ∂N∂x[i,3]
      sym∇[3+(i-1)*dim,5] = ∂N∂x[i,2]
      sym∇[1+(i-1)*dim,6] = ∂N∂x[i,3]
      sym∇[3+(i-1)*dim,6] = ∂N∂x[i,1]
    end
    #
    γ = w*Jac_det
    #
    for k = 1:dim*(dim-1)
      for l = 1:dim*(dim-1)
        @inbounds c = Dᵢⱼₖₗ[k,l]
        for j = 1:nn*dim
          for i = 1:nn*dim
            @inbounds Kₑ[i,j] += (sym∇[i,k]*c*sym∇[j,l])*γ
          end
        end
      end
    end
    #
    γ *= ρ
    #
    for j = 1:nn
      for i = j:nn
        @inbounds c = (N[i]*N[j])*γ
        for k = 1:dim
          @inbounds Mₑ[k+(i-1)*dim,k+(j-1)*dim] += c
        end
      end
    end
    #
  end
  #
  for j = 1:nn*dim
    for i = j+1:nn*dim
      @inbounds Mₑ[j,i] = Mₑ[i,j]
      @inbounds Kₑ[j,i] = Kₑ[i,j]
    end
  end
  #
  return nothing
end







"""
> metric!(X::Array{Float64},∂N∂a::Matrix{Float64},
          ∂N∂x::Matrix{Float64},nn::Int64,
          Jac::Matrix{Float64},Jac⁻¹::Matrix{Float64})
It computes element metric
- X  : nodal coordinates
- ∂N∂a : shape functions master derivatives
- ∂N∂x : shape functions sparial derivatives
- nn : nodes number
- Jac : ∂x∂a
- Jac⁻¹: ∂a∂x
"""
function metric!(X::Array{Float64},∂N∂a::Matrix{Float64},
                 ∂N∂x::Matrix{Float64},nn::Int64,
                 Jac::Matrix{Float64},Jac⁻¹::Matrix{Float64})
  #
  fill!(Jac,0.0)
  for i = 1:dim
    for j = 1:dim
      for k = 1:nn
        Jac[i,j] += ∂N∂a[k,j]*X[i+(k-1)*dim]
      end
    end
  end
  #
  @inbounds Jac_det =  Jac[1,1]*(Jac[2,2]*Jac[3,3]-Jac[3,2]*Jac[2,3])
  @inbounds Jac_det -= Jac[1,2]*(Jac[2,1]*Jac[3,3]-Jac[3,1]*Jac[2,3])
  @inbounds Jac_det += Jac[1,3]*(Jac[2,1]*Jac[3,2]-Jac[3,1]*Jac[2,2])
  #
  @inbounds Jac⁻¹[1,1] = (Jac[2,2]*Jac[3,3]-Jac[2,3]*Jac[3,2])/Jac_det
  @inbounds Jac⁻¹[2,1] = (Jac[2,3]*Jac[3,1]-Jac[2,1]*Jac[3,3])/Jac_det
  @inbounds Jac⁻¹[3,1] = (Jac[2,1]*Jac[3,2]-Jac[2,2]*Jac[3,1])/Jac_det
  @inbounds Jac⁻¹[1,2] = (Jac[1,3]*Jac[3,2]-Jac[1,2]*Jac[3,3])/Jac_det
  @inbounds Jac⁻¹[2,2] = (Jac[1,1]*Jac[3,3]-Jac[1,3]*Jac[3,1])/Jac_det
  @inbounds Jac⁻¹[3,2] = (Jac[1,2]*Jac[3,1]-Jac[1,1]*Jac[3,2])/Jac_det
  @inbounds Jac⁻¹[1,3] = (Jac[1,2]*Jac[2,3]-Jac[1,3]*Jac[2,2])/Jac_det
  @inbounds Jac⁻¹[2,3] = (Jac[1,3]*Jac[2,1]-Jac[1,1]*Jac[2,3])/Jac_det
  @inbounds Jac⁻¹[3,3] = (Jac[1,1]*Jac[2,2]-Jac[1,2]*Jac[2,1])/Jac_det
  #
  fill!(∂N∂x,0.0)  
    
  for j = 1:dim
    for i = 1:nn
      for k = 1:dim
        ∂N∂x[i,j] += ∂N∂a[i,k]*Jac⁻¹[k,j]
      end
    end
  end
  #
  return Jac_det
end




function tv_sym!(op_in,op_out)
  @inbounds op_out[1] = op_in[1,1]
  @inbounds op_out[2] = op_in[2,2]
  @inbounds op_out[3] = op_in[3,3]
  @inbounds op_out[4] = op_in[1,2]*2.0
  @inbounds op_out[5] = op_in[2,3]*2.0
  @inbounds op_out[6] = op_in[1,3]*2.0
  return nothing
end