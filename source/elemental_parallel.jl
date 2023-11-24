
function integrate_H_par!(Fₑ::Vector{ComplexF64},X::Vector{Float64},
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
    ∇U₁kth=zeros(dim,dim,nn)
    ∇U₂kth=zeros(dim,dim,nn)
    ∇U₃kth=zeros(dim,dim,nn)
#
    @inbounds for j = 1:dim
        for i = 1:dim
            Threads.@threads :dynamic   for k = 1:nn
                                            ∇U₁kth[i,j,k] = ∂N∂x[k,j]*Ψ₁ₑ[i+(k-1)*dim]
                                            ∇U₂kth[i,j,k] = ∂N∂x[k,j]*Ψ₂ₑ[i+(k-1)*dim]
                                            ∇U₃kth[i,j,k] = ∂N∂x[k,j]*Ψ₃ₑ[i+(k-1)*dim]
                                        end
        end
    end
    ∇U₁=sum(∇U₁kth;dims=3)
    ∇U₂=sum(∇U₂kth;dims=3)
    ∇U₃=sum(∇U₃kth;dims=3)
    #=
    @inbounds for j = 1:dim
        for i = 1:dim
            for k = 1:nn
                ∇U₁[i,j] += ∂N∂x[k,j]*Ψ₁ₑ[i+(k-1)*dim]
                ∇U₂[i,j] += ∂N∂x[k,j]*Ψ₂ₑ[i+(k-1)*dim]
                ∇U₃[i,j] += ∂N∂x[k,j]*Ψ₃ₑ[i+(k-1)*dim]
            end
        end
    end
    =#
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
