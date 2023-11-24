using Combinatorics


"""
mutable struct Parametrisation
"""
mutable struct Parametrisation
  nc::Int64
  ncindep::Int64
  Avector::Vector{Any}
  corresp::Vector{Int64}
  W::Matrix{ComplexF64}
  f::Matrix{ComplexF64}
  Wr::Matrix{ComplexF64}
  fr::Matrix{ComplexF64}
  rhs::Matrix{ComplexF64}
  Wf::Matrix{ComplexF64}
  Parametrisation() = new()
end

"""
fills parametrization structure
"""
function initParametrisation!(info::Infostruct)

  Cp = [Parametrisation() for i in 0:info.max_order]
  for p in 1:info.max_order
    Cp[p].Avector,Cp[p].corresp,Cp[p].nc,Cp[p].ncindep=indexset(info,p)
    Cp[p].W = zeros(ComplexF64,info.nA,Cp[p].nc)  
    Cp[p].Wr = zeros(ComplexF64,info.nA,Cp[p].nc)  
    Cp[p].f = zeros(ComplexF64,info.nrom,Cp[p].nc)
    Cp[p].Wf = zeros(ComplexF64,info.nA,Cp[p].nc)  
    Cp[p].fr = zeros(ComplexF64,info.nrom,Cp[p].nc)
    Cp[p].rhs = zeros(ComplexF64,info.nA,Cp[p].nc)
  end
  return Cp
  end 


"""
creates multiexponents and looks for conjugates  
"""
function indexset(info::Infostruct,p::Int64)

nz=info.nz
nzf=info.nzforce  
nzhalf=Int(nz/2)
nzfhalf=Int(nzf/2)
ndof=info.nrom

a=collect(multiexponents(ndof,p))
lena=length(a)
corresp=zeros(Int64,lena)
nc=0
#corresp=[i for i in 1:lena] 
for i in 1:lena
  orderna=sum(a[i][nz+1:ndof])
#  println(orderna)
  if orderna <= info.max_orderNA  
    b=[a[i][nzhalf+1:nz]; a[i][1:nzhalf]; a[i][nz+nzfhalf+1:ndof]; a[i][nz+1:nz+nzfhalf]]
    if corresp[i]==0
      nc+=1
      for j in i:lena
        if a[j]==b 
          corresp[j]=-i
          corresp[i]=j
          continue
        end  
      end 
    end  
  end 
end    

return a,corresp,lena,nc
end





