
function realification!(Cp::Vector{Parametrisation},info::Infostruct)

  for p in 1:info.max_order
    for i in 1:Cp[p].nc
      Avec=Cp[p].Avector[i][:]
      Ivec=zeros(Int64,p)
      counter=0
      for j in 1:info.nrom 
        Ivec[counter+1:counter+Avec[j]].=j
        counter+=Avec[j]
      end
      pos=1
      coeff=1.0+0.0im
      Avec.=0
      recursive_C2R!(Ivec,p,pos,i,Avec,coeff,Cp[p],info)
    end  
  end

end  



function recursive_C2R!(Ivec::Vector{Int64},p::Int64,pos::Int64,posinit::Int64,Avec::Vector{Int64},
                             coeff::ComplexF64,Cp::Parametrisation,info::Infostruct)

  nzhalf=Int(info.nz/2)
  nzfhalf=Int(info.nzforce/2)                           
  if pos==(p+1)
    pos1=findfirst(x->x==Avec,Cp.Avector)
    Cp.Wr[:,pos1]+=coeff*Cp.W[:,posinit]
    Cp.fr[:,pos1]+=coeff*Cp.f[:,posinit]
  else
    Avec1=Avec[:]   # new vectors
    Avec2=Avec[:]
    iz=Ivec[pos]    # z var 
   #= 
    if iz <= nzhalf         # ORIGINAL
      coeff1=0.5*coeff  
      Avec1[iz]+=1
      coeff2=0.5*im*coeff
      Avec2[iz+nzhalf]+=1
    elseif iz <= info.nz
      coeff1=0.5*coeff
      Avec1[iz-nzhalf]+=1
      coeff2=-0.5*im*coeff
      Avec2[iz]+=1
    elseif iz <= info.nz+nzfhalf  # to have cos and sin the 1/2 must avoided
      coeff1=coeff
      Avec1[iz]+=1
      coeff2=im*coeff
      Avec2[iz+nzfhalf]+=1
    else    
      coeff1=coeff
      Avec1[iz-nzfhalf]+=1
      coeff2=-im*coeff
      Avec2[iz]+=1
    end     
   =# 

    if iz <= nzhalf
      coeff1=0.5*coeff  
      Avec1[iz]+=1
      coeff2=-0.5*im*coeff
      Avec2[iz+nzhalf]+=1
    elseif iz <= info.nz
      coeff1=0.5*coeff
      Avec1[iz-nzhalf]+=1
      coeff2=0.5*im*coeff
      Avec2[iz]+=1
    elseif iz <= info.nz+nzfhalf  # to have cos and sin the 1/2 must avoided
      coeff1=im*coeff
      Avec1[iz]+=1
      coeff2=coeff
      Avec2[iz+nzfhalf]+=1
    else    
      coeff1=-im*coeff
      Avec1[iz-nzfhalf]+=1
      coeff2=coeff
      Avec2[iz]+=1
    end     

    pos+=1
    recursive_C2R!(Ivec,p,pos,posinit,Avec1,coeff1,Cp,info)
    recursive_C2R!(Ivec,p,pos,posinit,Avec2,coeff2,Cp,info)
  end

end    


"""
"""
function write_rdyn(info::Infostruct,Cp::Vector{Parametrisation})

  rdyn = ["" for i in 1:info.nz]
  for i = 1:info.nz
    rdyn[i] = "a"*string(i)*"' = "
  end
  # add unfolding parameter to first dofs identity-tangent to a modal displacement
  #nm = Int(ndofs/2)
  #for i = 1:nm
  #  rdyn[nm+i] *= "+mu*z"*string(nm+i)
  #end

  for p in 1:info.max_order
    for c = 1:Cp[p].nc
      Avector=Cp[p].Avector[c]   
      monomial = ""
      for d = 1:info.nrom        
        if (Avector[d]!=0)
          monomial *= "*a"*string(d)*"^"*string(Avector[d])
        end
      end
      for j in 1:info.nm
        rcoeff=2*real(Cp[p].fr[j,c])
#        
#        icoeff=2*imag(Cp[p].fr[j,c]) # ORIGINAL
        icoeff=-2*imag(Cp[p].fr[j,c])
#
        if abs(rcoeff)>1e-20
          rdyn[j] *= " + "*string(rcoeff)*monomial
        end
        if abs(icoeff)>1e-20
          rdyn[j+info.nm] *= " + "*string(icoeff)*monomial
        end
      end
    end
  end

  ofile = open("./output/equations.txt","w")
  for i = 1:info.nz
    write(ofile,rdyn[i]*";\n")
  end
  close(ofile)  

end

