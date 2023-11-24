
function assembler_dummy_MK(mesh::Grid,U::Field)

  dofs = zeros(Int64,nne_max*dim)
  neq=U.neq
  Mat = SparseMatrixLNK(Int64,neq,neq)

  for iΩ ∈ mesh.Ω
    for set = 1:iΩ.Sen
      nn = iΩ.Senn[set]
      dimg=nn*dim
      skip = iΩ.eskip[set]
      for e = 1:iΩ.ne[set]
        conn = @view iΩ.e2n[skip+1+(e-1)*nn:skip+e*nn]
        dofs!(U,nn,conn,dofs)
        for irow = 1:dimg
          if dofs[irow]>0
            for jcol = 1:dimg
              if dofs[jcol]>0
                Mat[dofs[irow],dofs[jcol]]=1
              end
            end
          end
        end
      end  
    end
  end

  Mat=SparseMatrixCSC(Mat)
  return Mat.colptr,Mat.rowval

end


#=
function assembler_dummy_Mat(mesh::Grid,U::Field,nz::Int64)

  dofs = zeros(Int64,nne_max*dim)
  neq=U.neq
  twoneq=2*neq
  Mat = SparseMatrixLNK(Int64,twoneq+nz,twoneq+nz)

  for iΩ ∈ mesh.Ω
    for set = 1:iΩ.Sen
      nn = iΩ.Senn[set]
      dimg=nn*dim
      skip = iΩ.eskip[set]
      for e = 1:iΩ.ne[set]
        conn = @view iΩ.e2n[skip+1+(e-1)*nn:skip+e*nn]
        dofs!(U,nn,conn,dofs)
        for irow = 1:dimg
          if dofs[irow]>0
            for jcol = 1:dimg
              if dofs[jcol]>0
                Mat[dofs[irow],dofs[jcol]]=1    # M U
                Mat[neq+dofs[irow],neq+dofs[jcol]]=1 # M U
                Mat[dofs[irow],neq+dofs[jcol]]=1   # -MV
                Mat[neq+dofs[irow],dofs[jcol]]=1  # KU
              end
            end
          end
        end
      end  
    end
  end
  
  Mat[twoneq+1:twoneq+nz,1:twoneq].=1
  Mat[1:twoneq,twoneq+1:twoneq+nz].=1
  Mat[twoneq+1:twoneq+nz,twoneq+1:twoneq+nz].=1

  Mat=SparseMatrixCSC(Mat)
  return Mat.colptr,Mat.rowval

end
=#
