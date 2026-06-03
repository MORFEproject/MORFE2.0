using SparseArrays
using ExtendableSparse
using LinearAlgebra
using FEMQuad
using Arpack
using MAT
using Dates

include("./source/defs.jl")
include("./source/param_struct.jl")
include("./source/shape_functions.jl")
include("./source/materials.jl")
include("./source/mesh.jl")
include("./source/field.jl")
include("./source/assembler.jl")
include("./source/assembler_dummy.jl")
include("./source/elemental.jl")
include("./source/utils.jl")
include("./source/dpim.jl")
include("./source/dpim_routines.jl")
include("./source/realification.jl")
include("./source/export_solution.jl")

root = "beam_h27_20x2x2"
mesh_file = root * ".mphtxt"
info_file = root * ".jl"

include("./input/" * info_file)
info.output_dir = compose_name_output_dir(root, info)

println("Reading mesh")
mesh = read_mesh(mesh_file, domains_list, materials_list, materials_dict,
	boundaries_list, constrained_dof, bc_vals)

U = Field(mesh, dim)

info.nm      = length(info.Φ)
info.nz      = 2 * info.nm
info.nzforce = 2
if info.Ffreq == 0
	info.nzforce = 0
end
info.nrom = info.nz + info.nzforce
info.nK   = U.neq
info.nA   = 2 * info.nK
info.nMat = info.nA + info.nz

println("Output directory: ", info.output_dir)
println("nK = ", info.nK, "  nA = ", info.nA, "  nMat = ", info.nMat)
println("max_order = ", info.max_order, "  nrom = ", info.nrom)

@time Cp = dpim(mesh, U, info)

@time realification!(Cp, info)

@time howmany = count_terms_dyn(Cp, info)

ofile = open(info.output_dir * "/outCFULL.txt", "w")
write(ofile, "\nParametrization f\n")
for p in 1:info.max_order
	write(ofile, "\nOrder " * string(p) * "\n")
	for i in 1:Cp[p].nc
		write(ofile, string(Cp[p].Avector[i]) * "    ")
		for j in 1:info.nz
			write(ofile, string(Cp[p].f[j, i]) * "    ")
		end
		write(ofile, "\n")
	end
end
#write(ofile, "\nParametrization W\n")
#for p in 1:info.max_order
#	write(ofile, "\nOrder " * string(p) * "\n")
#	for i in 1:Cp[p].nc
#		write(ofile, "\n" * string(Cp[p].Avector[i]) * "\n")
#		for j in 1:info.nA
#			write(ofile, string(Cp[p].W[j, i]) * "\n")
#		end
#	end
#end
#close(ofile)

#write_rdyn(info, Cp)

println("All done. Results in: ", info.output_dir)
