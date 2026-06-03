
domains_list = [collect(329:488)]   # hex2 file entity indices 329–488 (all 160 volume elements)

# quad2 face elements: file indices 0–3 → code 1–4 (x=0), file 4–7 → code 5–8 (x=1000)
boundaries_list = [collect(1:8)]
constrained_dof  = [[1, 1, 1]]
bc_vals          = [[0.0, 0.0, 0.0]]

materials_list = ["polysilicon"]
density = 2.32e-3
young_modulus = 160e3
poisson_ratio = 0.22
mat = MORFE_newmaterial("polysilicon", density, young_modulus, poisson_ratio)
materials_dict = Dict("polysilicon" => mat)

info.α          = 0.5369754008568333 / 500.0
info.β          = 0.0
info.Φ          = [1]
info.neig        = 10
info.Ffreq       = 1
info.Fmodes      = [1]
info.Fmult       = 0.5 * [5.0]
info.omega_mul   = 1.0        # primary resonance, matches benchmark_ferrite.jl
info.style       = 'c'
info.max_order   = 9        # matches max_degree = 9 in benchmark_ferrite.jl
info.max_orderNA = 9
dirout           = compose_name_output_dir("beam_h27", info)
info.output_dir  = dirout
