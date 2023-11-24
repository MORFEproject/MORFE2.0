
domains_list = [[1,2]]
    
boundaries_list = [[1,11]]
constrained_dof = [[1, 1, 1]]
bc_vals = [[0.0, 0.0, 0.0]]

materials_list= ["polysilicon"]   # one mat for each dmomain
density = 2.32e-3;
young_modulus = 160e3;
poisson_ratio = 0.22;
mat=MORFE_newmaterial("polysilicon",density,young_modulus,poisson_ratio)
materials_dict = Dict("polysilicon"=>mat)

# define parameter for the analysis
info.α = 1.0330847062639803/500.0  
info.β = 0.0
info.Φ = [1]   # list of master modes
info.neig=10    # number of computed modes
info.Ffreq=1  # mode number that will give the freq (only one but +iomegat and -iomegat)
#info.Fmodes=[1, 3]  # loading will be prop to sum of these modes  
#info.Fmult=0.9551328114213917^2*0.09*[0.5, 1]  # with these amplitudes
info.Fmodes=[1; 4]  # loading will be prop to sum of these modes  
info.Fmult=[0.03*0.5; 40.0*0.5]  # with these amplitudes; 0.5 --> forcing 1x
#info.Fmodes=[4]  # loading will be prop to sum of these modes  
#info.Fmult=[100.0*0.5]  # with these amplitudes; 0.5 --> forcing 1x
info.omega_mul=1.0 # moltiplicatore della frequenza delle forzante NON ATTIVO
info.style = 'c'
info.max_order = 9
info.max_orderNA = 9
# faccio una directory unica che contiene nel nome l'ordine della parametrizzazione e info
dirout=compose_name_output_dir("arch_hot_Q500_beta_003_beta_40",info)
info.output_dir=dirout
