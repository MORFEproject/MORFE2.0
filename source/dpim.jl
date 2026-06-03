

function dpim(mesh::Grid, U::Field, info::Infostruct)

	# initialize output folder
	println("Initializing directories")
	make_output_dir(info)

	#odir_C, odir_W, odir_f, odir_M = structure_odir(odir,0)

	# initialise sparseCSC matrices
	println("Assemblying M K")
	colptr, rowval = assembler_dummy_MK(mesh, U)
	val=zeros(Float64, length(rowval))
	K=SparseMatrixCSC(info.nK, info.nK, colptr, rowval, val)
	M = deepcopy(K)
	assembler_MK!(mesh, U, K, M)
	C = info.α*M+info.β*K
	println("Size of M: ", size(M))

	# check maximum number of computer eigenvalues
	neig = maximum([info.Φ; info.neig])
	println("Computing eigenvalues")
	λ, ϕ = eigs(K, M, nev = neig, which = :SM, maxiter = 3000)
	λ = real(λ)
	ϕ = real(ϕ)
	println(size(ϕ))
	#set the mode directions in a predictable way
	for ieig ∈ 1:neig
		if sum(ϕ[:, ieig])<0.0
			ϕ[:, ieig]=-ϕ[:, ieig]
		end
	end

	# mass-normalise the modes
	mass_normalization!(ϕ, M, neig)
	# export eigenvalues
	export_eig(mesh, U, λ, ϕ, info, 0)

	println("Init Parametrisation")
	Cp=initParametrisation!(info)

	println("Init System1")
	#Rhs_FULL = Array{ComplexF64}(undef, info.nMat)
	#Sol_FULL = Array{ComplexF64}(undef, info.nMat)
	#Mat_FULL=spzeros(ComplexF64, info.nMat, info.nMat)
	Rhs = Array{ComplexF64}(undef, info.nK+info.nz)
	Sol = Array{ComplexF64}(undef, info.nK+info.nz)
	Mat=spzeros(ComplexF64, info.nK+info.nz, info.nK+info.nz)

	println("Order 1")
	ω₀ = zeros(Float64, info.nm)
	ζ₀ = zeros(Float64, info.nm)
	AY = Array{ComplexF64}(undef, info.nA, info.nz)
	XTA = Array{ComplexF64}(undef, info.nz, info.nA)
	for i ∈ 1:info.nm
		ω₀[i] = sqrt(λ[info.Φ[i]])
		ζ₀[i] = 0.5*(info.α/ω₀[i]+info.β*ω₀[i])
		println(ζ₀[i])
		λ₁ = -ζ₀[i]*ω₀[i]+ω₀[i]*sqrt(Complex(1.0-ζ₀[i]^2.0))*im
		λ₂ = -ζ₀[i]*ω₀[i]-ω₀[i]*sqrt(Complex(1.0-ζ₀[i]^2.0))*im

		println(λ₁)
		println(λ₂)

		Cp[1].f[i, i] = λ₁
		Cp[1].f[i+info.nm, i+info.nm] = λ₂
		Cp[1].W[1:info.nK, i] = ϕ[:, info.Φ[i]]
		Cp[1].W[(info.nK+1):info.nA, i] = ϕ[:, info.Φ[i]]*λ₁
		Cp[1].W[1:info.nK, i+info.nm] = ϕ[:, info.Φ[i]]
		Cp[1].W[(info.nK+1):info.nA, i+info.nm] = ϕ[:, info.Φ[i]]*λ₂

		AY[1:info.nK, i] = M*Cp[1].W[1:info.nK, i]
		AY[1:info.nK, i+info.nm] = AY[1:info.nK, i]
		AY[(info.nK+1):info.nA, i] = AY[1:info.nK, i]*λ₁
		AY[(info.nK+1):info.nA, i+info.nm] = AY[1:info.nK, i]*λ₂
		XTA[i, 1:info.nK] = -AY[1:info.nK, i]*λ₂/(λ₁-λ₂)
		XTA[i+info.nm, 1:info.nK] = -AY[1:info.nK, i]*λ₁/(λ₂-λ₁)
		XTA[i, (info.nK+1):info.nA] = AY[1:info.nK, i]/(λ₁-λ₂)
		XTA[i+info.nm, (info.nK+1):info.nA] = AY[1:info.nK, i]/(λ₂-λ₁)
	end

	Λ=Vector{ComplexF64}(undef, info.nrom)
	for i in 1:info.nz
		Λ[i]=Cp[1].f[i, i]
	end
	if info.nzforce==2
		Λ[info.nz+1]=im*ω₀[info.Ffreq]*info.omega_mul
		Cp[1].f[info.nz+1, info.nz+1] = Λ[info.nz+1]
		Λ[info.nz+2]=-im*ω₀[info.Ffreq]*info.omega_mul
		Cp[1].f[info.nz+2, info.nz+2] = Λ[info.nz+2]
	end

	#println(norm(Cp[1].rhs))
	# new first order terms due to forcing
	for i in 1:info.nzforce
		for j in axes(info.Fmodes)[1]
			Cp[1].rhs[(info.nK+1):info.nA, info.nz+i]+=info.Fmult[j]*M*ϕ[:, info.Fmodes[j]]
		end
		#homological_FULL!(Sol_FULL, Rhs_FULL, Mat_FULL, Λ, Cp[1], info.nz+i, M, K, C, AY, XTA, info)
		homological_HALF!(Sol, Rhs, Mat, Λ, Cp, 1, info.nz+i, M, K, C, AY, info)
	end

	println("Higher orders")
	# --- benchmark CSV output ---
	f_order = open(info.output_dir * "/benchmark_per_order.csv", "w")
	f_mono  = open(info.output_dir * "/benchmark_per_monomial.csv", "w")
	write(f_order, "order,n_monomials,n_solve,n_conj," *
				   "fillrhsG_time_s,fillrhsG_alloc_bytes," *
				   "fillrhsH_time_s,fillrhsH_alloc_bytes," *
				   "fillWf_time_s,fillWf_alloc_bytes," *
				   "cumul_time_s,mem_live_bytes\n")
	write(
		f_mono,
		"order,monomial_idx,alpha_vector,corresp," *
		"fillWfnonaut_time_s,fillWfnonaut_alloc_bytes," *
		"homological_time_s,homological_alloc_bytes," *
		"conj_copy_time_s,conj_copy_alloc_bytes," *
		"monomial_total_time_s,monomial_total_alloc_bytes," *
		"cumul_time_s,mem_live_bytes\n",
	)

	t_loop_start = time()

	for p ∈ 2:info.max_order
		println("Order $p")

		resG  = @timed fillrhsG!(mesh, U, Cp, p)
		resH  = @timed fillrhsH!(mesh, U, Cp, p)
		resWf = @timed fillWf!(Cp, p, info)

		n_solve = count(x -> x > 0, Cp[p].corresp)
		n_conj  = count(x -> x < 0, Cp[p].corresp)

		for i in 1:Cp[p].nc
			corresp  = Cp[p].corresp[i]
			alpha    = string(Cp[p].Avector[i])
			t_fwnaut = NaN;
			b_fwnaut = Int64(0)
			t_hom    = NaN;
			b_hom    = Int64(0)
			t_conj   = NaN;
			b_conj   = Int64(0)

			if corresp > 0
				r1 = @timed fillWfnonaut!(Cp, p, Cp[p].Avector[i], i, info)
				t_fwnaut = r1.time;
				b_fwnaut = r1.bytes
				#r2 = @timed homological_FULL!(Sol_FULL, Rhs_FULL, Mat_FULL, Λ, Cp[p], i, M, K, C, AY, XTA, info)
				r2 = @timed homological_HALF!(Sol, Rhs, Mat, Λ, Cp, p, i, M, K, C, AY, info)
				t_hom = r2.time;
				b_hom = r2.bytes
			elseif corresp < 0
				r3 = @timed begin
					Cp[p].W[:, i]                       = conj(Cp[p].W[:, -corresp])
					Cp[p].f[1:info.nm, i]               = conj(Cp[p].f[(info.nm+1):(2*info.nm), -corresp])
					Cp[p].f[(info.nm+1):(2*info.nm), i] = conj(Cp[p].f[1:info.nm, -corresp])
				end
				t_conj = r3.time;
				b_conj = r3.bytes
			end

			t_tot = (isnan(t_fwnaut) ? 0.0 : t_fwnaut) +
					(isnan(t_hom) ? 0.0 : t_hom) +
					(isnan(t_conj) ? 0.0 : t_conj)
			b_tot = b_fwnaut + b_hom + b_conj

			t_cumul  = time() - t_loop_start
			mem_live = Base.gc_live_bytes()
			write(f_mono, "$p,$i,\"$alpha\",$corresp," *
						  "$(t_fwnaut),$(b_fwnaut)," *
						  "$(t_hom),$(b_hom)," *
						  "$(t_conj),$(b_conj)," *
						  "$(t_tot),$(b_tot)," *
						  "$(t_cumul),$(mem_live)\n")
		end

		t_cumul  = time() - t_loop_start
		mem_live = Base.gc_live_bytes()
		write(f_order, "$p,$(Cp[p].nc),$n_solve,$n_conj," *
					   "$(resG.time),$(resG.bytes)," *
					   "$(resH.time),$(resH.bytes)," *
					   "$(resWf.time),$(resWf.bytes)," *
					   "$(t_cumul),$(mem_live)\n")
	end

	close(f_order)
	close(f_mono)

	return Cp

end
