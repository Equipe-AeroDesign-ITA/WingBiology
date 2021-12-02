"""
Test code for aerodynamic coefficients of a flexible, straight wing
"""

@test begin
	U∞ = 20.0

	b = 5.0
	c = 1.0

	NW = 41

	# a beam of square cross-section, with dimensions 0.02 x 0.08
	# in material with 100GPa of elasticity modulus and density 1500
	ν = 0.3
	ρ = 1.5e3
	E = 1e11
	G = E / (2 * (ν + 1.0))

	L = 0.02
	Lx = 0.08

	EIx = E * L ^ 3 * Lx / 12
	EIy = E * Lx ^ 3 * L / 12
	EA = E * L * Lx
	GJ = G * (L ^ 2 + Lx ^ 2) * L * Lx / 12
	m = L * Lx * ρ
	Ixx = Lx ^ 2 * m / 6

	e = 0.1
	xcg = 0.1

	afl = WingBiology.flat_plate

	ys = sin.(
		collect(
			LinRange(
				- π / 2, π / 2, NW
			)
		)
	) .* (b / 2)

	points = zeros(3, NW)
	transpose!(
		points,
		hcat(
			zeros(NW), ys, zeros(NW)
		)
	)

	acft = Aircraft(points)

	for i = 1:(NW - 1)
		add_wing_strip!(
			acft,
			i,
			i + 1,
			c;
			incidence = 5.0,
			afl = afl,
			xCG = xcg,
			e = e,
			m = m,
			EA = EA,
			EIy = EIx,
			EIz = EIy,
			GJ = GJ,
			Ixx = Ixx
		)
	end

	fix = NW ÷ 2 + 1

	q0 = get_state(
		acft, U∞
	)
	q = copy(q0)

	qd, fcon = state_space(acft, q, U∞; Sref = b * c, cref = c, bref = b, fixed_points = [fix])

	# @show fcon.total

	A, eig = Arnoldi_solve!(
		acft, q, U∞; 
		Sref = b * c, bref = b, cref = c, fixed_points = [fix], n_Krylov = 3
	)

	# plot_aircraft(acft, q)

	q .= q0

	@info "Arnoldi solution"

	@time begin q .= q0; Arnoldi_solve!(acft, q, U∞; Sref = b * c, bref = b, cref = c, fixed_points = [fix], n_Krylov = 3) end
	@time begin q .= q0; Arnoldi_solve!(acft, q, U∞; Sref = b * c, bref = b, cref = c, fixed_points = [fix], n_Krylov = 3) end
	@time begin q .= q0; Arnoldi_solve!(acft, q, U∞; Sref = b * c, bref = b, cref = c, fixed_points = [fix], n_Krylov = 3) end
	@time begin q .= q0; Arnoldi_solve!(acft, q, U∞; Sref = b * c, bref = b, cref = c, fixed_points = [fix], n_Krylov = 3) end
	@time begin q .= q0; Arnoldi_solve!(acft, q, U∞; Sref = b * c, bref = b, cref = c, fixed_points = [fix], n_Krylov = 3) end

	qd, fcon = state_space(acft, q, U∞; Sref = b * c, cref = c, bref = b, fixed_points = [fix])

	ϕ, ω = get_eigenmodes(acft; fixed_points = [fix])

	@info "Oscillation mode obtention"

	@time ϕ, ω = get_eigenmodes(acft; fixed_points = [fix])
	@time ϕ, ω = get_eigenmodes(acft; fixed_points = [fix])
	@time ϕ, ω = get_eigenmodes(acft; fixed_points = [fix])
	@time ϕ, ω = get_eigenmodes(acft; fixed_points = [fix])
	@time ϕ, ω = get_eigenmodes(acft; fixed_points = [fix])

	#=
	for i = 1:size(ϕ, 2)
		@show ω[i]

		plot_aircraft(acft, q0 .+ ϕ[:, i])
	end
	=#

	@info "Solve with assumed modes method"

	q .= q0

	A, eig = assumed_modes_solve(acft, q, U∞; fixed_points = [fix])

	#=
	@show eig.λ, eig.ϕ

	for i = 1:size(eig.ϕ, 2)
		dq = A * real.(eig.ϕ[:, i])
		dq ./= norm(dq)

		plot_aircraft(acft, q .+ dq)
	end
	=#

	@time A, eig = assumed_modes_solve(acft, q, U∞; fixed_points = [fix])
	@time A, eig = assumed_modes_solve(acft, q, U∞; fixed_points = [fix])
	@time A, eig = assumed_modes_solve(acft, q, U∞; fixed_points = [fix])
	@time A, eig = assumed_modes_solve(acft, q, U∞; fixed_points = [fix])
	@time A, eig = assumed_modes_solve(acft, q, U∞; fixed_points = [fix])

	true
end
