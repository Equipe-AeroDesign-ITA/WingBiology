"""
Test code for aerodynamic coefficients of a straight wing
"""

using ForwardDiff

@test begin
	NW = 65
	NM = 10
	NF = 21

	missile_xs = collect(LinRange(-0.5, 0.0, NM))
	fus_xs = collect(LinRange(- 3.5, 3.5, NF))

	fus_Rs = sqrt.(1.0 .- (fus_xs ./ 3.5) .^ 2) .* 0.2
	fus_As = fus_Rs .^ 2 .* π

	fus_z = 0.1 

	pts = hcat(
		[
			[0.0, y, 0.0] for y in LinRange(- 3.5, 3.5, NW)
		]...
	)

	pts = hcat(
		pts,
		[
			[missile_xs[id], - 3.5, 0.0] for id = 1:(NM - 1)
		]...,
		[
			[missile_xs[id], 3.5, 0.0] for id = 1:(NM - 1)
		]...,
		[
			[fus_xs[id], 0.0, fus_z] for id = 1:NF
		]...
	)

	acft = Aircraft(pts)

	for i = 1:(NW - 1)
		add_wing_strip!(
			acft,
			i,
			i + 1,
			1.0;
			e = 0.25,
			m = 0.2
		)
	end

	AM = 0.01 ^ 2 * π

	for i = 1:(NM - 1)
		add_fuselage_cut!(
			acft,
			i + NW,
			(i == (NM - 1) ? 1 : i + NW + 1),
			(i == 1 ? 0.0 : AM),
			(i == NM - 1 ? 0.0 : AM)
		)

		add_fuselage_cut!(
			acft,
			i + NW + (NM - 1),
			(i == (NM - 1) ? NW : i + NW + NM),
			(i == 1 ? 0.0 : AM),
			(i == NM - 1 ? 0.0 : AM)
		)
	end

	for i = 1:(NF - 1)
		add_fuselage_cut!(
			acft,
			NW + 2 * (NM - 1) + i,
			NW + 2 * (NM - 1) + i + 1,
			fus_As[i],
			fus_As[i + 1]
		)
	end

	add_link!(
		acft,
		NW ÷ 2 + 1,
		NW + 2 * (NM - 1) + NF ÷ 2 + 1
	)

	add_point_drag!(
		acft,
		NW ÷ 2 + 1,
		4 * π * 0.025 ^ 3 / 3
	)

	add_engine!(
		acft,
		NW + 1,
		1.0
	)
	add_engine!(
		acft,
		NW + NM,
		1.0
	)

	fcon = get_state(
		acft,
		100.0;
		α = 5.0
	)

	FM = [
		100.0 100.0
		(- 10.0) 10.0
	]

	@btime fcon_derivative, coefs = state_space($acft, $fcon, 100.0; α = 5.0, Sref = 7.0, bref = 7.0, engine_loads = $FM)
	fcon_derivative, coefs = state_space(acft, fcon, 100.0; α = 5.0, Sref = 7.0, bref = 7.0, engine_loads = FM)

	# plot_aircraft(acft)

	true
end
