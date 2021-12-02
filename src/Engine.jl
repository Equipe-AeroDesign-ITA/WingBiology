const _Kϵ = 0.11

function engine_influence(
	ξ⃗::AbstractVector,
	R::Fg,
	pt::AbstractVector,
	colpt::AbstractVector,
	F::Fg,
	M::Fg,
	v⃗::AbstractVector,
	ρ::Fg
) where Fg

	x = - ξ⃗ ⋅ (colpt .- pt)

	if x < 0.0
		return zeros(Fg, 3)
	end

	r⃗ = begin
		d = colpt .- pt

		d .+ x .* ξ⃗
	end
	r = norm(r⃗)

	A = R ^ 2 * π

	V = - v⃗ ⋅ ξ⃗
	ΔV = sqrt(V ^ 2 + 2 * F / (ρ * A)) - V

	ϵc = _Kϵ * abs(ΔV) / (abs(V + ΔV) + 1e-5)
	ϵb = _Kϵ * abs(ΔV) / (abs(V + 0.5 * ΔV) + 1e-5)

	Rp = R * sqrt((V + ΔV / 2) / (abs(V + ΔV) + 1e-5))

	b = Rp + ϵb * x
	c = max(
		0.0,
		Rp - ϵc * x
	)

	k₁ = c ^ 2 + 0.9 * c * (b - c) + 9 * (b - c) ^ 2 / 35
	k₂ = c ^ 2 + (243 / 385) * c * (b - c) + (243 / 1820) * (b - c) ^ 2

	ΔV0 = sqrt(0.25 * (k₁ / k₂) ^ 2 * V ^ 2 + F / (ρ * π * k₂)) - 0.5 * (k₁ / k₂) * V
	ΔV = (
		r > b ?
		0.0 :
		(
			r < c ?
			ΔV0 :
			ΔV0 * (1.0 - ((r - c) / (b - c)) ^ 1.5) ^ 2
		)
	)

	k₃ = 2 * c ^ 3 / 3 + 0.9 * c ^ 2 * (b - c) + 9 * c * (b - c) ^ 2 / 35 + (b - c) ^ 3 / 9
	k₄ = 2 * c ^ 3 / 3 + 243 * c ^ 2 * (b - c) / 385 + 243 * c * (b - c) ^ 2 / 1820 + 2 * (b - c) ^ 3 / 45

	ΔW0 = - M / (π * ρ * (abs(k₃ * V + k₄ * ΔV0) + 1e-5))
	ΔW = (
		r > b ?
		0.0 :
		(
			r < c ?
			ΔW0 :
			ΔW0 * (1.0 - ((r - c) / (b - c)) ^ 1.5) ^ 2
		)
	)

	cross_dir = cross(ξ⃗, r⃗ ./ r)
	@. - ξ⃗ * ΔV + cross_dir * ΔW

end

function add_engine_influence!(
	acft,
	q::AbstractVector,
	FM::AbstractMatrix,
	cols,
	us,
	ρ::Real
)

	Fg = eltype(q)

	Δus = [
		zeros(Fg, 3) for _ in us
	]

	for (ie, eng) in enumerate(acft.engines)
		θ = [
			get_θx(q, eng.ipt),
			get_θy(q, eng.ipt),
			get_θz(q, eng.ipt)
		]

		pt = [
			get_u(q, eng.ipt),
			get_v(q, eng.ipt),
			get_w(q, eng.ipt)
		]

		ξ⃗ = cross(θ, eng.ξ⃗) .+ eng.ξ⃗
		ξ⃗ ./= norm(ξ⃗)

		F, M = FM[:, ie]

		for (colpt, Δu, u) in zip(cols, Δus, us)
			Δu .+= engine_influence(
				ξ⃗,
				eng.R,
				pt,
				colpt,
				F,
				M,
				u,
				ρ
			)
		end
	end

	for (u, Δu) in zip(us, Δus)
		u .+= Δu
	end

end
