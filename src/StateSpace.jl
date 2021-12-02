export state_space, ndofs, get_state,
	get_u, get_̇u, get_v, get_̇v, get_w, get_̇w,
	get_θx, get_θy, get_θz, get_̇θx, get_̇θy, get_̇θz,
	set_u, set_v, set_w, set_̇u, set_̇v, set_̇w,
	set_θx, set_θy, set_θz, set_̇θx, set_̇θy, set_̇θz

"""
```
	ndofs(acft::Aircraft)
```

Get number of degrees of freedom for an aircraft
"""
@inline ndofs(acft::Aircraft{Fg}) where Fg = size(acft.points, 2) * 12

"""
```
	get_u(q::AbstractVector, i::Int64)
```

Get x axis position of a point in the flexible structure
"""
@inline get_u(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 1]
"""
```
	get_v(q::AbstractVector, i::Int64)
```

Get y axis position of a point in the flexible structure
"""
@inline get_v(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 2]
"""
```
	get_w(q::AbstractVector, i::Int64)
```

Get z axis position of a point in the flexible structure
"""
@inline get_w(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 3]

"""
```
	get_θx(q::AbstractVector, i::Int64)
```

Get x axis angular defflection of a point in the flexible structure
"""
@inline get_θx(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 4]
"""
```
	get_θy(q::AbstractVector, i::Int64)
```

Get y axis angular defflection of a point in the flexible structure
"""
@inline get_θz(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 5]
"""
```
	get_θz(q::AbstractVector, i::Int64)
```

Get z axis angular defflection of a point in the flexible structure
"""
@inline get_θz(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 6]

"""
```
	get_̇u(q::AbstractVector, i::Int64)
```

Get derivative in time of x axis position of a point in the flexible structure
"""
@inline get_̇u(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 1 + length(q) ÷ 2]
"""
```
	get_̇v(q::AbstractVector, i::Int64)
```

Get derivative in time of y axis position of a point in the flexible structure
"""
@inline get_̇v(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 2 + length(q) ÷ 2]
"""
```
	get_̇w(q::AbstractVector, i::Int64)
```

Get derivative in time of z axis position of a point in the flexible structure
"""
@inline get_̇w(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 3 + length(q) ÷ 2]

"""
```
	get_̇̇θx(q::AbstractVector, i::Int64)
```

Get derivative in time of x axis angular defflection of a point in the flexible structure
"""
@inline get_̇θx(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 4 + length(q) ÷ 2]
"""
```
	get_̇̇θy(q::AbstractVector, i::Int64)
```

Get derivative in time of y axis angular defflection of a point in the flexible structure
"""
@inline get_̇θy(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 5 + length(q) ÷ 2]
"""
```
	get_̇̇θz(q::AbstractVector, i::Int64)
```

Get derivative in time of z axis angular defflection of a point in the flexible structure
"""
@inline get_̇θz(q::AbstractVector, i::Int64) = q[(i - 1) * 6 + 6 + length(q) ÷ 2]

"""
```
	set_u(q::AbstractVector, i::Int64, v)
```

Set x axis defflection of a point in the flexible structure
"""
@inline set_u(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 1] = v end
"""
```
	set_v(q::AbstractVector, i::Int64, v)
```

Set y axis defflection of a point in the flexible structure
"""
@inline set_v(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 2] = v end
"""
```
	set_w(q::AbstractVector, i::Int64, v)
```

Set z axis defflection of a point in the flexible structure
"""
@inline set_w(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 3] = v end

"""
```
	set_θx(q::AbstractVector, i::Int64, v)
```

Set x axis angular defflection of a point in the flexible structure
"""
@inline set_θx(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 4] = v end
"""
```
	set_θy(q::AbstractVector, i::Int64, v)
```

Set y axis angular defflection of a point in the flexible structure
"""
@inline set_θy(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 5] = v end
"""
```
	set_θz(q::AbstractVector, i::Int64, v)
```

Set z axis angular defflection of a point in the flexible structure
"""
@inline set_θz(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 6] = v end

"""
```
	set_̇u(q::AbstractVector, i::Int64, v)
```

Set derivative in time of x axis defflection of a point in the flexible structure
"""
@inline set_̇u(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 1 + length(q) ÷ 2] = v end
"""
```
	set_̇v(q::AbstractVector, i::Int64, v)
```

Set derivative in time of y axis defflection of a point in the flexible structure
"""
@inline set_̇v(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 2 + length(q) ÷ 2] = v end
"""
```
	set_̇w(q::AbstractVector, i::Int64, v)
```

Set derivative in time of z axis defflection of a point in the flexible structure
"""
@inline set_̇w(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 3 + length(q) ÷ 2] = v end

"""
```
	set_̇θx(q::AbstractVector, i::Int64, v)
```

Set derivative in time of x axis angular defflection of a point in the flexible structure
"""
@inline set_̇θx(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 4 + length(q) ÷ 2] = v end
"""
```
	set_̇θy(q::AbstractVector, i::Int64, v)
```

Set derivative in time of y axis angular defflection of a point in the flexible structure
"""
@inline set_̇θy(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 5 + length(q) ÷ 2] = v end
"""
```
	set_̇θz(q::AbstractVector, i::Int64, v)
```

Set derivative in time of z axis angular defflection of a point in the flexible structure
"""
@inline set_̇θz(q::AbstractVector, i::Int64, v) = begin q[(i - 1) * 6 + 6 + length(q) ÷ 2] = v end

"""
```
	get_u(q::AbstractVector)
```

Get vector of x axis defflection of all points in the flexible structure
"""
@inline get_u(q::AbstractVector) = q[1:6:(length(q) ÷ 2)]
"""
```
	get_v(q::AbstractVector)
```

Get vector of y axis defflection of all points in the flexible structure
"""
@inline get_v(q::AbstractVector) = q[2:6:(length(q) ÷ 2)]
"""
```
	get_w(q::AbstractVector)
```

Get vector of z axis defflection of all points in the flexible structure
"""
@inline get_w(q::AbstractVector) = q[3:6:(length(q) ÷ 2)]

"""
```
	get_θx(q::AbstractVector)
```

Get vector of z axis defflection of all points in the flexible structure
"""
@inline get_θx(q::AbstractVector) = q[4:6:(length(q) ÷ 2)]
"""
```
	get_θy(q::AbstractVector)
```

Get vector of y axis defflection of all points in the flexible structure
"""
@inline get_θy(q::AbstractVector) = q[5:6:(length(q) ÷ 2)]
"""
```
	get_θz(q::AbstractVector)
```

Get vector of z axis defflection of all points in the flexible structure
"""
@inline get_θz(q::AbstractVector) = q[6:6:(length(q) ÷ 2)]

"""
```
	get_̇u(q::AbstractVector)
```

Get vector of derivatives in time of x axis defflection of all points in the flexible structure
"""
@inline get_̇u(q::AbstractVector) = q[(1 + length(q) ÷ 2):6:end]
"""
```
	get_̇v(q::AbstractVector)
```

Get vector of derivatives in time of y axis defflection of all points in the flexible structure
"""
@inline get_̇v(q::AbstractVector) = q[(2 + length(q) ÷ 2):6:end]
"""
```
	get_̇w(q::AbstractVector)
```

Get vector of derivatives in time of z axis defflection of all points in the flexible structure
"""
@inline get_̇w(q::AbstractVector) = q[(3 + length(q) ÷ 2):6:end]

"""
```
	get_̇θx(q::AbstractVector)
```

Get vector of derivatives in time of x axis angular defflection of all points in the flexible structure
"""
@inline get_̇θx(q::AbstractVector) = q[(4 + length(q) ÷ 2):6:end]
"""
```
	get_̇θy(q::AbstractVector)
```

Get vector of derivatives in time of y axis angular defflection of all points in the flexible structure
"""
@inline get_̇θy(q::AbstractVector) = q[(5 + length(q) ÷ 2):6:end]
"""
```
	get_̇θz(q::AbstractVector)
```

Get vector of derivatives in time of z axis angular defflection of all points in the flexible structure
"""
@inline get_̇θz(q::AbstractVector) = q[(6 + length(q) ÷ 2):6:end]

"""
```
	set_u(q::AbstractVector, v)
```

Set value of all x axis defflections in the structure
"""
@inline set_u(q::AbstractVector, v) = begin q[1:6:(length(q) ÷ 2)] .= v end
"""
```
	set_v(q::AbstractVector, v)
```

Set value of all y axis defflections in the structure
"""
@inline set_v(q::AbstractVector, v) = begin q[2:6:(length(q) ÷ 2)] .= v end
"""
```
	set_w(q::AbstractVector, v)
```

Set value of all z axis defflections in the structure
"""
@inline set_w(q::AbstractVector, v) = begin q[3:6:(length(q) ÷ 2)] .= v end

"""
```
	set_θx(q::AbstractVector, v)
```

Set value of all x axis angular defflections in the structure
"""
@inline set_θx(q::AbstractVector, v) = begin q[4:6:(length(q) ÷ 2)] .= v end
"""
```
	set_θy(q::AbstractVector, v)
```

Set value of all y axis angular defflections in the structure
"""
@inline set_θy(q::AbstractVector, v) = begin q[5:6:(length(q) ÷ 2)] .= v end
"""
```
	set_θz(q::AbstractVector, v)
```

Set value of all z axis angular defflections in the structure
"""
@inline set_θz(q::AbstractVector, v) = begin q[6:6:(length(q) ÷ 2)] .= v end

"""
```
	set_̇u(q::AbstractVector, v)
```

Set time derivative of all x axis defflections in the structure
"""
@inline set_̇u(q::AbstractVector, v) = begin q[(1 + length(q) ÷ 2):6:end] .= v end
"""
```
	set_̇v(q::AbstractVector, v)
```

Set time derivative of all y axis defflections in the structure
"""
@inline set_̇v(q::AbstractVector, v) = begin q[(2 + length(q) ÷ 2):6:end] .= v end
"""
```
	set_̇w(q::AbstractVector, v)
```

Set time derivative of all z axis defflections in the structure
"""
@inline set_̇w(q::AbstractVector, v) = begin q[(3 + length(q) ÷ 2):6:end] .= v end

"""
```
	set_̇θx(q::AbstractVector, v)
```

Set time derivative of all x axis angular defflections in the structure
"""
@inline set_̇θx(q::AbstractVector, v) = begin q[(4 + length(q) ÷ 2):6:end] .= v end
"""
```
	set_̇θy(q::AbstractVector, v)
```

Set time derivative of all y axis angular defflections in the structure
"""
@inline set_̇θy(q::AbstractVector, v) = begin q[(5 + length(q) ÷ 2):6:end] .= v end
"""
```
	set_̇θz(q::AbstractVector, v)
```

Set time derivative of all z axis angular defflections in the structure
"""
@inline set_̇θz(q::AbstractVector, v) = begin q[(6 + length(q) ÷ 2):6:end] .= v end

"""
Function to get incidences for each wing strip
"""
get_incidences(acft::Aircraft{Fg}, q::AbstractVector) where Fg = [
	(
		strp.Mtouni[1, 1] * (get_θx(q, strp.ipt1) + get_θx(q, strp.ipt2)) / 2 +
		strp.Mtouni[2, 1] * (get_θy(q, strp.ipt1) + get_θy(q, strp.ipt2)) / 2 +
		strp.Mtouni[3, 1] * (get_θz(q, strp.ipt1) + get_θz(q, strp.ipt2)) / 2
	) for strp in acft.wing_strips
]

function get_freestream_velocities(acft::Aircraft{Fg}, q::AbstractVector) where Fg

	#=
	[
		- (
			[
				(get_̇u(q, strp.ipt1) + get_̇u(q, strp.ipt2)) / 2,
				(get_̇v(q, strp.ipt1) + get_̇v(q, strp.ipt2)) / 2,
				(get_̇w(q, strp.ipt1) + get_̇w(q, strp.ipt2)) / 2
			] .+ cross(
				[
					(get_̇θx(q, strp.ipt1) + get_̇θx(q, strp.ipt2)) / 2,
					(get_̇θy(q, strp.ipt1) + get_̇θy(q, strp.ipt2)) / 2,
					(get_̇θz(q, strp.ipt1) + get_̇θz(q, strp.ipt2)) / 2
				],
				[strp.c / 2, 0.0, 0.0]
			)
		) for strp in acft.wing_strips
	]
	=#

	[
		(
			Real[
				- (get_̇u(q, strp.ipt1) + get_̇u(q, strp.ipt2)) / 2,
				- (get_̇v(q, strp.ipt1) + get_̇v(q, strp.ipt2)) / 2,
				- (get_̇w(q, strp.ipt1) + get_̇w(q, strp.ipt2)) / 2
			]
		) for strp in acft.wing_strips
	]

end

"""
Function to get vortex kins, collocation points and normals
"""
function vortex_info(acft::Aircraft{Fg}, q::AbstractVector, θ::AbstractVector) where Fg

	kins = [
		(
			[
				get_u(q, strp.ipt1),
				get_v(q, strp.ipt1),
				get_w(q, strp.ipt1)
			],
			[
				get_u(q, strp.ipt2),
				get_v(q, strp.ipt2),
				get_w(q, strp.ipt2)
			]
		) for strp in acft.wing_strips
	]

	cols = [
		(
			(k[1] .+ k[2]) ./ 2 # .+ [strp.c / 2, 0.0, 0.0]
		) for (k, strp) in zip(kins, acft.wing_strips)
	]

	norms = [
		strp.Mtouni * [
			0.0,
			sin(θt),
			- cos(θt)
		] for (θt, strp) in zip(θ, acft.wing_strips)
	]

	kins, cols, norms

end

"""
Get fuselage section singularities
"""
function fuselage_singularities(acft::Aircraft{Fg}, q::AbstractVector) where Fg

	ls = [
		(
			acft.points[:, sect.ipt2] .- acft.points[:, sect.ipt1]
		) for sect in acft.fuselage_cuts
	]

	Ls = [
		begin
			L = norm(l)

			l ./= L

			L
		end for l in ls
	]

	σ = [
		begin
			V = [
				- (get_̇u(q, sect.ipt1) + get_̇u(q, sect.ipt2)) / 2,
				- (get_̇v(q, sect.ipt1) + get_̇v(q, sect.ipt2)) / 2,
				- (get_̇w(q, sect.ipt1) + get_̇w(q, sect.ipt2)) / 2
			]

			(
				(sect.A2 - sect.A1) * (
					ls[i] ⋅ V
				) + Ls[i] * (sqrt(sect.A2 * π) + sqrt(sect.A1 * π)) * sect.Cf * norm(V) / 2
			) / Ls[i]
		end for (i, sect) in enumerate(acft.fuselage_cuts)
	]

	μ = [
		(sect.A1 + sect.A2) .* (
			begin
				ξ = [
					(get_̇u(q, sect.ipt1) + get_̇u(q, sect.ipt2)) / 2,
					(get_̇v(q, sect.ipt1) + get_̇v(q, sect.ipt2)) / 2,
					(get_̇w(q, sect.ipt1) + get_̇w(q, sect.ipt2)) / 2
				]

				ξ .- (ξ ⋅ ls[i]) .* ls[i]
			end
		) for (i, sect) in enumerate(acft.fuselage_cuts)
	]

	σ, μ

end

"""
Find singularity values for punctual drag components
"""
function point_drag_singularities(acft::Aircraft{Fg}, q::AbstractVector) where Fg

	vs = [
		[
			get_̇u(q, pmass.ipt),
			get_̇v(q, pmass.ipt),
			get_̇w(q, pmass.ipt)
		] for pmass in acft.point_drags
	]

	σ = [
		norm(v) * pmass.A * pmass.CD for (v, pmass) in zip(vs, acft.point_drags)
	]

	μ = [
		v .* pmass.V for (v, pmass) in zip(vs, acft.point_drags)
	]

	σ, μ

end

"""
```
	function get_state(
		acft::Aircraft{Fg},
		U∞::Fg;
		α::Fg = 0.0,
		β::Fg = 0.0,
		p::Fg = 0.0,
		q::Fg = 0.0,
		r::Fg = 0.0
	) where Fg
```

Function to instantiate state vector with a given flight condition

* `acft`: aircraft at hand
* `U∞`: freestream velocity
* `α`: angle of attack
* `β`: sideslip angle
* `p`: x axis angular velocity (in the stability coordinate system)
* `q`: y axis angular velocity (in the stability coordinate system)
* `q`: z axis angular velocity (in the stability coordinate system)
"""
function get_state(
	acft::Aircraft{Fg},
	U∞::Fg;
	α::Fg = 0.0,
	β::Fg = 0.0,
	p::Fg = 0.0,
	q::Fg = 0.0,
	r::Fg = 0.0
) where Fg

	V = [
		- cosd(α) * cosd(β) * U∞,
		- cosd(α) * sind(β) * U∞,
		- sind(α) * U∞
	]
	ω = [
		- p,
		q,
		- r
	]

	vcat(
		[
			[
				acft.points[:, ipt];
				0.0;
				0.0;
				0.0
			] for ipt in 1:size(acft.points, 2)
		]...,
		[
			[
				V .+ ω × acft.points[:, ipt];
				ω
			] for ipt = 1:size(acft.points, 2)
		]...
	)

end

"""
Get coefficients for wing sections
"""
function get_section_forces(
	acft::Aircraft{Fg},
	Γs, U∞,
	us, u0s, ρ, kins,
	∂Γ!∂t
) where Fg

	Γr = (
		isnothing(∂Γ!∂t) ?
		Γs :
		[
			s.c * Γt / norm(u) for (i, (u, s, Γt)) in enumerate(zip(us, acft.wing_strips, ∂Γ!∂t))
		] .+ Γs
	)

	Fs = [
		cross(
			u, k[2] .- k[1]
		) .* ρ .* Γ for (u, k, Γ) in zip(us, kins, Γr)
	]

	CLs = [
		2 * Γs[i] / (U∞ * sect.c) for (i, sect) in enumerate(acft.wing_strips)
	]
	Cms = [
		sect.afl.Cm0 + sect.afl.Cmα * (CL - sect.afl.CL0) / sect.afl.CLα for (sect, CL) in zip(acft.wing_strips, CLs)
	]
	CDs = [
		sect.afl.polar[3] + sect.afl.polar[2] * CL + sect.afl.polar[1] * CL ^ 2 for (sect, CL) in zip(acft.wing_strips, CLs)
	]

	for (i, (k, u, CD, F, sect)) in enumerate(zip(kins, u0s, CDs, Fs, acft.wing_strips))
		nu = norm(u)

		q = (sect.c * sqrt((k[2][2] - k[1][2]) ^ 2 + (k[2][3] - k[1][3]) ^ 2)) * ρ * nu ^ 2 / 2

		@. F += u * CD * q / nu
		CDs[i] = F ⋅ u / (nu * q)
	end

	Ms = [
		begin
			nu = norm(u)

			q = sect.c ^ 2 * ρ * nu ^ 2 / 2

			M = @. Cm * q * (k[2] - k[1])

			M
		end for (k, u, Cm, sect) in zip(kins, us, Cms, acft.wing_strips)
	]

	Fs, Ms, CLs, CDs, Cms

end

"""
Obtain fuselage loads
"""
function fuselage_loads(acft::Aircraft{Fg}, q::AbstractVector, ρ, σs, PGβ) where Fg

	[
		begin
			V = [
				- (get_̇u(q, sect.ipt1) + get_̇u(q, sect.ipt2)) / 2,
				- (get_̇v(q, sect.ipt1) + get_̇v(q, sect.ipt2)) / 2,
				- (get_̇w(q, sect.ipt1) + get_̇w(q, sect.ipt2)) / 2
			]

			(σ * sect.L * PGβ * ρ) .* V
		end for (sect, σ) in zip(acft.fuselage_cuts, σs)
	]

end

"""
Obtain loads at point drag elements
"""
function point_loads(acft::Aircraft{Fg}, q::AbstractVector, ρ, σs, PGβ) where Fg

	[
		begin
			V = [
				- get_̇u(q, pmass.ipt),
				- get_̇v(q, pmass.ipt),
				- get_̇w(q, pmass.ipt)
			]

			(σ * PGβ * ρ) .* V
		end for (pmass, σ) in zip(acft.point_drags, σs)
	]

end

"""
```
	function state_space(
		acft::Aircraft{Fg}, x::AbstractVector, U∞::Fg;
		a::Fg = 343.0,
		ρ::Fg = 1.225,
		α::Fg = 0.0,
		β::Fg = 0.0,
		p::Fg = 0.0,
		q::Fg = 0.0,
		r::Fg = 0.0,
		Sref::Fg = 1.0,
		cref::Fg = 1.0,
		bref::Fg = 1.0,
		∂Γ!∂t::Union{AbstractVector, Nothing} = nothing,
		fixed_points::Vector{Fi} = Int64[],
		engine_loads::Matrix{Fg} = Matrix{Fg}(undef, 2, 0),
		amplification_factor::Float64 = 1.0,
		aerodynamic_forces::Bool = true
	) where {Fg, Fi <: Int}
```

Function to get a dynamic analysis of an aircraft, given a state vector

* `acft`: aircraft object
* `q`: variables of interest
* `U∞`: freestream velocity
* `a`: sound speed
* `ρ`: density
* `α`: angle of attack
* `β`: sideslip angle
* `p`: x axis angular velocity (in the stability coordinate system)
* `q`: y axis angular velocity (in the stability coordinate system)
* `q`: z axis angular velocity (in the stability coordinate system)
* `Sref`: reference surface
* `cref`: reference chord
* `bref`: reference span
* `∂Γ!∂t`: derivative of circulation over time, for unsteady simulations
* `engine_loads`: forces and moments at the engines, as in:

```
engine_loads = [
	F1 F2 ...;
	M1 M2 ...
]
```

* output: derivatives of variables of interest (altered in place)
* output: named tuple as in example:

```
(
	Γ = Γs,
	fs = (
		[ # elastic forces and moments acting at each beam segment
			[fx1, fy1, fz1, mx1, my1, mz1, fx2, fy2, fz2, mx2, my2, mz2],
			[fx1, fy1, fz1, mx1, my1, mz1, fx2, fy2, fz2, mx2, my2, mz2],
			...
		]
	)
	sectional = (
		CL = CLs, 
		CD = CDs, 
		Cm = Cms
	), # sectional coefficients for each strip
	total = (
		CL = CL,
		CD = CD,
		Cl = Cl,
		Cm = Cm,
		Cn = Cn,
		CX = CX,
		CZ = CY,
		CZ = CZ
	)
)
```
"""
function state_space(
	acft::Aircraft{Fg}, x::AbstractVector, U∞::Fg;
	a::Fg = 343.0,
	ρ::Fg = 1.225,
	α::Fg = 0.0,
	β::Fg = 0.0,
	p::Fg = 0.0,
	q::Fg = 0.0,
	r::Fg = 0.0,
	Sref::Fg = 1.0,
	cref::Fg = 1.0,
	bref::Fg = 1.0,
	∂Γ!∂t::Union{AbstractVector, Nothing} = nothing,
	fixed_points::Vector{Fi} = Int64[],
	engine_loads::Matrix{Fg} = Matrix{Fg}(undef, 2, 0),
	amplification_factor::Float64 = 1.0,
	aerodynamic_forces::Bool = true
) where {Fg, Fi <: Int}

	xhat = [
		1.0, 0.0, 0.0
	]
	θ = get_incidences(acft, x)
	kins, cols, norms = vortex_info(acft, x, θ)

	_engine_loads = hcat(
		engine_loads,
		zeros(Fg, 2, length(acft.engines) - size(engine_loads, 2))
	)

	us = get_freestream_velocities(acft, x)

	add_engine_influence!(
		acft,
		x,
		_engine_loads,
		cols,
		us,
		ρ
	)

	σ, μ = fuselage_singularities(acft, x)
	σ_point, μ_point = point_drag_singularities(acft, x)

	if aerodynamic_forces
		for (ipt, pmass) in enumerate(acft.point_drags)
			pt = [
				get_u(x, pmass.ipt),
				get_v(x, pmass.ipt),
				get_w(x, pmass.ipt)
			]

			for (u, colpt) in zip(us, cols)
				u .+= σ_point[ipt] .* point_source_infl(
					colpt, pt
				)
				u .+= begin
					nμ = norm(μ_point[ipt])
					ξ = μ_point[ipt] ./ (nμ + 1e-10)

					nμ .* point_doublet_infl(colpt, pt, ξ)
				end
			end
		end

		for (ifus, sect) in enumerate(acft.fuselage_cuts)
			pt1 = [
				get_u(x, sect.ipt1),
				get_v(x, sect.ipt1),
				get_w(x, sect.ipt1)
			]
			pt2 = [
				get_u(x, sect.ipt2),
				get_v(x, sect.ipt2),
				get_w(x, sect.ipt2)
			]

			for (u, colpt) in zip(us, cols)
				u .+= σ[ifus] .* source_line(
					colpt, pt1, pt2
				)
				u .+= begin
					nμ = norm(μ[ifus])
					η = μ[ifus] ./ (nμ + 1e-10)

					nμ .* doublet_line(
						colpt, pt1, pt2, η
					)
				end
			end
		end
	end

	M∞ = U∞ / a
	PGβ = sqrt(1.0 - M∞ ^ 2)

	Γ, u0s = (
		aerodynamic_forces ? 
		begin
			A, b, As = linear_problem(acft, kins, cols, norms, us, xhat, PGβ)

			Γ = (A \ b) ./ PGβ
			for i = 1:length(acft.fuselage_cuts)
				σ[i] /= PGβ
				μ[i] ./= PGβ
			end
			for i = 1:length(acft.point_drags)
				σ_point[i] /= PGβ
				μ_point[i] ./= PGβ
			end

			u0s = deepcopy(us)
			add_downwash!(us, Γ, As)

			Γ, u0s
		end :
		(zeros(eltype(x), length(acft.wing_strips)), us)
	)

	Fs, Ms, CLs, CDs, Cms = get_section_forces(acft, Γ, U∞, us, u0s, ρ, kins, ∂Γ!∂t)

	fFs = fuselage_loads(acft, x, ρ, σ, PGβ)
	pFs = point_loads(acft, x, ρ, σ_point, PGβ)

	f = zeros(Real, length(x) ÷ 2)

	if aerodynamic_forces
		for (strip, F, M) in zip(acft.wing_strips, Fs, Ms)
			f1inds = (6 * (strip.ipt1 - 1) + 1):(6 * (strip.ipt1 - 1) + 3)
			f2inds = (6 * (strip.ipt2 - 1) + 1):(6 * (strip.ipt2 - 1) + 3)

			m1inds = (6 * (strip.ipt1 - 1) + 4):(6 * strip.ipt1)
			m2inds = (6 * (strip.ipt2 - 1) + 4):(6 * strip.ipt2)

			Ft = F ./ 2
			Mt = cross(
				Ft,
				[strip.c * strip.xCG, 0.0, 0.0]
			) .+ M ./ 2

			m_sect = (acft.masses[strip.ipt1] + acft.masses[strip.ipt2]) / 2
			inv_inertia = inv(m_sect.I) .* m_sect.m
			df = [strip.xCG * strip.c, 0.0, 0.0] × (inv_inertia * Mt)

			f[f1inds] .+= Ft .+ df
			f[f2inds] .+= Ft .+ df

			f[m1inds] .+= Mt
			f[m2inds] .+= Mt
		end

		for (cut, F) in zip(acft.fuselage_cuts, fFs)
			f1inds = (6 * (cut.ipt1 - 1) + 1):(6 * (cut.ipt1 - 1) + 3)
			f2inds = (6 * (cut.ipt2 - 1) + 1):(6 * (cut.ipt2 - 1) + 3)

			f[f1inds] .+= F ./ 2
			f[f2inds] .+= F ./ 2
		end

		for (pmass, F) in zip(acft.point_drags, pFs)
			f[(6 * (pmass.ipt - 1) + 1):(6 * (pmass.ipt - 1) + 3)] .+= F
		end
	end

	is_calculated = zeros(Bool, length(acft.beams))
	elastic_effort = Vector{Any}(undef, length(acft.beams))

	for sect in acft.wing_strips
		is_calculated[sect.ibeam] = true

		θ1 = [
			get_θx(x, sect.ipt1),
			get_θy(x, sect.ipt1),
			get_θz(x, sect.ipt1)
		]
		θ2 = [
			get_θx(x, sect.ipt2),
			get_θy(x, sect.ipt2),
			get_θz(x, sect.ipt2)
		]

		x1 = [
			get_u(x, sect.ipt1),
			get_v(x, sect.ipt1),
			get_w(x, sect.ipt1)
		] .- acft.points[:, sect.ipt1]
		x2 = [
			get_u(x, sect.ipt2),
			get_v(x, sect.ipt2),
			get_w(x, sect.ipt2)
		] .- acft.points[:, sect.ipt2]

		x1 .+= cross(
			θ1,
			[sect.c * sect.e, 0.0, 0.0]
		)
		x2 .+= cross(
			θ2,
			[sect.c * sect.e, 0.0, 0.0]
		)

		b = acft.beams[sect.ibeam]

		pl = b.K * [
			b.Mtosys * x1;
			b.Mtosys * θ1;
			b.Mtosys * x2;
			b.Mtosys * θ2
		]

		elastic_effort[sect.ibeam] = pl

		f1inds = (6 * (sect.ipt1 - 1) + 1):(6 * (sect.ipt1 - 1) + 3)
		f2inds = (6 * (sect.ipt2 - 1) + 1):(6 * (sect.ipt2 - 1) + 3)

		m1inds = (6 * (sect.ipt1 - 1) + 4):(6 * sect.ipt1)
		m2inds = (6 * (sect.ipt2 - 1) + 4):(6 * sect.ipt2)

		m_sect = (acft.masses[sect.ipt1] + acft.masses[sect.ipt2]) / 2
		inv_inertia = inv(m_sect.I) .* m_sect.m

		fb = b.Mtouni * pl[1:3]
		mb = cross(
			fb, [sect.c * (sect.xCG - sect.e), 0.0, 0.0]
		) .+ b.Mtouni * pl[4:6]
		f[f1inds] .+= fb .+ [sect.c * sect.xCG, 0.0, 0.0] × (inv_inertia * mb)
		f[m1inds] .+= mb

		fb .= b.Mtouni * pl[7:9]
		mb .= cross(
			fb, [sect.c * (sect.xCG - sect.e), 0.0, 0.0]
		) .+ b.Mtouni * pl[10:12]
		f[f2inds] .+= fb .+ [sect.c * sect.xCG, 0.0, 0.0] × (inv_inertia * mb)
		f[m2inds] .+= mb
	end

	for (ibeam, beam) in enumerate(acft.beams)
		if !is_calculated[ibeam]
			is_calculated[ibeam] = true

			θ1 = [
				get_θx(x, beam.ipt1),
				get_θy(x, beam.ipt1),
				get_θz(x, beam.ipt1)
			]
			θ2 = [
				get_θx(x, beam.ipt2),
				get_θy(x, beam.ipt2),
				get_θz(x, beam.ipt2)
			]

			x1 = [
				get_u(x, beam.ipt1),
				get_v(x, beam.ipt1),
				get_w(x, beam.ipt1)
			] .- acft.points[:, beam.ipt1]
			x2 = [
				get_u(x, beam.ipt2),
				get_v(x, beam.ipt2),
				get_w(x, beam.ipt2)
			] .- acft.points[:, beam.ipt2]

			pl = beam.K * [
				beam.Mtosys * x1;
				beam.Mtosys * θ1;
				beam.Mtosys * x2;
				beam.Mtosys * θ2
			]

			elastic_effort[ibeam] = pl

			f1inds = (6 * (beam.ipt1 - 1) + 1):(6 * (beam.ipt1 - 1) + 3)
			f2inds = (6 * (beam.ipt2 - 1) + 1):(6 * (beam.ipt2 - 1) + 3)

			m1inds = (6 * (beam.ipt1 - 1) + 4):(6 * beam.ipt1)
			m2inds = (6 * (beam.ipt2 - 1) + 4):(6 * beam.ipt2)

			f[f1inds] .+= beam.Mtouni * pl[1:3]
			f[f2inds] .+= beam.Mtouni * pl[7:9]

			f[m1inds] .+= beam.Mtouni * pl[4:6]
			f[m2inds] .+= beam.Mtouni * pl[10:12]
		end
	end

	for (ie, eng) in enumerate(acft.engines)
		inds = (6 * (eng.ipt - 1) + 1):(6 * eng.ipt)

		θ = [
			get_θx(x, eng.ipt),
			get_θy(x, eng.ipt),
			get_θz(x, eng.ipt)
		]
		ξ⃗ = eng.ξ⃗ .+ cross(θ, eng.ξ⃗)
		ξ⃗ ./= norm(ξ⃗)

		f[inds] .+= [
			ξ⃗ .* _engine_loads[1, ie];
			ξ⃗ .* _engine_loads[2, ie]
		]
	end

	qd = similar(f)

	V = [
		- cosd(α) * cosd(β) * U∞,
		- cosd(α) * sind(β) * U∞,
		- sind(α) * U∞
	]
	ω = [
		- p,
		q,
		- r
	]

	for (im, mass) in enumerate(acft.masses)
		finds = (6 * (im - 1) + 1):(6 * (im - 1) + 3)
		minds = (6 * (im - 1) + 4):(6 * im)

		rp = [
			get_u(x, im),
			get_v(x, im),
			get_w(x, im)
		]
		vp = [
			get_̇u(x, im),
			get_̇v(x, im),
			get_̇w(x, im)
		]

		qd[finds] .= (
			f[finds] .- 
			mass.m .* cross(ω, cross(ω, rp)) .- 
			(2 * mass.m) .* cross(ω, vp)
		) ./ mass.m
		qd[minds] .= inv(mass.I) * f[minds]
	end

	qdot = [
		begin
			vs = x[(length(x) ÷ 2 + 1):end]

			for i = 1:6:length(vs)
				vs[i:(i + 2)] .-= V
				vs[(i + 3):(i + 5)] .-= ω
			end

			vs
		end;
		qd
	]

	amp = amplification_factor
	for f in fixed_points
		set_u(qdot, f, (acft.points[1, f] - get_u(x, f)) * amp)
		set_v(qdot, f, (acft.points[2, f] - get_v(x, f)) * amp)
		set_w(qdot, f, (acft.points[3, f] - get_w(x, f)) * amp)
		set_θx(qdot, f, - get_θx(x, f) * amp)
		set_θy(qdot, f, - get_θy(x, f) * amp)
		set_θz(qdot, f, - get_θz(x, f) * amp)

		set_̇u(qdot, f, (V[1] - get_̇u(x, f)) * amp)
		set_̇v(qdot, f, (V[2] - get_̇v(x, f)) * amp)
		set_̇w(qdot, f, (V[3] - get_̇w(x, f)) * amp)
		set_̇θx(qdot, f, (ω[1] - get_̇θx(x, f)) * amp)
		set_̇θy(qdot, f, (ω[2] - get_̇θy(x, f)) * amp)
		set_̇θz(qdot, f, (ω[3] - get_̇θz(x, f)) * amp)
	end

	F = zeros(Fg, 3)
	M = zeros(Fg, 3)
	if length(acft.wing_strips) > 0
		F = copy(Fs[1])
		M = Ms[1] .+ cross(
			(kins[1][1] .+ kins[1][2]) ./ 2, Fs[1]
		)
		for i = 2:length(acft.wing_strips)
			F .+= Fs[i]
			M .+= Ms[i] .+ cross(
				(kins[i][1] .+ kins[i][2]) ./ 2, Fs[i]
			)
		end
	end

	for (i, pF) in enumerate(pFs)
		F .+= pF
		M .+= cross(
			begin
				ipt = acft.point_drags[i].ipt

				[
					get_u(x, ipt),
					get_v(x, ipt),
					get_w(x, ipt)
				]
			end,
			pF
		)
	end

	for (i, fF) in enumerate(fFs)
		F .+= fF
		M .+= cross(
			begin
				ipt1 = acft.fuselage_cuts[i].ipt1
				ipt2 = acft.fuselage_cuts[i].ipt2

				[
					(get_u(x, ipt1) + get_u(x, ipt2)) / 2,
					(get_v(x, ipt1) + get_v(x, ipt2)) / 2,
					(get_w(x, ipt1) + get_w(x, ipt2)) / 2
				]
			end,
			fF
		)
	end

	q = U∞ ^ 2 * ρ / 2

	CX = - F[1] / (Sref * q)
	CY = F[2] / (Sref * q)
	CZ = - F[3] / (Sref * q)

	xv = [
		cosd(α) * cosd(β),
		cosd(α) * sind(β),
		sind(α)
	]

	CF = [
		- CX,
		CY,
		- CZ
	]

	CD = CF ⋅ xv
	CL = norm(
		CF .- CD .* xv
	)

	Cl = - M[1] / (q * Sref * bref)
	Cm = M[2] / (q * Sref * cref)
	Cn = - M[3] / (q * Sref * bref)

	is_stalled = [
		cl < s.afl.CLmin || cl > s.afl.CLmax for (cl, s) in zip(CLs, acft.wing_strips)
	]

	qdot, (
		Γ = Γ,
		sectional = (
			CL = CLs,
			CD = CDs,
			Cm = Cms,
			is_stalled = is_stalled
		),
		total = (
			CX = CX,
			CY = CY,
			CZ = CZ,
			CL = CL,
			CD = CD,
			Cl = Cl,
			Cm = Cm,
			Cn = Cn
		),
		f = elastic_effort
	)

end
