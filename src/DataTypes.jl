export WingStrip, FuselageCut, LumpedMass, Beam, Aircraft, Airfoil, PointDrag, 
	add_beam!, add_mass!, add_wing_strip!, add_point_drag!, add_fuselage_cut!, flat_plate,
	Engine, add_engine!, add_point!, interpolate_wing!, interpolate_fuselage!,
	plain_flap, slotted_flap, slat, fixed_slat, Kruger_slat,
	get_massic_properties,
	WingSectInfo, FuselageCutInfo, interpolate!

"""
```
	struct LumpedMass{Fg <: Real} 
		m::Fg
		I::Matrix{Fg}
	end
```

Lumped mass struct
"""
struct LumpedMass{Fg <: Real} 
	m::Fg
	I::Matrix{Fg}
end

import Base: *
"""
Matrix multiplications for lumped mass
"""
function *(LM::LumpedMass{Fg}, M::Matrix{Fc}) where {Fg, Fc <: Real} 
	
	LM.I .= LM.I * M

	LM

end

"""
Matrix multiplications for lumped mass
"""
function *(M::Matrix{Fc}, LM::LumpedMass{Fg}) where {Fg, Fc <: Real} 
	
	LM.I .= M * LM.I

	LM

end

"""
Multiplication between lumped mass and scalar
"""
*(Sc::Fc, M::LumpedMass{Fg}) where {Fg, Fc <: Real} = LumpedMass{Fg}(
	m * Sc,
	M.I .* Sc
)

"""
Multiplication between lumped mass and scalar
"""
*(M::LumpedMass{Fg}, Sc::Fc) where {Fg, Fc <: Real} = Sc * M

import Base: /
"""
Division between lumped mass and scalar
"""
/(M::LumpedMass{Fg}, Sc::Fc) where {Fg, Fc <: Real} = LumpedMass{Fg}(
	M.m / Sc,
	M.I ./ Sc
)

import Base: +
"""
Addition between two lumped masses
"""
+(M1::LumpedMass{Fg}, M2::LumpedMass{Fg}) where Fg = LumpedMass{Fg}(
	M1.m + M2.m,
	M1.I .+ M2.I
)

"""
Abstract type for a fuselage cut or wing strip
"""
abstract type AbstractAeroModel end

"""
```
	struct Airfoil{Fg <: Real}
		CL0::Fg
		CLα::Fg
		Cm0::Fg
		Cmα::Fg
		polar::Vector{Fg}
		CLmin::Fg
		CLmax::Fg
	end
```

Data type defining an airfoil

* `polar` should be an array `[a, b, c]` such that `CD = CL ^ 2 * a + CL * b + c`
"""
struct Airfoil{Fg <: Real}
	CL0::Fg
	CLα::Fg
	Cm0::Fg
	Cmα::Fg
	polar::Vector{Fg}
	CLmin::Fg
	CLmax::Fg
end

Base.:+(afl1::Airfoil, afl2::Airfoil) = Airfoil(
    (afl1.CL0 + afl2.CL0),
    (afl1.CLα + afl2.CLα),
    (afl1.Cm0 + afl2.Cm0),
    (afl1.Cmα + afl2.Cmα),
    (afl1.polar .+ afl2.polar),
	(afl1.CLmin + afl2.CLmin),
	(afl1.CLmax + afl2.CLmax)
)
Base.:*(η::Number, afl::Airfoil) = Airfoil(
    η * afl.CL0,
    η * afl.CLα,
    η * afl.Cm0,
    η * afl.Cmα,
    η .* afl.polar,
	η * afl.CLmin,
	η * afl.CLmax
)
Base.:*(afl::Airfoil, η::Number) = η * afl

Airfoil(
	;
	CL0::Fg = 0.0,
	CLα::Fg = 2 * π,
	Cm0::Fg = 0.0,
	Cmα::Fg = 0.0,
	polar::Vector{Fg} = zeros(3),
	CLmin::Fg = -0.7,
	CLmax::Fg = 0.7
) where Fg = Airfoil{Fg}(
	CL0, 
	CLα, 
	Cm0, 
	Cmα, 
	polar, 
	CLmin, 
	CLmax
)

"""
```
	struct WingStrip{Fg <: Real} <: AbstractAeroModel
		ipt1::Int64
		ipt2::Int64
		ibeam::Int64
		Mtosys::Matrix{Fg}
		Mtouni::Matrix{Fg}
		L::Fg
		c::Fg
		xCG::Fg
		e::Fg
		afl::Airfoil{Fg}
	end
```

Wing strip data type
"""
struct WingStrip{Fg <: Real} <: AbstractAeroModel
	ipt1::Int64
	ipt2::Int64
	ibeam::Int64
	Mtosys::Matrix{Fg}
	Mtouni::Matrix{Fg}
	L::Fg
	c::Fg
	xCG::Fg
	e::Fg
	afl::Airfoil{Fg}
end

"""
```
	struct FuselageCut{Fg <: Real} <: AbstractAeroModel
		ipt1::Int64
		ipt2::Int64
		ibeam::Int64
		Mtosys::Matrix{Fg}
		Mtouni::Matrix{Fg}
		A1::Fg
		A2::Fg
		L::Fg
		Cf::Fg
	end
```

Fuselage cut data type
"""
struct FuselageCut{Fg <: Real} <: AbstractAeroModel
	ipt1::Int64
	ipt2::Int64
	ibeam::Int64
	Mtosys::Matrix{Fg}
	Mtouni::Matrix{Fg}
	A1::Fg
	A2::Fg
	L::Fg
	Cf::Fg
end

"""
```
	struct Beam{Fg <: Real}
		ipt1::Int64
		ipt2::Int64
		K::Matrix{Fg}
		Mtosys::Matrix{Fg}
		Mtouni::Matrix{Fg}
	end
```

Struct to designate a simple beam element
"""
struct Beam{Fg <: Real}
	ipt1::Int64
	ipt2::Int64
	K::Matrix{Fg}
	Mtosys::Matrix{Fg}
	Mtouni::Matrix{Fg}
end

"""
```
	struct PointDrag{Fg <: Real} 
		ipt::Int64
		CD::Fg
		A::Fg
		V::Fg
	end
```

Struct to designate a punctual drag singularity
"""
struct PointDrag{Fg <: Real} 
	ipt::Int64
	CD::Fg
	A::Fg
	V::Fg
end

"""
```
	struct Engine{Fg <: Real}
		ipt::Int64
		ξ⃗::AbstractVector
		R::Fg
	end
```
 
Struct to designate an engine
"""
struct Engine{Fg <: Real}
	ipt::Int64
	ξ⃗::AbstractVector
	R::Fg
end

"""
```
	mutable struct Aircraft{Fg <: Real}
		points::Matrix{Fg}
		wing_strips::Vector{WingStrip{Fg}}
		fuselage_cuts::Vector{FuselageCut{Fg}}
		point_drags::Vector{PointDrag{Fg}}
		engines::Vector{Engine{Fg}}
		beams::Vector{Beam{Fg}}
		masses::Vector{LumpedMass{Fg}}
	end
```

Mutable struct to designate an aircraft
"""
mutable struct Aircraft{Fg <: Real}
	points::Matrix{Fg}
	wing_strips::Vector{WingStrip{Fg}}
	fuselage_cuts::Vector{FuselageCut{Fg}}
	point_drags::Vector{PointDrag{Fg}}
	engines::Vector{Engine{Fg}}
	beams::Vector{Beam{Fg}}
	masses::Vector{LumpedMass{Fg}}
end

"""
```
	Aircraft{Fg}(points::Matrix{Fg})
```

Constructor for an aircraft, initially devoid of wings, fuselage, beams or mass points.
Takes a matrix of points (shape `(3, npts)`)
"""
Aircraft{Fg}(points::Matrix{Fg}) where Fg = Aircraft{Fg}(
	points,
	WingStrip{Fg}[],
	FuselageCut{Fg}[],
	PointDrag{Fg}[],
	Engine{Fg}[],
	Beam{Fg}[],
	[
		LumpedMass{Fg}(
			0.0,
			zeros(Fg, 3, 3)
		) for i = 1:size(points, 2)
	]
)
Aircraft(points::Matrix{Fg}) where {Fg <: Real} = Aircraft{Fg}(points)
"""
```
	Aircraft()
```

Constructor for an aircraft initially devoid of all elements
"""
Aircraft() = Aircraft(Matrix{Float64}(undef, 3, 0))

"""
```
	function add_mass!(
		acft::Aircraft{Fg},
		ipt::Int64,
		m::Fg;
		Ixx::Fg = 0.0,
		Iyy::Fg = 0.0,
		Izz::Fg = 0.0,
		Ixy::Fg = 0.0,
		Ixz::Fg = 0.0,
		Iyz::Fg = 0.0
	) where {Fg}
```

Add point mass to an aircraft.

* `acft`: aircraft at hand
* `ipt`: index for point in which the mass will be added
* `m`: mass
* `Ixx, Iyy, Izz, Ixy, Ixz, Iyz`: moments of inertia.
	Default to zero if absent
"""
function add_mass!(
	acft::Aircraft{Fg},
	ipt::Int64,
	m::Fg;
	Ixx::Fg = 0.0,
	Iyy::Fg = 0.0,
	Izz::Fg = 0.0,
	Ixy::Fg = 0.0,
	Ixz::Fg = 0.0,
	Iyz::Fg = 0.0
) where {Fg}

	acft.masses[ipt] += LumpedMass{Fg}(
		m,
		[
			Ixx Ixy Ixz;
			Ixy Iyy Iyz;
			Ixz Iyz Izz
		]
	)

end

"""
```
	function add_beam!(
		acft::Aircraft{Fg},
		ipt1::Int64,
		ipt2::Int64;
		m::Fg = 1e-3,
		ŷ::Vector{Fg} = [1.0, 0.0, 0.0],
		EA::Fg = Fmax,
		EIy::Fg = Fmax,
		EIz::Fg = Fmax,
		GJ::Fg = Fmax,
		Ixx::Fg = 1e-3
	) where Fg
```

Add a beam to the aircraft model

* `acft`: aircraft struct
* `ipt1`: index for first point
* `ipt2`: index for second point
* `m`: linear mass density
* `EA`: normal stiffness
* `EIy`: y axis flexural stiffness
* `EIz`: z axis flexural stiffness
* `GJ`: x axis torsional stiffness
* `Ixx`: moment of inertia around beam axis (per unit length)
* `ŷ`: y direction used to define the sectional coordinate system
"""
function add_beam!(
	acft::Aircraft{Fg},
	ipt1::Int64,
	ipt2::Int64;
	m::Fg = 1e-3,
	ŷ::Vector{Fg} = [1.0, 0.0, 0.0],
	EA::Fg = Fmax,
	EIy::Fg = Fmax,
	EIz::Fg = Fmax,
	GJ::Fg = Fmax,
	Ixx::Fg = 1e-3
) where Fg

	x̂ = acft.points[:, ipt2] .- acft.points[:, ipt1]

	L = norm(x̂)
	x̂ ./= L

	ŷ .-= x̂ * (x̂ ⋅ ŷ)
	ny = norm(ŷ)
	if ny < 1e-5
		ŷ = [0.0, 1.0, 0.0]
		
		ŷ .-= x̂ * (x̂ ⋅ ŷ)
		ny = norm(ŷ)
	end
	ŷ ./= ny

	ẑ = cross(x̂, ŷ)

	Mtouni = [
		x̂ ŷ ẑ
	]
	Mtosys = inv(Mtouni)

	K = beam_get_K(
		L,
		EA,
		GJ,
		EIy,
		EIz
	)

	push!(acft.beams, Beam{Fg}(ipt1, ipt2, K, Mtosys, Mtouni))

	Ilat = m * L ^ 3 / 3

	LM = LumpedMass{Fg}(
		m * L,
		Mtouni * [
			(Ixx * L) 0.0 0.0;
			0.0 Ilat 0.0;
			0.0 0.0 Ilat
		] * Mtosys
	)

	acft.masses[ipt1] += LM / 2
	acft.masses[ipt2] += LM / 2

	acft.beams[end]

end

"""
```
	function add_point_drag!(
		acft::Aircraft{Fg},
		ipt::Int64,
		V::Fg,
		A::Union{Fg, Nothing} = nothing;
		CD::Fg = 2.0
	) where Fg
```

Function to add a punctual drag source to an aerodynamic model

* `acft`: aircraft model
* `ipt`: point number in aircraft mesh
* `V`: volume
* `A`: surface area (optional, assumed spherical if absent)
* `CD`: drag coeficient in relationship to `A`. Defaults to 2
	(detached flow all throughout the surface)
"""
function add_point_drag!(
	acft::Aircraft{Fg},
	ipt::Int64,
	V::Fg,
	A::Union{Fg, Nothing} = nothing;
	CD::Fg = 2.0
) where Fg

	push!(
		acft.point_drags,
		PointDrag{Fg}(
			ipt,
			CD,
			(isnothing(A) ? 4 * π * (3 * V / (π * 4)) ^ (2 / 3) : A),
			V
		)
	)

	acft.point_drags[end]

end

"""
```
	function add_engine!(
		acft::Aircraft{Fg},
		ipt::Int64,
		R::Fg;
		ξ⃗::AbstractVector = [- 1.0, 0.0, 0.0]
	) where Fg
```

Add engine to aircraft model.

* `acft`: aircraft model
* `ipt`: index for point in aircraft model
* `R`: actuator disk radius
* `ξ⃗`: vector pointing in the direction of the traction and moment vectors.
	Defaults to the opposite of the x vector
"""
function add_engine!(
	acft::Aircraft{Fg},
	ipt::Int64,
	R::Fg;
	ξ⃗::AbstractVector = [- 1.0, 0.0, 0.0]
) where Fg

	push!(
		acft.engines,
		Engine{Fg}(
			ipt,
			ξ⃗ ./ norm(ξ⃗),
			R
		)
	)

	length(acft.engines)

end

"""
Default, flat plate airfoil
"""
flat_plate = Airfoil{Float64}(0.0, 2 * π, 0.0, 0.0, [0.0, 0.0, 0.0], -0.7, 0.7)

"""
```
	function add_wing_strip!(
		acft::Aircraft{Fg},
		ipt1::Int64,
		ipt2::Int64,
		c::Fg;
		incidence::Fg = 0.0,
		afl::Airfoil{Fg} = flat_plate,
		xCG::Fg = 0.25,
		e::Fg = 0.0,
		m::Fg = 1e-3,
		ŷ::Vector{Fg} = [1.0, 0.0, 0.0],
		EA::Fg = Fmax,
		EIy::Fg = Fmax,
		EIz::Fg = Fmax,
		GJ::Fg = Fmax,
		Ixx::Fg = 1e-3
	) where Fg
```

Add wing strip to aircraft model

* `acft`: aircraft struct
* `ipt1`: index for point number 1
* `ipt2`: index for point number 2
* `c`: chord
* `incidence`: incidence (degrees)
* `afl`: airfoil
* `xCG`: position of the center of gravity (perc. of chord from aerodynamic center)
* `e`: position of the elastic center (perc. of chord from aerodynamic center)
* `m`: linear density of the wing section
* `ŷ`: y vector, perpendicular to the axis of the wing,
with which the flexural resistances `EIy` and `EIz` are calculated.
Defaults to x axis
* `EA`: axial resistance for the beam
* `EIy`: flexural resistance in the axis indicated by `ŷ`
* `EIz`: flexural resistance in the axis perpendicular to the wing and 
the `ŷ` axis
* `GJ`: torsional resistance
* `Ixx`: inertia of the wing section around its elastic center (per unit length)
"""
function add_wing_strip!(
	acft::Aircraft{Fg},
	ipt1::Int64,
	ipt2::Int64,
	c::Fg;
	incidence::Fg = 0.0,
	afl::Airfoil{Fg} = flat_plate,
	xCG::Fg = 0.25,
	e::Fg = 0.0,
	m::Fg = 1e-3,
	ŷ::Vector{Fg} = [1.0, 0.0, 0.0],
	EA::Fg = Fmax,
	EIy::Fg = Fmax,
	EIz::Fg = Fmax,
	GJ::Fg = Fmax,
	Ixx::Fg = 1e-3
) where Fg

	_afl = Airfoil{Fg}(
		afl.CL0 + afl.CLα * deg2rad(incidence),
		afl.CLα,
		afl.Cm0 + afl.Cmα * deg2rad(incidence),
		afl.Cmα,
		copy(afl.polar),
		afl.CLmin,
		afl.CLmax
	)

	b = add_beam!(
		acft,
		ipt1,
		ipt2;
		m = m,
		Ixx = Ixx,
		ŷ = ŷ,
		EA = EA,
		EIy = EIy,
		EIz = EIz,
		GJ = GJ
	)

	sect = WingStrip{Fg}(
		ipt1, ipt2, length(acft.beams),
		b.Mtosys,
		b.Mtouni,
		norm(acft.points[:, ipt2] .- acft.points[:, ipt1]),
		c,
		xCG,
		e,
		_afl
	)

	push!(acft.wing_strips, sect)

	sect

end

"""
```
	function add_fuselage_cut!(
		acft::Aircraft{Fg},
		ipt1::Int64,
		ipt2::Int64,
		A1::Fg,
		A2::Fg;
		m::Fg = 1e-3,
		ŷ::Vector{Fg} = [0.0, 1.0, 0.0],
		EA::Fg = Fmax,
		EIy::Fg = Fmax,
		EIz::Fg = Fmax,
		GJ::Fg = Fmax,
		Ixx::Fg = 1e-3,
		Cf::Fg = 0.0
	) where Fg
```

Add fuselage section to aerodynamic model

* `acft`: aircraft struct
* `ipt1`: index for point number 1
* `ipt2`: index for point number 2
* `A1`: sectional area at point 1
* `A2`: sectional area at point 2
* `m`: linear density of the wing section
* `ŷ`: y vector, perpendicular to the axis of the fuselage,
with which the flexural resistances `EIy` and `EIz` are calculated.
Defaults to y axis
* `EA`: axial resistance for the beam
* `EIy`: flexural resistance in the axis indicated by `ŷ`
* `EIz`: flexural resistance in the axis perpendicular to the fuselage and 
the `ŷ` axis
* `GJ`: torsional resistance
* `Ixx`: inertia of the fuselage section around its axis (per unit length)
* `Cf`: friction coefficient
"""
function add_fuselage_cut!(
	acft::Aircraft{Fg},
	ipt1::Int64,
	ipt2::Int64,
	A1::Fg,
	A2::Fg;
	m::Fg = 1e-3,
	ŷ::Vector{Fg} = [0.0, 1.0, 0.0],
	EA::Fg = Fmax,
	EIy::Fg = Fmax,
	EIz::Fg = Fmax,
	GJ::Fg = Fmax,
	Ixx::Fg = 1e-3,
	Cf::Fg = 0.0
) where Fg

	b = add_beam!(
		acft,
		ipt1,
		ipt2;
		m = m,
		ŷ = ŷ,
		EA = EA,
		EIy = EIy,
		EIz = EIz,
		GJ = GJ,
		Ixx = Ixx
	)

	L = norm(acft.points[:, ipt1] .- acft.points[:, ipt2])

	cut = FuselageCut{Fg}(
		ipt1,
		ipt2,
		length(acft.beams),
		b.Mtosys,
		b.Mtouni,
		A1,
		A2,
		L,
		Cf
	)

	push!(acft.fuselage_cuts, cut)

	cut

end

"""
```
	function add_point!(
		acft::Aircraft{Fg},
		point::AbstractVector,
		tol::Real = 1e-10
	) where Fg
```

Add a point to an aircraft.

* `acft`: aircraft at hand
* `point`: point at hand (vector of length 3)
* `tol`: tolerance for point merging
"""
function add_point!(
	acft::Aircraft{Fg},
	point::AbstractVector,
	tol::Real = 1e-10
) where Fg

	if tol != 0.0 && size(acft.points, 2) != 0
		indmin, distmin = let dp = acft.points .- point
			dists = sqrt.(dp[1, :] .^ 2 .+ dp[2, :] .^ 2 .+ dp[3, :] .^ 2)

			indmin = argmin(dists)

			indmin, dists[indmin]
		end

		if distmin < tol
			return indmin
		end
	end

	acft.points = [acft.points point]
	push!(acft.masses, LumpedMass{Fg}(zero(Fg), zeros(Fg, 3, 3)))

	size(acft.points, 2)

end

"""
```
	function interpolate_wing!(
		acft::Aircraft,
		dist::AbstractVector,
		point1::AbstractVector,
		point2::AbstractVector,
		tol::Real = 1e-10;
		args1 = [],
		args2 = [],
		kwargs1 = Dict(),
		kwargs2 = Dict()
	)
```

Add a wing by interpolating properties between two points

* `acft`: the aircraft at hand
* `dist`: array with points in a 0 to 1 interval, used 
	as base for interpolation
* `point1`: coordinates for first point
* `point2`: coordinates for second point
* `tol`: optional argument for point merging tolerance
* `args1`: arguments for first point
* `args2`: arguments for second point
* `kwargs1`: keyword arguments for first point
* `kwargs2`: keyword arguments for second point

Returns: array of point indices and array of wing section indices
"""
function interpolate_wing!(
	acft::Aircraft,
	dist::AbstractVector,
	point1::AbstractVector,
	point2::AbstractVector,
	tol::Real = 1e-10;
	args1 = [],
	args2 = [],
	kwargs1 = Dict(),
	kwargs2 = Dict()
)

	point_inds = [
		add_point!(
			acft,
			(@. point1 * (1.0 - η) + point2 * η),
			tol
		) for η in dist
	]

	distmed = (dist[1:(end - 1)] .+ dist[2:end]) ./ 2

	sect_inds = [
		begin
			args = [
				(a1 * (1.0 - η) + a2 * η) for (a1, a2) in zip(args1, args2)
			]
			kwargs = [
				k1 => (v1 * (1.0 - η) + kwargs2[k1] * η) for (
					k1, v1
				) in kwargs1
			]

			add_wing_strip!(
				acft,
				point_inds[i],
				point_inds[i + 1],
				args...;
				kwargs...
			)

			length(acft.wing_strips)
		end for (i, η) in enumerate(distmed)
	]

	point_inds, sect_inds

end

Base.length(afl::Airfoil) = 1
Base.iterate(afl::Airfoil) = (afl, false)
Base.iterate(afl::Airfoil, state::Bool) = nothing

"""
```
	function interpolate_fuselage!(
		acft::Aircraft,
		dist::AbstractVector,
		point1::AbstractVector,
		point2::AbstractVector,
		A1::Real,
		A2::Real,
		tol::Real = 1e-10;
		args1 = [],
		args2 = [],
		kwargs1 = Dict(),
		kwargs2 = Dict()
	)
```

Add a fuselage by interpolating properties between two points

* `acft`: the aircraft at hand
* `dist`: array with points in a 0 to 1 interval, used 
	as base for interpolation
* `point1`: coordinates for first point
* `point2`: coordinates for second point
* `A1`: cross-sectional area at first section
* `A2`: cross-sectional area at second section
* `tol`: optional argument for point merging tolerance
* `args1`: arguments for first point
* `args2`: arguments for second point
* `kwargs1`: keyword arguments for first point
* `kwargs2`: keyword arguments for second point

Returns: array of point indices and array of wing section indices
"""
function interpolate_fuselage!(
	acft::Aircraft,
	dist::AbstractVector,
	point1::AbstractVector,
	point2::AbstractVector,
	A1::Real,
	A2::Real,
	tol::Real = 1e-10;
	args1 = [],
	args2 = [],
	kwargs1 = Dict(),
	kwargs2 = Dict()
)

	point_inds = [
		add_point!(
			acft,
			(@. point1 * (1.0 - η) + point2 * η),
			tol
		) for η in dist
	]

	distmed = (dist[1:(end - 1)] .+ dist[2:end]) ./ 2

	sect_inds = [
		begin
			args = [
				(@. a1 * (1.0 - η) + a2 * η) for (a1, a2) in zip(args1, args2)
			]
			kwargs = [
				k1 => (@. v1 * (1.0 - η) + kwargs2[k1] * η) for (
					k1, v1
				) in kwargs1
			]

			add_fuselage_cut!(
				acft,
				point_inds[i],
				point_inds[i + 1],
				A1 * (1.0 - dist[i]) + A2 * dist[i],
				A1 * (1.0 - dist[i + 1]) + A2 * dist[i + 1],
				args...;
				kwargs...
			)

			length(acft.fuselage_cuts)
		end for (i, η) in enumerate(distmed)
	]

	point_inds, sect_inds

end

const _default_flap_deflection = 0.698
const _default_slat_deflection = 0.698

"""
```
	function plain_flap(
		afl::Airfoil{Fg}, 
		cperc::Real, 
		deflection::Real;
		sweep::Real = 0.0
	) where Fg
```

Convert an airfoil to its approximate, plain flap equivalent, 
according to Raymer's "Aircraft Design: a Conceptual Approach"

* `afl`: original airfoil
* `cperc`: percentage of chord
* `deflection`: deflection in degrees
* `sweep`: sweep in degrees (optional)
"""
function plain_flap(
	afl::Airfoil{Fg}, 
	cperc::Real, 
	deflection::Real;
	sweep::Real = 0.0
) where Fg

	deflection = deg2rad(deflection)

	Δα = cperc * deflection

	dCLm = 0.9 * deflection * cosd(sweep) / _default_flap_deflection

	Airfoil{Fg}(
		afl.CL0 + Δα * afl.CLα,
		afl.CLα,
		afl.Cm0 + Δα * afl.Cmα,
		afl.Cmα,
		[
			afl.polar[1],
			afl.polar[2],
			afl.polar[3] + rad2deg(deflection) * 0.0023
		],
		afl.CLmin + dCLm,
		afl.CLmax + dCLm
	)

end

"""
```
	function slotted_flap(
		afl::Airfoil{Fg}, 
		cperc::Real, 
		deflection::Real;
		sweep::Real = 0.0,
		n_slots::Int64 = 1
	) where Fg
```

Convert an airfoil to its approximate, slotted flap equivalent, 
according to Raymer's "Aircraft Design: a Conceptual Approach"

* `afl`: original airfoil
* `cperc`: percentage of chord
* `deflection`: deflection in degrees
* `sweep`: sweep in degrees (optional)
* `n_slots`: number of slots. Ideally between 1 and 3. Defaults
	to 1
"""
function slotted_flap(
	afl::Airfoil{Fg}, 
	cperc::Real, 
	deflection::Real;
	sweep::Real = 0.0,
	n_slots::Int64 = 1
) where Fg

	deflection = deg2rad(deflection)

	Δα = cperc * deflection

	mult = 1.0 + 0.3 * n_slots

	dCLm = mult * deflection * cosd(sweep) * cperc / _default_flap_deflection

	Airfoil{Fg}(
		afl.CL0 + Δα * afl.CLα,
		afl.CLα,
		afl.Cm0 + Δα * afl.Cmα,
		afl.Cmα,
		[
			afl.polar[1],
			afl.polar[2],
			afl.polar[3] + rad2deg(deflection) * 0.0023
		],
		afl.CLmin + dCLm,
		afl.CLmax + dCLm
	)

end

"""
```
	function fixed_slat(
		afl::Airfoil{Fg}, 
		cperc::Real, 
		deflection::Real;
		sweep::Real = 0.0
	) where Fg
```

Convert an airfoil to its approximate, fixed slat equivalent, 
according to Raymer's "Aircraft Design: a Conceptual Approach"

* `afl`: original airfoil
* `cperc`: percentage of chord
* `deflection`: deflection in degrees
* `sweep`: sweep in degrees (optional)
"""
function fixed_slat(
	afl::Airfoil{Fg}, 
	cperc::Real, 
	deflection::Real;
	sweep::Real = 0.0
) where Fg

	deflection = deg2rad(deflection)

	Δα = 0.0 # cperc * deflection

	dCLm = 0.2 * deflection * cosd(sweep) / _default_slat_deflection

	Airfoil{Fg}(
		afl.CL0 + Δα * afl.CLα,
		afl.CLα,
		afl.Cm0 + Δα * afl.Cmα,
		afl.Cmα,
		[
			afl.polar[1],
			afl.polar[2],
			afl.polar[3] + rad2deg(deflection) * 0.0023
		],
		afl.CLmin + dCLm,
		afl.CLmax + dCLm
	)

end

"""
```
	function slat(
		afl::Airfoil{Fg}, 
		cperc::Real, 
		deflection::Real;
		sweep::Real = 0.0
	) where Fg
```

Convert an airfoil to its approximate, slat equivalent, 
according to Raymer's "Aircraft Design: a Conceptual Approach"

* `afl`: original airfoil
* `cperc`: percentage of chord
* `deflection`: deflection in degrees
* `sweep`: sweep in degrees (optional)
"""
function slat(
	afl::Airfoil{Fg}, 
	cperc::Real, 
	deflection::Real;
	sweep::Real = 0.0
) where Fg

	deflection = deg2rad(deflection)

	Δα = 0.0 # cperc * deflection

	dCLm = 0.4 * cperc * deflection * cosd(sweep) / _default_slat_deflection

	Airfoil{Fg}(
		afl.CL0 + Δα * afl.CLα,
		afl.CLα,
		afl.Cm0 + Δα * afl.Cmα,
		afl.Cmα,
		[
			afl.polar[1],
			afl.polar[2],
			afl.polar[3] + rad2deg(deflection) * 0.0023
		],
		afl.CLmin + dCLm,
		afl.CLmax + dCLm
	)

end

"""
```
	function Kruger_slat(
		afl::Airfoil{Fg}, 
		cperc::Real, 
		deflection::Real;
		sweep::Real = 0.0
	) where Fg
```

Convert an airfoil to its approximate, Kruger slat equivalent, 
according to Raymer's "Aircraft Design: a Conceptual Approach"

* `afl`: original airfoil
* `cperc`: percentage of chord
* `deflection`: deflection in degrees
* `sweep`: sweep in degrees (optional)
"""
function Kruger_slat(
	afl::Airfoil{Fg}, 
	cperc::Real, 
	deflection::Real;
	sweep::Real = 0.0
) where Fg

	deflection = deg2rad(deflection)

	Δα = 0.0 # cperc * deflection

	dCLm = 0.3 * deflection * cosd(sweep) / _default_slat_deflection

	Airfoil{Fg}(
		afl.CL0 + Δα * afl.CLα,
		afl.CLα,
		afl.Cm0 + Δα * afl.Cmα,
		afl.Cmα,
		[
			afl.polar[1],
			afl.polar[2],
			afl.polar[3] + rad2deg(deflection) * 0.0023
		],
		afl.CLmin + dCLm,
		afl.CLmax + dCLm
	)

end

"""
```
	function get_massic_properties(
		acft::Aircraft{Fg}
	) where Fg
```

Get mass properties for an aircraft.

Returns in format:

```
(
	I = [
		1.0 0.0 0.0;
		0.0 1.0 0.0;
		0.0 0.0 1.0 # example
	],
	CG = [0.0, 0.0, 0.0],
	M = 1.0
)
```
"""
function get_massic_properties(
	acft::Aircraft{Fg}
) where Fg

	I = zeros(3, 3)
	CG = zeros(3)
	M = 0.0

	for (i, m) in enumerate(acft.masses)
		r = acft.points[:, i]

		M += m.m
		I += m.I + m.m * [
			(r[2] ^ 2 + r[3] ^ 2) (r[1] * r[2]) (r[1] * r[3]);
			(r[1] * r[2]) (r[1] ^ 2 + r[3] ^ 2) (r[2] * r[3]);
			(r[1] * r[3]) (r[2] * r[3]) (r[1] ^ 2 + r[2] ^ 2)
		]
		CG += m.m .* r
	end

	(
		M = M,
		CG = CG / M,
		I = I
	)

end

struct WingSectInfo
	pt::AbstractVector
	args
	kwargs::Dict{Symbol, Any}
end

"""
```
	WingSectInfo(
		pt::AbstractVector,
		c::Real,
		args...;
		kwargs...
	)
```

Struct to contain wing section information.

* `pt`: vector with quarter chord location
* `c`: local chord

Pass any other arguments and keyword arguments to be passed
to `interpolate_wing!` when `interpolate!` is called.
"""
WingSectInfo(
	pt::AbstractVector,
	args...;
	kwargs...
) = WingSectInfo(
	pt,
	args,
	Dict{Symbol, Any}(
		kwargs...
	)
)

"""
```
	interpolate!(
		acft,
		distribution,
		s1::WingSectInfo,
		s2::WingSectInfo,
		tol::Real = 1e-10
	)
```

Interpolate two wing sections, as defined by `WingSectInfo`.

Returns the same as `interpolate_wing!`
"""
interpolate!(
	acft::Aircraft,
	distribution,
	s1::WingSectInfo,
	s2::WingSectInfo,
	tol::Real = 1e-10
) = interpolate_wing!(
	acft,
	(
		isa(distribution, Int) ?
		collect(LinRange(0.0, 1.0, distribution)) :
		collect(distribution)
	),
	s1.pt,
	s2.pt,
	tol;
	args1 = s1.args,
	args2 = s2.args,
	kwargs1 = s1.kwargs,
	kwargs2 = s2.kwargs
)

struct FuselageCutInfo
	pt::AbstractVector
	A::Real
	args
	kwargs::Dict{Symbol, Any}
end

"""
```
	FuselageCutInfo(
		pt::AbstractVector,
		A::Real,
		args...;
		kwargs...
	)
```

Struct to contain fuselage cut information.

* `pt`: vector with center location
* `A`: cross-sectional area

Pass any other arguments and keyword arguments to be passed
to `interpolate_fuselage!` when `interpolate!` is called.
"""
FuselageCutInfo(
	pt::AbstractVector,
	A::Real,
	args...;
	kwargs...
) = FuselageCutInfo(
	pt,
	A,
	args,
	Dict{Symbol, Any}(
		kwargs...
	)
)

"""
```
	interpolate!(
		acft,
		distribution,
		s1::FuselageCutInfo,
		s2::FuselageCutInfo,
		tol::Real = 1e-10
	)
```

Interpolate two wing sections, as defined by `FuselageCutInfo`.

Returns the same as `interpolate_fuselage!`
"""
interpolate!(
	acft::Aircraft,
	distribution,
	s1::FuselageCutInfo,
	s2::FuselageCutInfo,
	tol::Real = 1e-10
) = interpolate_fuselage!(
	acft,
	(
		isa(distribution, Int) ?
		collect(LinRange(0.0, 1.0, distribution)) :
		collect(distribution)
	),
	s1.pt,
	s2.pt,
	s1.A,
	s2.A,
	tol;
	args1 = s1.args,
	args2 = s2.args,
	kwargs1 = s1.kwargs,
	kwargs2 = s2.kwargs
)
