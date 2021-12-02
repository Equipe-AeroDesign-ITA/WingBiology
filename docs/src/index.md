# WingBiology

## Basic data types

```@docs
Aircraft
WingStrip
FuselageCut
PointDrag
Engine
Beam
LumpedMass
Airfoil
```

## Construction functions

```@docs
Aircraft()
Aircraft(pts)
add_point!
add_mass!
add_beam!
add_engine!
add_wing_strip!
add_fuselage_cut!
add_point_drag!
```

## Shortcuts for geometry creation

```@docs
WingSectInfo
WingSectInfo(r::AbstractVector, c::Real)
FuselageCutInfo
FuselageCutInfo(r::AbstractVector, A::Real)
interpolate!(acft::Aircraft, distro, s1::FuselageCutInfo, s2::FuselageCutInfo)
interpolate!(acft::Aircraft, distro, s1::WingSectInfo, s2::WingSectInfo)
interpolate_fuselage!(
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
interpolate_wing!(
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

## State space and force calculations

```@docs
get_state(
    acft::Aircraft{Fg},
    U∞::Fg;
    α::Fg = 0.0,
    β::Fg = 0.0,
    p::Fg = 0.0,
    q::Fg = 0.0,
    r::Fg = 0.0
) where Fg
state_space(
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

## Aerostructural solutions

**Use `Arnoldi_solve!` whenever possible!** It's far faster than the full Newton method

```
solve!(
    acft::Aircraft{Fg},
    x::AbstractVector,
    U∞::Fg;
    n_iter::Int64 = 1,
    verbose::Bool = false,
    kwargs...
) where Fg
Arnoldi_solve!(
    acft::Aircraft{Fg},
    x::AbstractVector,
    U∞::Fg;
    n_Krylov::Int64 = 6,
    first_vector::Union{AbstractVector, Nothing} = nothing,
    preconditioner::Union{LU, Nothing} = nothing,
    rtol::Real = 1e-2,
    kwargs...
) where Fg
```

## Plotting and obtaining results

```@docs
plot_aircraft(acft::Aircraft{Fg}, q::Union{AbstractVector, Nothing} = nothing) where Fg
```

For analysis of the state-space vector, see docstrings for `get_u, get_v, get_w, set_u...`

## Aeroelasticity (experimental)

```@docs
get_eigenmodes(
    aircraft::Aircraft{Fg},
    n::Int64 = 6,
    MinvK::Union{AbstractMatrix, Nothing} = nothing;
    fixed_points::Vector{Fi} = Int64[],
    kwargs... # keyword arguments for Arpack
) where {Fg, Fi <: Int}
assumed_modes_solve(
    acft::Aircraft{Fg},
    x::AbstractVector,
    U∞::Fg;
    n_eig::Int64 = 6,
    kwargs...
) where {Fg}
``` 
