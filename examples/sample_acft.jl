using WingBiology

using LinearAlgebra

using DelimitedFiles

include("airfoils.jl")
include("structural_properties.jl")

root_afl = afl350
root_afl_Arel = afl350_A

tip_afl = afl1
tip_afl_Arel = afl1_A

acft = Aircraft()

# geometry
b = 3.0
cr = 0.5
ct = 0.25
incidence = 8.0
twist = -5.0
sweep = 20.0

R_tube = 15e-3
t_tube = 3e-3

# root and tip structural properties
EAt, EIt, GJt, mt, Ixxt, mxt = structural_tube(
    R_tube, R_tube + t_tube, ct
) .+ structural_skin(
    tip_afl_Arel,
    3e-3,
    ct
)
EAr, EIr, GJr, mr, Ixxr, mxr = structural_tube(
    R_tube, R_tube + t_tube, cr
) .+ structural_skin(
    root_afl_Arel,
    3e-3,
    cr
)

s1 = WingSectInfo(
    [tand(sweep) * b / 2, - b / 2, 0.0],
    ct;
    EIy = EIt, # flexural rigidity. EIz can also be used for rigidity around the vertical axis
    GJ = GJt, # local torsional rigidity
    m = mt, # local mass per unit span
    xCG = mxt / mt, # chord perc. of CG from quarter chord
    Ixx = Ixxt, # inertia for rotation around the wing axis
    e = mxt / mt, # distance from quarter chord to elastic axis. I didn't calculate this, but you get the idea
    incidence = incidence + twist
)
s2 = WingSectInfo(
    [0.0, 0.0, 0.0],
    cr;
    EIy = EIr,
    GJ = GJr,
    m = mr,
    xCG = mxr / mr,
    Ixx = Ixxr,
    e = mxr / mr,
    incidence = incidence
)
s3 = WingSectInfo(
    [tand(sweep) * b / 2, b / 2, 0.0],
    ct;
    EIy = EIt,
    GJ = GJt,
    m = mt,
    xCG = mxt / mt,
    Ixx = Ixxt,
    e = mxt / mt,
    incidence = incidence + twist
)

interp_points, interp_wing_strips = interpolate!(acft, 20, s1, s2) # discretize with 20 strips
interpolate!(acft, 20, s2, s3) # discretize with 20 strips

fixed_points = [interp_points[end]]

U∞ = 25.0

q = get_state(acft, U∞)
q0 = copy(q)

@info "Using Arnoldi method to solve for static aerostructural solution"

A, eig = Arnoldi_solve!(acft, q, U∞; fixed_points = fixed_points)

plot_aircraft(acft, q)
qd, fcon = state_space(acft, q, U∞; fixed_points = fixed_points)

#=
@info "Obtaining assumed modes from structural matrices"

ϕ, ωs = get_eigenmodes(acft; fixed_points = fixed_points)

for (i, ω) in enumerate(ωs)
    @show i, ω

    plot_aircraft(acft, ϕ[:, i] .+ q0)
end
=#

# I haven't yet corrected the flutter prediction feature!

@info "Using assumed modes method to deduce assumed mode eigenvalues"
@warn "The flutter prediction feature is as of yet untested"

A, eig = assumed_modes_solve(acft, q, U∞; fixed_points = fixed_points, n_eig = 10)
for (i, λ) in enumerate(eig.λ)
    @show i, λ
end

writedlm("eig_basis.dat", [q0 real.(A) imag.(A)])

writedlm("assumed_eigs.dat", [real.(eig.ϕ) imag.(eig.ϕ) real.(eig.λ) imag.(eig.λ)])

@info "Performing dynamic simulation"

A, ts, qs = backward_Euler(acft, q0, U∞; fixed_points = fixed_points, transient_aerodynamics = true, dt = 1e-3, t = 0.2)

writedlm("Arnoldi_basis.dat", [q0 A])
writedlm("results.dat", [ts qs'])

#=
@info "Using full Jacobian to deduce eigenvalues"

eig = get_eigen(acft, q0, U∞) # ; fixed_points = fixed_points)

for (i, λ) in enumerate(eig.values)
    if abs(λ) < 40.0 * (2 * π)
        @show i, λ
    end
end
=#

@info "Calculating masses and inertias"

@show get_massic_properties(acft)
