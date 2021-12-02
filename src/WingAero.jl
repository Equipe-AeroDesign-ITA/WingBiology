geps = 1e-5

"""
Function to define inf. velocity of a horseshoe vortex

* `a`: first kin (relative to collocation point: `colpt .- kin`)
* `b`: second kin (relative to collocation point: `colpt .- kin`)
* `x`: wind direction (unit vector for adimensionalized results, equal to vectorial
freestream velocity if dimensional results are desired)

* return: the influence velocity over the collocation point
"""
function hshoe_vinf(a::AbstractVector, b::AbstractVector, x::AbstractVector)
    na=norm(a)
    nb=norm(b)
    ab=dot(a, b)

    d1=(na*nb+ab)
    if d1>-geps && d1<=0.0
        d1=-geps
    elseif d1<geps && d1>0.0
        d1=geps
    end

    d2=(na-dot(a, x))
    if d2>-geps && d2<=0.0
        d2=-geps
    elseif d2<geps && d2>0.0
        d2=geps
    end

    d3=(nb-dot(b, x))
    if d3>-geps && d3<=0.0
        d3=-geps
    elseif d3<geps && d3>0.0
        d3=geps
    end

    vbound = (cross(a, b)/d1)*(1.0/na+1.0/nb)
    vinf = cross(a, x)/(na*d2).-cross(b, x)/(nb*d3)

    vinf./=(4*pi)
    vbound./=(4*pi)

    return vinf, vbound
end

"""
Function to define inf. velocity of a finite vortex line segment

* `a`: first kin (relative to collocation point: `colpt .- kin`)
* `b`: second kin (relative to collocation point: `colpt .- kin`)

* return: the influence velocity over the collocation point
"""
function vline_segment_vinf(a::AbstractVector, b::AbstractVector)
    na=norm(a)
    nb=norm(b)
    ab=dot(a, b)

    d1=(na*nb+ab)
    if d1>-geps && d1<=0.0
        d1=-geps
    elseif d1<geps && d1>0.0
        d1=geps
    end

    vinf=(cross(a, b)/d1)*(1.0/na+1.0/nb)

    vinf/=(4*pi)

    return vinf
end

"""
Build linear aerodynamic problem for an aircraft
"""
function linear_problem(acft, kins, cols, norms, us, xhat, PGβ)

    b = [
        begin
            nu = norm(u)

            (
                - nu * strip.afl.CL0 - u ⋅ n * strip.afl.CLα
            )
        end * PGβ for (i, (u, n, strip)) in enumerate(zip(us, norms, acft.wing_strips))
    ]

    Fg = eltype(b)

    As = [
        Matrix{Fg}(undef, length(b), length(b)) for i = 1:3
    ]

    A = [
        begin
            vinfl, vbound = hshoe_vinf(cols[i] .- kins[j][1], cols[i] .- kins[j][2], xhat)

            for (id, A) in enumerate(As)
                A[i, j] = vinfl[id] + vbound[id]
            end

            C = acft.wing_strips[i].afl.CLα * PGβ * (
                (vinfl .+ vbound) ⋅ 
                norms[i]
            )

            if i == j
                C -= 2.0 / acft.wing_strips[i].c
            end

            C
        end for i = 1:length(acft.wing_strips), j = 1:length(acft.wing_strips)
    ]

    A, b, As

end

"""
Add downwash to velocities at collocation points
"""
function add_downwash!(us, Γs, As)

    dVs = [
        A * Γs for A in As
    ]

    for (id, dV) in enumerate(dVs)
        for (i, u) in enumerate(us)
            u[id] += dV[i]
        end
    end

end
