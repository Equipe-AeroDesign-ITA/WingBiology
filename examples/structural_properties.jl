"""
Obtain structural properties of a carbon-fiber, tubular spar,
given internal and external diameters, local chord and 
displacement from quarter chord
"""
function structural_tube(
    Ri::Real,
    Ro::Real,
    c::Real = 1.0;
    e::Real = 0.0,
    E::Real = 57.9e9, # AeroDesign ITA: carbono em tecido
    G::Real = 2.710e9,
    ρ::Real = 1.8e3
)

    R = (Ri + Ro) / 2
    t = (Ro - Ri)

    A = π * R ^ 2

    GJ = 4 * A ^ 2 * t * G / (2 * π * R)

    A = π * (Ro ^ 2 - Ri ^ 2)

    EA = E * A

    I = π * (Ro ^ 4 - Ri ^ 4) / 4

    EI = E * I

    m = A * ρ

    Ixx = ρ * π * (Ro ^ 4 - Ri ^ 4) / 2 + c ^ 2 * e ^ 2 * m

    mx = e * m

    EA, EI, GJ, m, Ixx, mx

end

"""
Get carbon-fiber fuselage structural data, given
radius and skin thickness
"""
structural_fuselage(
    R::Real,
    t::Real;
    E::Real = 19.3e9, # AeroDesign ITA: carbono em tecido
    G::Real = (1.0 + 0.35) * 19.3e9 / 2,
    ρ::Real = 1.4e3
) = structural_tube(
    R - t, R;
    E = E,
    G = G,
    ρ = ρ
)

"""
Get approximate structural properties of a balsa wood structural skin,
given internal airfoil area (for unit chord) and local chord
"""
function structural_skin(
    Arel::Real,
    t::Real,
    c::Real = 1.0;
    G::Real = 3.4e9 / (2 * (1.0 + 0.009)),
    ρ::Real = 0.16e3
)

    GJ = 4 * Arel ^ 2 * t * G * c / 2

    m = t * 2 * ρ * c
    Ixx = m * c/ 12 + 0.25 ^ 2 * m * c

    EA = EI = 0.0

    mx = 0.25 * m

    EA, EI, GJ, m, Ixx, mx

end
