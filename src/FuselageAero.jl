function _source_f(l::AbstractVector, r::AbstractVector, ϵ::Real)

    L = dot(r, l)

    non_normal = r - L * l

    norm_non_normal = sqrt(norm(non_normal) ^ 2 + ϵ ^ 2)
    nr = sqrt(norm(r) ^ 2 + ϵ ^ 2)

    (l - non_normal * (L / norm_non_normal ^ 2)) / nr

end

"""
Source line influence coefficient
"""
function source_line(colpt::AbstractVector, r1::AbstractVector, r2::AbstractVector, ϵ::Real)

    l = r2 - r1
    nl = norm(l) + 1e-10
    @. l /= nl

    (_source_f(l, colpt - r2, ϵ) - _source_f(l, colpt - r1, ϵ)) / (4 * pi)

end

function _gradPhi(r⃗::AbstractVector, ξ⃗::AbstractVector, η⃗::AbstractVector, ϵ::Real)

    rxi = dot(r⃗, ξ⃗)
    reta = dot(r⃗, η⃗)

    r = sqrt(norm(r⃗) ^ 2 + ϵ ^ 2)

    den = (r ^ 2 - rxi ^ 2 + ϵ ^ 2) * r

    (
        (rxi * η⃗ + reta * ξ⃗) * den - (den * r⃗ / r ^ 2 + 2 * (r⃗ - rxi * ξ⃗) * r) * rxi * reta
    ) / den ^ 2

end

"""
Doublet line influence coefficient
"""
function doublet_line(colpt::AbstractVector, r1::AbstractVector, r2::AbstractVector, η⃗::AbstractVector, ϵ::Real)

    l = r2 - r1
    nl = norm(l) + 1e-10
    @. l /= nl

    (_gradPhi(colpt - r2, l, η⃗, ϵ) - _gradPhi(colpt - r1, l, η⃗, ϵ)) / (4 * pi)

end

"""
Point doublet influence coefficient
"""
function point_doublet_infl(colpt::AbstractVector, pt::AbstractVector, ξ::AbstractVector, ϵ::Real)

    r = colpt .- pt
    nr = sqrt(norm(r) ^ 2 + ϵ ^ 2)

    z = ξ ⋅ r
    nprojr = r .- z .* ξ

    (3 * z) .* nprojr ./ (nr ^ 5 * 4 * π) .- ξ .* ((nr ^ 2 - 3 * z ^ 2) / (4 * π * nr ^ 5))

end

"""
Point source influence coefficient
"""
function point_source_infl(colpt::AbstractVector, pt::AbstractVector, ϵ::Real)

    r = colpt .- pt

    r ./ (4 * π * sqrt(norm(r) + ϵ ^ 2) ^ 3)

end
