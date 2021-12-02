function _source_f(l::AbstractVector, r::AbstractVector)

    L = dot(r, l)

    non_normal = r - L * l

    norm_non_normal = norm(non_normal)
    nr = norm(r)

    (l - non_normal * (L / norm_non_normal ^ 2)) / nr

end

"""
Source line influence coefficient
"""
function source_line(colpt::AbstractVector, r1::AbstractVector, r2::AbstractVector)

    l = r2 - r1
    nl = norm(l)
    @. l /= nl

    (_source_f(l, colpt - r2) - _source_f(l, colpt - r1)) / (4 * pi)

end

function _gradPhi(r⃗::AbstractVector, ξ⃗::AbstractVector, η⃗::AbstractVector)

    rxi = dot(r⃗, ξ⃗)
    reta = dot(r⃗, η⃗)

    r = norm(r⃗)

    den = (r ^ 2 - rxi ^ 2) * r

    (
        (rxi * η⃗ + reta * ξ⃗) * den - (den * r⃗ / r ^ 2 + 2 * (r⃗ - rxi * ξ⃗) * r) * rxi * reta
    ) / den ^ 2

end

"""
Doublet line influence coefficient
"""
function doublet_line(colpt::AbstractVector, r1::AbstractVector, r2::AbstractVector, η⃗::AbstractVector)

    l = r2 - r1
    nl = norm(l)
    @. l /= nl

    (_gradPhi(colpt - r2, l, η⃗) - _gradPhi(colpt - r1, l, η⃗)) / (4 * pi)

end

"""
Point doublet influence coefficient
"""
function point_doublet_infl(colpt::AbstractVector, pt::AbstractVector, ξ::AbstractVector)

    r = colpt .- pt
    nr = norm(r)

    z = ξ ⋅ r
    nprojr = r .- z .* ξ

    (3 * z) .* nprojr ./ (nr ^ 5 * 4 * π) .- ξ .* ((nr ^ 2 - 3 * z ^ 2) / (4 * π * nr ^ 5))

end

"""
Point source influence coefficient
"""
function point_source_infl(colpt::AbstractVector, pt::AbstractVector)

    r = colpt .- pt

    r ./ (4 * π * norm(r) ^ 3)

end
