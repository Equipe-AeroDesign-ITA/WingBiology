export get_eigenmodes, assumed_modes_solve

"""
    ```
    function get_eigenmodes(
        aircraft::Aircraft{Fg},
        n::Int64 = 6,
        MinvK::Union{AbstractMatrix, Nothing} = nothing;
        fixed_points::Vector{Fi} = Int64[],
        kwargs... # keyword arguments for Arpack
    ) where {Fg, Fi <: Int}
    ```

Get the `n` first eigenmodes for an aircraft,
in order of increasing frequency

Returns a matrix in which each column is an eigenmode, 
and a vector of corresponding frequencies
"""
function get_eigenmodes(
    aircraft::Aircraft{Fg},
    n::Int64 = 10,
    MinvK::Union{AbstractMatrix, Nothing} = nothing;
    fixed_points::Vector{Fi} = Int64[],
    kwargs... # keyword arguments for Arpack
) where {Fg, Fi <: Int}

    MinvK = (
        isnothing(MinvK) ?
        get_elastic_matrix(aircraft; fixed_points = fixed_points) :
        MinvK
    )

    N = ndofs(aircraft) ÷ 2
    ϕ, ω = let A = MinvK[(N + 1):end, 1:N]
        for f in fixed_points
            for idof = (6 * (f - 1) + 1):(6 * f)
                A[idof, :] .= 0.0

                A[idof, idof] = - Fmax
            end
        end

        dropzeros!(A)

        ω, ϕ = eigs(A; nev = n, which = :SM, kwargs...)

        ϕ = [
            real.(ϕ);
            zeros(Fg, N, n)
        ]

        ϕ, sqrt.(abs.(ω))
    end

    ϕ, ω

end

"""
```
	function assumed_modes_solve(
		acft::Aircraft{Fg},
		x::AbstractVector,
		U∞::Fg;
		n_eig::Int64 = 6,
        transient_aerodynamics::Bool = true,
		kwargs...
	) where Fg
```

Solve aircraft aerostructural problem using Krylov subspace,
Proper Orthogonal Decomposition method

* `acft`: aircraft object
* `q`: variables of interest
* `U∞`: freestream velocity
* `n_eig`: number of eigenvalues to consider
* `transient_aerodynamics`: whether to use transient 
aerodynamics for calculations
* `kwargs`: keyword arguments for `state_space`
"""
function assumed_modes_solve(
    acft::Aircraft{Fg},
    x::AbstractVector,
    U∞::Fg;
    n_eig::Int64 = 6,
    transient_aerodynamics::Bool = true,
    kwargs...
) where {Fg}

    fpoints = Int64[]
    for k in kwargs
        if k[1] == :fixed_points
            fpoints = k[2]
        end
    end

    ϕ, _ = get_eigenmodes(acft, n_eig; fixed_points = fpoints)

    r, fcon = state_space(
        acft,
        x,
        U∞;
        kwargs...
    )

    f = (q; Γt = nothing) -> begin
        qd, fcon = state_space(
            acft,
            q,
            U∞;
            ∂Γ!∂t = Γt,
            kwargs...
        )

        qd, fcon.Γ
    end
    df = qd -> begin
        qduals = [
            ForwardDiff.Dual(qq, qqd) for (qq, qqd) in zip(x, qd)
        ]

        qd_d, Γ_d = f(qduals)

        [
            v.partials[1] for v in qd_d
        ], [
            v.partials[1] for v in Γ_d
        ]
    end

    A = Matrix{Fg}(undef, ndofs(acft), n_eig)
    JV = similar(A)

    dΓ!dx = Matrix{Fg}(undef, length(acft.wing_strips), n_eig)

    for i = 1:n_eig
        A[:, i] .= ϕ[:, i]

        #=
        for j = 1:(i-1)
            A[:, i] .-= (A[:, i] ⋅ A[:, j]) .* A[:, j]
        end
        =#

        A[:, i] ./= norm(A[:, i])

        qd, Γ = df(A[:, i])

        JV[:, i] .= qd
        dΓ!dx[:, i] .= Γ
    end

    for i = 1:n_eig
        v = [
            zeros(Fg, ndofs(acft) ÷ 2);
            A[1:(ndofs(acft) ÷ 2), i]
        ]

        #=
        for j = 1:size(A, 2)
            v .-= (v ⋅ A[:, j]) .* A[:, j]
        end
        =#

        v ./= norm(v)

        qd, Γ = df(v)

        A = [A v]
        JV = [JV qd]
        dΓ!dx = [dΓ!dx Γ]
    end

    D = (
        transient_aerodynamics ?
        hcat(
            [
                f(x; Γt = dΓ!dx[:, i])[1] .- r for i = 1:size(dΓ!dx, 2)
            ]...
        ) :
        nothing
    )

    #=
    A, JV, D = map(
        mat -> mat[:, 1:n_eig] .+ mat[:, (n_eig + 1):end] .* 1im,
        [
            A, JV, D
        ]
    )
    =#

    eig = let M = (
        transient_aerodynamics ?
        (A - D) :
        A
    ), K = JV
        eigen(
            _linv_reg(M) * K
        )
    end

    A, (
        ϕ = eig.vectors,
        λ = eig.values,
        JA = JV,
        B = D
    )

end
