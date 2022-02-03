export Arnoldi_solve!, backward_Euler

"""
```
    function Arnoldi_solve!(
        acft::Aircraft{Fg},
        x::AbstractVector,
        U∞::Fg;
        n_Krylov::Int64 = 6,
        first_vector::Union{AbstractVector,Nothing} = nothing,
        preconditioner::Union{LU,Nothing} = nothing,
        add_non_stationary::Bool = false,
        warning_tol::Real = Inf,
        kwargs...
    ) where {Fg}
```

Solve aircraft aerostructural problem using Krylov subspace,
Proper Orthogonal Decomposition method

* `acft`: aircraft object
* `q`: variables of interest
* `U∞`: freestream velocity
* `n_Krylov`: number of dimensions of Krylov subspace
* `first_vector`: first vector with which to populate the
	Krylov subspace basis. Taken as the initial residual if
	not provided
* `preconditioner`: elastic preconditioner as returned by 
	`get_elastic_preconditioner`. Calculated in place if not provided
* `add_non_stationary`: derivate the modes of the deduced subspace in respect
    to time and add them to the considered base. Necessary if the preconditioner 
    favors a stationary solution - as is the case with the default
    structural preconditioner - if a transient solution is desired
* `warning_tol`: relative tolerance for residual reduction above which a warning
    will be issued
* `kwargs`: keyword arguments for `state_space`
"""
function Arnoldi_solve!(
    acft::Aircraft{Fg},
    x::AbstractVector,
    U∞::Fg;
    n_Krylov::Int64 = 6,
    first_vector::Union{AbstractVector,Nothing} = nothing,
    preconditioner::Union{LU,Nothing} = nothing,
    add_non_stationary::Bool = false,
    warning_tol::Real = Inf,
    kwargs...
) where {Fg}

    fpoints = Int64[]
    for k in kwargs
        if k[1] == :fixed_points
            fpoints = k[2]
        end
    end

    M = (
        isnothing(preconditioner) ?
        get_elastic_preconditioner(acft; fixed_points = fpoints) :
        preconditioner
    )

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

    A = Matrix{Fg}(undef, ndofs(acft), n_Krylov)
    JV = similar(A)

    dΓ!dx = Matrix{Fg}(undef, length(acft.wing_strips), n_Krylov)

    for i = 1:n_Krylov
        A[:, i] .= (
            i == 1 ?
            (
                isnothing(first_vector) ?
                M \ r :
                first_vector
            ) :
            M \ JV[:, i-1]
        )

        for j = 1:(i-1)
            A[:, i] .-= (A[:, i] ⋅ A[:, j]) .* A[:, j]
        end

        A[:, i] ./= norm(A[:, i])

        qd, Γ = df(A[:, i])

        JV[:, i] .= qd
        dΓ!dx[:, i] .= Γ
    end

    if add_non_stationary
        for i = 1:n_Krylov
            v = [
                zeros(Fg, ndofs(acft) ÷ 2);
                A[1:(ndofs(acft) ÷ 2), i]
            ]

            for j = 1:size(A, 2)
                v .-= (v ⋅ A[:, j]) .* A[:, j]
            end

            v ./= norm(v)

            qd, Γ = df(v)

            A = [A v]
            JV = [JV qd]
            dΓ!dx = [dΓ!dx Γ]
        end
    end

    D = hcat(
        [
            f(x; Γt = dΓ!dx[:, i])[1] .- r for i = 1:size(dΓ!dx, 2)
        ]...
    )

    Ainv = _linv_reg(A)

    eig = let (M1, M2) = (
        (
            add_non_stationary ?
            Matrix{Fg}(I, size(Ainv, 1), size(Ainv, 1)) .- Ainv * D : 
            Matrix{Fg}(I, size(Ainv, 1), size(Ainv, 1))
        ),
        Ainv * JV
    )
        eigen(M2, M1)
    end

    let dα = - (JV \ r)
        x .+= A * dα

        if norm(r .+ JV * dα) > warning_tol * norm(r)
            @warn "WingBiology: desired tolerance not reached with Arnoldi method"
        end
    end

    A, (
        ϕ = eig.vectors,
        λ = eig.values,
        JA = JV,
        B = D
    )

end

"""
Function to perform a structural dynamics simulation using Backward Euler
method.

* `acft`: aircraft at hand
* `x`: state space vector
* `U∞`: freestream velocity
* `t`: total time
* `dt`: time step
* `transient_aerodynamics`: whether to use Drela's approximation
    of Theodorsen's theory for load calculation
* `kwargs`: any keyword arguments to be passed to `Arnoldi_solve!` to 
    produce a vector
"""
function backward_Euler(
    acft::Aircraft{Fg},
    x::AbstractVector,
    U∞::Real;
    t::Real = 1.0,
    dt::Real = 1e-5,
    transient_aerodynamics::Bool = true,
    kwargs...
) where {Fg}

    t = collect(0.0:dt:t)

    x0 = copy(x)
    A, info = Arnoldi_solve!(acft, x, U∞; kwargs...)

    ident = Matrix{Fg}(I, size(A, 2), size(A, 2))

    r = info.JA * (A' * (x0 .- x))

    qs = Vector{Vector{eltype(r)}}(undef, length(t))
    qs[1] = zeros(eltype(r), size(A, 2))

    x .= x0

    Δt = t[2] - t[1]

    dR = info.JA

    mat = let M = (
        transient_aerodynamics ?
        (A .- dR .* Δt .- info.B) :
        (A .- dR .* Δt)
    )
        inv(M' * M) * (M')
    end

    for i = 2:length(t)
        δx = mat * (Δt .* r)

        r .+= dR * δx

        qs[i] = qs[i - 1] .+ δx
    end

    A, t, hcat(qs...)

end
