@test begin
    acft = Aircraft()

    i1 = add_point!(acft, zeros(3))
    i2 = add_point!(acft, [1.0, 0.0, 0.0])
    i3 = add_point!(acft, [1.0, 1e-2, 0.0])
    i4 = add_point!(acft, [1.0, -1e-2, 0.0])

    add_beam!(
        acft, i1, i2;
        m = 1.0,
        Ixx = 1.0,
        EA = 10.0
    )

    add_mass!(acft, i3, 1.0)
    add_mass!(acft, i4, 1.0)

    add_link!(acft, i1)

    u = get_state(acft, 0.0)
    set_u(u, i2, 1.5)
    qd, _ = state_space(acft, u, 0.0)

    acc_free = get_̇u(qd)

    add_link!(acft, i2, i3)
    add_link!(acft, i2, i4)

    u = get_state(acft, 0.0)
    set_u(u, i2, 1.5)
    qd, _ = state_space(acft, u, 0.0)

    acc_const = get_̇u(qd)

    ms = [
        m.m for m in acft.masses
    ]

    @assert ms ⋅ acc_const ≈ ms ⋅ acc_free

    K = get_elastic_matrix(acft)
    K = Array(K)

    Ndof2 = ndofs(acft) ÷ 2
    K = (K[(Ndof2 + 1):end, 1:Ndof2])

    eig = eigen(K)

    zero_eigs = abs.(eig.values) .< 1e-5

    zero_vals = eig.values[zero_eigs]
    zero_vecs = eig.vectors[:, zero_eigs]

    @show zero_vecs

    true
end
