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

    ϕ, ω = get_eigenmodes(acft, 1)
    
    @assert isapprox(ω[1], 2.0; rtol = 1e-2)

    true
end
