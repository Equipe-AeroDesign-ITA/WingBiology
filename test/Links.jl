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

    add_link!(acft, i1)
    add_link!(acft, i2, i3)
    add_link!(acft, i2, i4)

    add_mass!(acft, i3, 1.0)
    add_mass!(acft, i4, 1.0)

    K = get_elastic_matrix(acft)

    u = get_state(acft, 0.0)

    set_u(u, i2, 1.5)

    qd, _ = state_space(acft, u, 0.0)

    @show get_̇u(qd)

    true
end
