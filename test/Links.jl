@test begin
    acft = Aircraft()

    i1 = add_point!(acft, zeros(3))
    i2 = add_point!(acft, [1.0, 0.0, 0.0])

    add_beam!(
        acft, i1, i2;
        m = 1.0,
        Ixx = 1.0,
        EA = 10.0
    )

    add_link!(acft, i1)

    K = get_elastic_matrix(acft)

    true
end
