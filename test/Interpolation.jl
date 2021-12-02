@assert begin
    acft = Aircraft()

    dist = collect(LinRange(0.0, 1.0, 50))

    s1 = WingSectInfo(
        [1.0, -5.0, 0.0],
        1.0;
        xCG = 0.5
    )
    s2 = WingSectInfo(
        zeros(3),
        2.0;
        xCG = 0.25
    )
    s3 = WingSectInfo(
        [1.0, 5.0, 0.0],
        1.0;
        xCG = 0.5
    )
    
    interpolate!(acft, dist, s1, s2)
    interpolate!(acft, dist, s2, s3)

    # plot_aircraft(acft)
    
    true
end
