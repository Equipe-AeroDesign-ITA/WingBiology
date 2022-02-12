using WingBiology

using JSON

results = Dict()


for g in (4e-3, 7e-3, 1e-2)
    c = 1.8288
    b = 2 * 6.096

    e = 0.33 - 0.25
    xCG = 0.43 - 0.25

    m = 35.71
    Ixx = 8.64

    GJ = 0.99e6
    EIxx = 9.77e6

    ρ = 1.020
    α = 0.05

    acft = Aircraft()

    ls = WingSectInfo(
        [0.0, - b / 2, 0.0],
        c;
        e = e,
        xCG = xCG,
        m = m,
        Ixx = Ixx,
        EIy = EIxx,
        GJ = GJ,
        g = g
    )
    cs = WingSectInfo(
        [0.0, 0.0, 0.0],
        c;
        e = e,
        xCG = xCG,
        m = m,
        Ixx = Ixx,
        EIy = EIxx,
        GJ = GJ,
        g = g
    )
    rs = WingSectInfo(
        [0.0, b / 2, 0.0],
        c;
        e = e,
        xCG = xCG,
        m = m,
        Ixx = Ixx,
        EIy = EIxx,
        GJ = GJ,
        g = g
    )

    lwing_pts, lwing_sects = interpolate!(acft, 30, ls, cs)
    rwing_pts, rwing_sects = interpolate!(acft, 30, cs, rs)

    add_link!(acft, rwing_pts[1])


    gres = Dict()
    results[
        "$g"
    ] = gres

    for U∞ = 100.0:5.0:220.0
        @show g, U∞

        vres = Dict()
        gres[
            "$U∞"
        ] = vres

        q = get_state(acft, U∞; α = α)

        #=
        Arnoldi_solve!(acft, q, U∞; α = α, ρ = ρ)

        plot_aircraft(acft, q)
        =#

        A, info = assumed_modes_solve(acft, q, U∞; ρ = ρ, α = α, n_eig = 4)

        λs = info.λ

        #=
        for (i, λ) in enumerate(λs)
            @show i, λ
        end
        =#

        vres["real_parts"] = real.(λs)
        vres["imaginary_parts"] = imag.(λs)
    end
end

let fobj = open("goland_results.json", "w")
    JSON.print(fobj, results)

    close(fobj)
end
