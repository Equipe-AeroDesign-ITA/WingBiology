@assert begin
    acft = Aircraft()

    U∞ = 100.0
    ρ = 1.2

    b = 2.0
    c = 0.5

    xCG = 0.1
    e = 0.2

    θc = 1e-2
    hc = 1e-2
    θdc = 1.0
    hdc = 2.0

    EI = 1e4
    GJ = 1e4

    m = 2.0
    Ixx = m * c ^ 2 / 12

    s1 = WingSectInfo(
        [0.0, - b / 2, 0.0],
        c;
        e = e,
        xCG = xCG,
        EIy = EI,
        GJ = GJ,
        m = m,
        Ixx = Ixx
    )
    s2 = WingSectInfo(
        [0.0, 0.0, 0.0],
        c;
        e = e,
        xCG = xCG,
        EIy = EI,
        GJ = GJ,
        m = m,
        Ixx = Ixx
    )
    s3 = WingSectInfo(
        [0.0, b / 2, 0.0],
        c;
        e = e,
        xCG = xCG,
        EIy = EI,
        GJ = GJ,
        m = m,
        Ixx = Ixx
    )

    interpolate!(acft, [0.0, 1.0], s1, s2)
    interpolate!(acft, [0.0, 1.0], s2, s3)

    q = get_state(acft, U∞)

    set_θy(q, 2, θc)
    set_w(q, 2, hc)
    set_̇θy(q, 2, θdc)
    set_̇w(q, 2, hdc)

    qd, fcon = state_space(acft, q, U∞; ρ = ρ, structural_damping = false)

    α_AC = get_̇θy(qd, 2)
    a_AC = get_̇w(qd, 2)

    α_CG_num = α_AC
    a_CG_num = a_AC - α_AC * c * xCG

    # analytical

    m = (b / 2) * m
    Ixx = (b / 2) * Ixx

    L = b / 2

    Kh = 2 * 12 * EI / L ^ 3
    Kθ = 2 * GJ / L

    S = b * c / 2

    pdyn = ρ * (U∞ ^ 2 + hdc ^ 2) / 2
    Laed = ((fcon.sectional.CL[1] + fcon.sectional.CL[2]) / 2) * S * pdyn
    Maed = ((fcon.sectional.Cm[1] + fcon.sectional.Cm[2]) / 2) * S * pdyn

    h_e = - c * e * θc + hc

    M_e = - Kθ * θc
    L_e = - Kh * h_e

    M_CG = Maed + Laed * c * xCG + M_e + L_e * (xCG - e) * c
    L_CG = Laed + L_e

    a_CG = L_CG / m
    α_CG = M_CG / Ixx

    # Gotta make these values match!
    @assert isapprox(a_CG, a_CG_num; rtol = 1e-3)
    @assert isapprox(α_CG, α_CG_num; rtol = 1e-3) 

    # plot_aircraft(acft, q)

    true
end
