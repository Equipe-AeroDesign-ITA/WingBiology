@test begin
	points = [
		0.0 1.0;
		0.0 0.0;
		0.0 0.0
	]

	acft = Aircraft(
		points
	)

	add_beam!(
		acft,
		1,
		2;
		m = 2.0,
		EA = 1e3,
		Ixx = 1.0,
		EIy = 1.0,
		EIz = 1.0,
		GJ = 1.0,
		ŷ = [0.0, 1.0, 0.0]
	)

	q = zeros(ndofs(acft))
	set_u(q, 2, 2.0)

	m = acft.masses[2].m
	qd, _ = state_space(acft, q, 0.0)
	
	@assert get_̇u(qd, 2) == - 1e3

	true
end
