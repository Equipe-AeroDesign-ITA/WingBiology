export solve!, get_eigen, get_Jacobian, get_elastic_preconditioner, get_elastic_matrix

"""
```
	function solve!(
		acft::Aircraft{Fg},
		x::AbstractVector,
		U∞::Fg;
		n_iter::Int64 = 1,
		verbose::Bool = false,
		kwargs...
	) where Fg
```

Solve aircraft aerostructural problem using Newton method

* `acft`: aircraft object
* `q`: variables of interest
* `U∞`: freestream velocity
* `n_iter`: number of iterations. Defaults to 1
* `verbose`: display residual at each iteration
* `kwargs`: keyword arguments for `state_space`
"""
function solve!(
	acft::Aircraft{Fg},
	x::AbstractVector,
	U∞::Fg;
	n_iter::Int64 = 1,
	verbose::Bool = false,
	kwargs...
) where Fg

	f = x -> begin
		qd, _ = state_space(
			acft,
			x,
			U∞;
			kwargs...
		)

		qd
	end

	for nit = 1:n_iter
		r = f(x)

		if verbose
			println("Iteration $nit - residual norm $(norm(r))")
		end

		J = ForwardDiff.jacobian(
			f,
			x
		)

		x .-= J \ r
	end

end

"""
Get eigenvalues and eigenvectors of the state-space description
of the aircraft

* `acft`: aircraft object
* `q`: variables of interest
* `U∞`: freestream velocity
* `kwargs`: keyword arguments for `state_space`
"""
function get_eigen(
	acft::Aircraft{Fg},
	x::AbstractVector,
	U∞::Fg;
	kwargs...
) where Fg

	f = x -> begin
		qd, _ = state_space(
			acft,
			x,
			U∞;
			kwargs...
		)

		qd
	end

	J = ForwardDiff.jacobian(
		f,
		x
	)

	eigen(J)

end

"""
Get Jacobian of the state-space description
of the aircraft

* `acft`: aircraft object
* `q`: variables of interest
* `U∞`: freestream velocity
* `kwargs`: keyword arguments for `state_space`
"""
function get_Jacobian(
	acft::Aircraft{Fg},
	x::AbstractVector,
	U∞::Fg;
	kwargs...
) where Fg

	f = x -> begin
		qd, _ = state_space(
			acft,
			x,
			U∞;
			kwargs...
		)

		qd
	end

	J = q -> ForwardDiff.jacobian(
		f,
		x .+ q
	)

	J

end

"""
	```
	function get_elastic_matrix(
		acft::Aircraft{Fg};
		fixed_points::Vector{Int64} = Int64[]
	) where Fg
	```

Get matrix `inv(M) * K` in sparse format
"""
function get_elastic_matrix(
	acft::Aircraft{Fg};
	fixed_points::Vector{Int64} = Int64[]
) where Fg

	q = get_state(acft, 1.0)
	
	f = (y, x) -> begin
		qd, _ = state_space(
			acft,
			x,
			1.0;
			aerodynamic_forces = false,
			structural_damping = false,
			fixed_points = fixed_points
		)

		y .= qd

		y
	end

	HD = ndofs(acft) ÷ 2

	n_elements = 0

	for b in acft.beams
		for j = (6 * (b.ipt1 - 1) + 1):(6 * b.ipt1)
			for i = (6 * (b.ipt1 - 1) + 1 + HD):(6 * b.ipt1 + HD)
				n_elements += 1
			end

			for i = (6 * (b.ipt2 - 1) + 1 + HD):(6 * b.ipt2 + HD)
				n_elements += 1
			end
		end

		for j = (6 * (b.ipt2 - 1) + 1):(6 * b.ipt2)
			for i = (6 * (b.ipt1 - 1) + 1 + HD):(6 * b.ipt1 + HD)
				n_elements += 1
			end

			for i = (6 * (b.ipt2 - 1) + 1 + HD):(6 * b.ipt2 + HD)
				n_elements += 1
			end
		end
	end

	for i = 1:HD
		n_elements += 1
	end
	
	for f in fixed_points
		n_elements += 12
	end

	data = ones(Fg, n_elements)
	n_rows = Vector{Int64}(undef, n_elements)
	n_cols = Vector{Int64}(undef, n_elements)

	n_elements = 0

	for b in acft.beams
		for j = (6 * (b.ipt1 - 1) + 1):(6 * b.ipt1)
			for i = (6 * (b.ipt1 - 1) + 1 + HD):(6 * b.ipt1 + HD)
				n_elements += 1

				n_rows[n_elements] = i
				n_cols[n_elements] = j
			end

			for i = (6 * (b.ipt2 - 1) + 1 + HD):(6 * b.ipt2 + HD)
				n_elements += 1

				n_rows[n_elements] = i
				n_cols[n_elements] = j
			end
		end

		for j = (6 * (b.ipt2 - 1) + 1):(6 * b.ipt2)
			for i = (6 * (b.ipt1 - 1) + 1 + HD):(6 * b.ipt1 + HD)
				n_elements += 1

				n_rows[n_elements] = i
				n_cols[n_elements] = j
			end

			for i = (6 * (b.ipt2 - 1) + 1 + HD):(6 * b.ipt2 + HD)
				n_elements += 1

				n_rows[n_elements] = i
				n_cols[n_elements] = j
			end
		end
	end

	for i = 1:HD
		n_elements += 1

		n_rows[n_elements] = i
		n_cols[n_elements] = i + HD
	end

	for f in fixed_points
		for i = (6 * (f - 1) + 1):(6 * f)
			n_elements += 1

			n_rows[n_elements] = i
			n_cols[n_elements] = i

			n_elements += 1

			n_rows[n_elements] = i + HD
			n_cols[n_elements] = i + HD
		end
	end

	jac = sparse(n_rows, n_cols, data)

	forwarddiff_color_jacobian!(
		jac,
		f,
		q;
		colorvec = matrix_colors(jac)
	)

	jac

end	

"""
	```
	function get_elastic_preconditioner(
		acft::Aircraft{Fg};
		fixed_points::Vector{Int64} = Int64[]
	) where Fg
	```

Get elastic matrix `inv(K) * M` preconditioner (LU decomposition
of Jacobian when disregarded elastic forces)
"""
function get_elastic_preconditioner(
	acft::Aircraft{Fg};
	fixed_points::Vector{Int64} = Int64[]
) where Fg

	jac = get_elastic_matrix(acft; fixed_points = fixed_points)

	lu(jac)

end
