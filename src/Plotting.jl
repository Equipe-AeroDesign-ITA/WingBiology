export plot_aircraft

_Ncirc = 20

function draw_circle(pt, R)

	θs = LinRange(0.0, 2 * π, _Ncirc)

	xs = [
		pt[1] for θ in θs
	]

	ys = [
		pt[2] + R * cos(θ) for θ in θs
	]

	zs = [
		pt[3] + R * sin(θ) for θ in θs
	]

	(xs, ys, zs)

end

"""
```
	plot_aircraft(acft::Aircraft{Fg}, q::Union{AbstractVector, Nothing} = nothing) where Fg
```

Plot an aircraft geometry with matplotlib binders. Deforms the structure if a state space vector is provided
"""
function plot_aircraft(acft::Aircraft{Fg}, q::Union{AbstractVector, Nothing} = nothing) where Fg

	points = (
		isnothing(q) ? acft.points : begin
			tA = Matrix{Fg}(undef, 3, size(acft.points, 2))

			transpose!(tA, hcat(get_u(q), get_v(q), get_w(q)))

			tA
		end
	)

	mv = maximum(abs.(points))

	_plt = pyimport("matplotlib.pyplot")
	_mpl3 = pyimport("mpl_toolkits")

	ax = _plt.axes(; projection = "3d")

	ax.scatter(points[1, :], points[2, :], points[3, :]; s = 1, color = "b")

	ax.set_xlim(- mv, mv)
	ax.set_ylim(- mv, mv)
	ax.set_zlim(- mv, mv)

	for sect in acft.fuselage_cuts
		A = (sect.A1 + sect.A2) / 2
		R = sqrt(A / π)

		xs, ys, zs = draw_circle(
			(points[:, sect.ipt1] .+ points[:, sect.ipt2]) ./ 2,
			R
		)

		ax.plot(xs, ys, zs; color = "r")
	end

	for eng in acft.engines
		R = eng.R

		xs, ys, zs = draw_circle(
			points[:, eng.ipt],
			R
		)

		ax.plot(xs, ys, zs; color = "b")
	end

	ipmass = [
		pmass.ipt for pmass in acft.point_drags
	]
	mass_pts = points[:, ipmass]
	xs = mass_pts[1, :]
	ys = mass_pts[2, :]
	zs = mass_pts[3, :]

	ax.scatter(xs, ys, zs; color = "k")

	θs = [
		begin
			θ = (
				isnothing(q) ?
				zeros(Fg, 3) :
				[
					(get_θx(q, sect.ipt1) + get_θx(q, sect.ipt2)) / 2,
					(get_θy(q, sect.ipt1) + get_θy(q, sect.ipt2)) / 2,
					(get_θz(q, sect.ipt1) + get_θz(q, sect.ipt2)) / 2
				]
			)

			acft.beams[sect.ibeam].Mtosys[1, :] ⋅ θ
		end for sect in acft.wing_strips
	]

	plot_beams = ones(Bool, length(acft.beams))

	for (θ, strip) in zip(θs, acft.wing_strips)
		cv = acft.beams[strip.ibeam].Mtouni[:, 3] .* strip.c * sin(θ)
		cv[1] += strip.c

		x1, y1, z1 = points[:, strip.ipt1]
		x2, y2, z2 = points[:, strip.ipt2]

		ax.plot(
			[
				x1 + cv[1] * 3 / 4, 
				x2 + cv[1] * 3 / 4, 
				x2 - cv[1] / 4, 
				x1 - cv[1] / 4, 
				x1 + cv[1] * 3 / 4
			],
			[
				y1 + cv[2] * 3 / 4, 
				y2 + cv[2] * 3 / 4, 
				y2 - cv[2] / 4, 
				y1 - cv[2] / 4, 
				y1 + cv[2] * 3 / 4
			],
			[
				z1 + cv[3] * 3 / 4, 
				z2 + cv[3] * 3 / 4, 
				z2 - cv[3] / 4, 
				z1 - cv[3] / 4, 
				z1 + cv[3] * 3 / 4
			];
			color = "y"
		)

		plot_beams[strip.ibeam] = false

		ax.plot(
			[x1 + cv[1] * strip.e, x2 + cv[1] * strip.e],
			[y1 + cv[2] * strip.e, y2 + cv[2] * strip.e],
			[z1 + cv[3] * strip.e, z2 + cv[3] * strip.e];
			color = "k"
		)
		ax.plot(
			[x1 + cv[1] * strip.xCG, x2 + cv[1] * strip.xCG],
			[y1 + cv[2] * strip.xCG, y2 + cv[2] * strip.xCG],
			[z1 + cv[3] * strip.xCG, z2 + cv[3] * strip.xCG];
			color = "b"
		)
	end

	for (ib, b) in enumerate(acft.beams)
		if plot_beams[ib]
			x1, y1, z1 = points[:, b.ipt1]
			x2, y2, z2 = points[:, b.ipt2]

			ax.plot(
				[x1, x2],
				[y1, y2],
				[z1, z2];
				color = "k"
			)
		end
	end

	_plt.show()

end
