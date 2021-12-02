using Base: Float64
"""
File with definitions of functions for structural analysis of aircraft configurations using beam elements
"""

Fmax = 1e15

"""
Get elasticity matrix for a single beam
"""
beam_get_K(
	L::Fg,
	EA::Fg,
	GJ::Fg,
	EIy::Fg,
	EIz::Fg
) where {Fg <: Real} = - Fg[
	(EA / L) 0.0 0.0 0.0 0.0 0.0 (- EA / L) 0.0 0.0 0.0 0.0 0.0;
	0.0 (12.0 * EIz / L ^ 3) 0.0 0.0 0.0 (6.0 * EIz / L ^ 2) 0.0 (- 12.0 * EIz / L ^ 3) 0.0 0.0 0.0 (6.0 * EIz / L ^ 2);
	0.0 0.0 (12.0 * EIy / L ^ 3) 0.0 (- 6.0 * EIy / L ^ 2) 0.0 0.0 0.0 (- 12.0 * EIy / L ^ 3) 0.0 (- 6.0 * EIy / L ^ 2) 0.0;
	0.0 0.0 0.0 (GJ / L) 0.0 0.0 0.0 0.0 0.0 (- GJ / L) 0.0 0.0;
	0.0 0.0 (- 6.0 * EIy / L ^ 2) 0.0 (4.0 * EIy / L) 0.0 0.0 0.0 (6.0 * EIy / L ^ 2) 0.0 (2.0 * EIy / L) 0.0;
	0.0 (6.0 * EIz / L ^ 2) 0.0 0.0 0.0 (4.0 * EIy / L) 0.0 (- 6.0 * EIz / L ^ 2) 0.0 0.0 0.0 (2.0 * EIz / L);
	(- EA / L) 0.0 0.0 0.0 0.0 0.0 (EA / L) 0.0 0.0 0.0 0.0 0.0;
	0.0 (- 12.0 * EIz / L ^ 3) 0.0 0.0 0.0 (- 6.0 * EIz / L ^ 2) 0.0 (12.0 * EIz / L ^ 3) 0.0 0.0 0.0 (- 6.0 * EIz / L ^ 2);
	0.0 0.0 (- 12.0 * EIy / L ^ 3) 0.0 (6.0 * EIy / L ^ 2) 0.0 0.0 0.0 (12.0 * EIy / L ^ 3) 0.0 (6.0 * EIy / L ^ 2) 0.0;
	0.0 0.0 0.0 ( - GJ / L) 0.0 0.0 0.0 0.0 0.0 (GJ / L) 0.0 0.0;
	0.0 0.0 (- 6.0 * EIy / L ^ 2) 0.0 (2.0 * EIy / L) 0.0 0.0 0.0 (6.0 * EIy / L ^ 2) 0.0 (4.0 * EIy / L) 0.0;
	0.0 (6.0 * EIz / L ^ 2) 0.0 0.0 0.0 (2.0 * EIz / L) 0.0 (- 6.0 * EIz / L ^ 2) 0.0 0.0 0.0 (4.0 * EIz / L)
]
