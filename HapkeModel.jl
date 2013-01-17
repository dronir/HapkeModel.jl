
# Hapke's reflection model from: 
#  Hapke, Bruce: Bidirectional Reflectance Spectroscopy 5. The coherent backscatter oppsition effect
#  and anisotropic scattering (2002), Icarus 157, 523-534, doi:10.1006/icar.2002.6853
# References to equations are to this paper.

module HapkeModel

require("legendre.jl")
using LegendrePolynomial

export BDRF

# Hapke's notation, followed in the code:
# w = single-scattering albedo
# i = incidence angle
# e = emergence angle
# mu0 = cos(i)
# mu = cos(e)
# g = phase angle
# xi = Henyey-Greenstein asymmetry parameter

# Define some types for phase functions
abstract PhaseFunction

type Isotropic <: PhaseFunction end
type Rayleigh <: PhaseFunction end

type HenyeyGreenstein <: PhaseFunction 
	xi::Float64
end

type DoubleHenyeyGreenstein <: PhaseFunction 
	c::Float64
	xi::Float64
end


# Hapke's H-function approximation (eq. 13)
function H(x::Real, w::Real)
	gamma = sqrt(1-w)
	r0 = (1-gamma) / (1+gamma)
	return 1/(1 - w*x*(r0 + (1 - 2*r0*x)/2 * log((1+x)/x)))
end

# Coefficient for hapke P function series expansion
global const MAX_ITER = 20
A(n::Integer) = n%2==0 ? 0.0 : (-1)^((n+1)/2)/n * reduce(*, 1:2:n) / reduce(*, 2:2:n+1) # Eq. 26-27
b(n::Integer, xi::Real) = (2n-1)*(-xi)^n # For HG function

# Series definitions (eq. 23-25)
Hapke_P(xi::Real, mu::Real) = 1.0 + sum(i->A(i)*b(i, xi)*P(i, mu), 1:MAX_ITER)
Hapke_P_const(xi::Real) = 1.0 - sum(i->A(i)^2*b(i,xi), 1:MAX_ITER)

# Hapke M function for HG phase function (eq. 17)
function M(xi::Real, mu::Real, mu0::Real) 
	h = H(mu) - 1
	h0 = H(mu0) - 1
	return h*Hapke_P(xi, mu0) + h0*Hapke_P(xi, mu) + Hapke_P_const(xi) * h * h0
end


# Shadow hiding opposition effect term
B_S(g::Real, hS::Real) = 1/(1 + (1/hS)*tan(g/2)) # eq. 29
hS(E::Real, a::Real, phi::Real) = -E * a * ln(1-phi) / (2phi) # eq. 30

phase(p::Isotropic, g::Real) = 1.0
phase(p::Rayleigh, g::Real) = 1.0 + 0.5 * P(2, cos(g))
phase(p::HenyeyGreenstein, g::Real) = (1 - p.xi^2) / (1 + p.xi^2 + 2*p.xi*cos(g))^1.5
function phase(p::DoubleHenyeyGreenstein, g::Real)
	a = (1 - p.xi^2)
	HG1 = a / (1 + p.xi^2 + 2*p.xi*cos(g))^1.5
	HG2 = a / (1 + p.xi^2 - 2*p.xi*cos(g))^1.5
	return (1+c)/2 * HG1 + (1-c)/2 * HG2
end

function BDRF(mu::Real, mu0::Real, g::Real, w::Real, p::PhaseFunction)
	return w/(4pi) * mu0/(mu+mu0) * (p(g) + M(mu, mu0))
end


end