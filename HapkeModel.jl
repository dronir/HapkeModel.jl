
# Hapke's reflection model from: 
#  Hapke, Bruce: Bidirectional Reflectance Spectroscopy 5. The coherent backscatter oppsition effect
#  and anisotropic scattering (2002), Icarus 157, 523-534, doi:10.1006/icar.2002.6853
# References to equations are to this paper.

include("legendre.jl")

module HapkeModel

using LegendrePolynomial

export Isotropic, Rayleigh, HenyeyGreenstein, DoubleHenyeyGreenstein
export UnitP, SeriesP, PhaseFunction
export ScatteringModel, BDRF

# Hapke's notation, followed in the code:
# w = single-scattering albedo
# i = incidence angle
# e = emergence angle
# mu0 = cos(i)
# mu = cos(e)
# g = phase angle
# xi = Henyey-Greenstein asymmetry parameter

# Define some types for phase functions
abstract type PhaseFunction end

# Group phase functions into two types: those for which P(mu) and the P-constant
# are unity, and those for which it is expressed as a Legendre polynomial series
abstract type UnitP <: PhaseFunction end
abstract type SeriesP <: PhaseFunction end

struct Isotropic <: UnitP end
struct Rayleigh <: UnitP end

struct HenyeyGreenstein <: SeriesP
	xi::Float64
end

struct DoubleHenyeyGreenstein <: SeriesP
	c::Float64
	xi::Float64
end

struct ScatteringModel{P<:PhaseFunction}
	p::P
	w::Float64
	shadows::Bool
	hS::Float64
	P_const::Float64
end

ScatteringModel(p::PhaseFunction, w::Real) = ScatteringModel(p,w,false,0.0,Hapke_P_const(p))
ScatteringModel(p::PhaseFunction, w::Real, h::Real) = ScatteringModel(p,w,true,h,Hapke_P_const(p))
ScatteringModel(p::PhaseFunction, w::Real, E::Real, a::Real, phi::Real) = ScatteringModel(p,w,true,hS(E,a,phi),Hapke_P_const(p))


# Hapke's H-function approximation (eq. 13)
function H(x::Real, w::Real)
	gamma = sqrt(1.0 - w)
	r0 = (1-gamma) / (1.0 + gamma)
	return 1/(1 - w*x*(r0 + (1 - 2*r0*x)/2 * log((1+x)/x)))
end


# Hapke P function and 'P' constant are unity for the isotropic and Rayleigh phase functions
Hapke_P(p::UnitP, mu::Real) = 1.0
Hapke_P_const(p::UnitP) = 1.0

const MAX_ITER = 20

# A coefficients for P function series representations, eq. 26-27
A(n::Integer) = n%2==0 ? 0.0 : (-1)^((n+1)/2)/n * reduce(*, 1:2:n) / reduce(*, 2:2:n+1)

# b coefficients for P function series representations
b(n::Integer, p::HenyeyGreenstein) = (2n-1)*(-p.xi)^n
b(n::Integer, p::DoubleHenyeyGreenstein) = p.c * (2n+1) * p.xi^n

# Hapke P function and 'P' constant for the single and double HG phase functions
Hapke_P(p::SeriesP, mu::Real) = 1.0 + sum(i -> A(i) * b(i, p) * P(i, mu), 1:MAX_ITER)
Hapke_P_const(p::SeriesP) = 1.0 - sum(i -> A(i)^2 * b(i, p), 1:MAX_ITER)

# Hapke M function (eq. 17)
function M(model::ScatteringModel, mu0::Real, mu::Real) 
	h = H(mu, model.w) - 1
	h0 = H(mu0, model.w) - 1
	return h*Hapke_P(model.p, mu0) + h0*Hapke_P(model.p, mu) + model.P_const * h * h0
end


# Shadow hiding opposition effect term
B_S(g::Real, hS::Real) = 1/(1 + (1/hS)*tan(g/2)) # eq. 29
hS(E::Real, a::Real, phi::Real) = -E * a * log(1-phi) / (2phi) # eq. 30

# The different phase functions
phase(p::Isotropic, g::Real) = 1.0
phase(p::Rayleigh, g::Real) = 1.0 + 0.5 * P(2, cos(g))
phase(p::HenyeyGreenstein, g::Real) = (1 - p.xi^2) / (1 + p.xi^2 + 2*p.xi*cos(g))^1.5
function phase(p::DoubleHenyeyGreenstein, g::Real)
	a = (1 - p.xi^2)
	HG1 = a / (1 + p.xi^2 + 2*p.xi*cos(g))^1.5
	HG2 = a / (1 + p.xi^2 - 2*p.xi*cos(g))^1.5
	return (1+p.c)/2 * HG1 + (1-p.c)/2 * HG2
end


function BDRF(model::ScatteringModel, mu0::Real, mu::Real, g::Real)
	if model.shadows
		return model.w/(4pi) * mu0/(mu+mu0) * (phase(model.p, g)*B_S(g, model.hS) + M(model, mu0, mu))
	else
		return model.w/(4pi) * mu0/(mu+mu0) * (phase(model.p, g) + M(model, mu0, mu))
	end
end


end
