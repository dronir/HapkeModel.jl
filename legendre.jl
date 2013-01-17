module LegendrePolynomial

# This is a very quick and dirty implementation of Legendre polynomials

export P

# Associate Legendre polynomials
function P(n::Integer, m::Integer, t::Real)
	# Negative order
	if m < 0
		return (-1.0)^m * factorial(n-m) / factorial(n+m) * P(n,m,t)
	end
	
	if n==0
		return 1.0
	end
	
	if n==1 && m==0
		return t
	end

	if n==2 && m==0
		return 1.5 * t^2 - 0.5
	end

	u = -sqrt(1.0 - t^2)

	if n==1 && m==1
		return u
	end
	
	if n==2 && m==1
		return 3u*t
	end
	if n==2 && m==2
		return 3u^2
	end
	
	# Compute with recursion formulas
	if m == 0
		return ((2*(n-1) + 1) * t * P(n-1, 0, t) - (n-1)*P(n-2, 0, t)) / n
	end
	if n == m+1
		return t * (2m + 1) * P(m, m, t)
	end
	if n == m
		return (2n - 1) * u * P(n-1, n-1, t)
	end
	return P(n-2, m, t) + (2*n - 1) * u * P(n-1, m-1, t)
end

# Legendre polynomials
P(n::Integer, t::Real) = P(n, 0, t)

end # module