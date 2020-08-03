### A Pluto.jl notebook ###
# v0.11.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 9860665a-d56b-11ea-01c2-e1703fb29f5d
md"# Beammmmms"

# ╔═╡ c429ada6-d4e3-11ea-25ff-57d5831d9494
## Import
begin
	import Pkg
	Pkg.activate(".")
	
	using ForwardDiff
	using Plots
	using PlutoUI
end

# ╔═╡ 6ddca854-d4e3-11ea-0632-7900cf82dec6
pwd()

# ╔═╡ 5ff1e7b4-d587-11ea-2bad-a3948c7418c3
@bind N Slider(11:10:101)

# ╔═╡ c429d1aa-d588-11ea-1c32-b133861aced7
N

# ╔═╡ d39fc69c-d4e5-11ea-17de-bb83220b1625
@bind L html"<input type=range min = 0.5 max = 11 step = 0.1>"

# ╔═╡ 92f4bb3e-d4e3-11ea-0424-f59ce379b160
function basis1D(xh)

	# Define form function
	form1(x) = 1 - 3*x^2 + 2*x^3
	form2(x) = x*(x-1)^2
	form3(x) = 3*x^2 - 2*x^3
	form4(x) = x^2*(x-1)

	# Get step size
	h = xh[2]-xh[1]
	n = length(xh)


	phi = Array{Function}(undef,2*n,1)

	phi[1] = x -> form1(x./h).*(x<=h)
	phi[2] = x -> h*form2(x./h).*(x<=h)
	for i = 2:(n-1)
		phi[2*i-1] = x ->
					form3((x-xh[i-1])/h).*
					((x >= xh[i-1]) & (x < xh[i]))+
					form1((x-xh[i])/h).*
					((x >= xh[i]) & (x <= xh[i+1]))

		phi[2*i] = x ->
					h*form4((x-xh[i-1])/h).*
					((x >= xh[i-1]) & (x < xh[i]))+
					h*form2((x-xh[i])/h).*
					((x >= xh[i]) & (x <= xh[i+1]))
	end
	phi[end-1] = x -> form3((x - xh[end-1])/h).*((x >= xh[end-1]) & (x <= xh[end]))
	phi[end] = x -> h*form4((x - xh[end-1])/h).*((x >= xh[end-1]) & (x <= xh[end]))
	return phi
end


# ╔═╡ ecd0e724-d4e3-11ea-2d64-0d231a130b99
xh = LinRange(0,L,21)

# ╔═╡ f28245b4-d4e3-11ea-31e4-bdac63887998
phik = basis1D(xh);

# ╔═╡ ec7ee97c-d4e5-11ea-2528-af2a46f9004d
scatter(LinRange(0,L,N),phik[5].(LinRange(0,L,N)), label = 'ϕ')

# ╔═╡ 969bf7a4-d4e3-11ea-1a6e-25a799e92a6f
phidd = map(phik) do f
   firstderivative = x -> ForwardDiff.derivative(f, x)
   x -> ForwardDiff.derivative(firstderivative, x)
end

# ╔═╡ 0ae99762-d580-11ea-337b-b1cb1e556277
Pkg.instantiate()

# ╔═╡ Cell order:
# ╠═9860665a-d56b-11ea-01c2-e1703fb29f5d
# ╠═6ddca854-d4e3-11ea-0632-7900cf82dec6
# ╠═5ff1e7b4-d587-11ea-2bad-a3948c7418c3
# ╠═c429d1aa-d588-11ea-1c32-b133861aced7
# ╠═d39fc69c-d4e5-11ea-17de-bb83220b1625
# ╠═ec7ee97c-d4e5-11ea-2528-af2a46f9004d
# ╠═92f4bb3e-d4e3-11ea-0424-f59ce379b160
# ╠═ecd0e724-d4e3-11ea-2d64-0d231a130b99
# ╠═f28245b4-d4e3-11ea-31e4-bdac63887998
# ╠═969bf7a4-d4e3-11ea-1a6e-25a799e92a6f
# ╠═c429ada6-d4e3-11ea-25ff-57d5831d9494
# ╠═0ae99762-d580-11ea-337b-b1cb1e556277
