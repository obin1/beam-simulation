### A Pluto.jl notebook ###
# v0.11.14

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

# ╔═╡ e125c7d6-ef97-11ea-39fd-b176103f48ba
begin 
	using ForwardDiff
	using Plots
	using PlutoUI
	using SparseArrays
	pyplot()
end

# ╔═╡ 219fcfc0-ef96-11ea-3baa-85eb1e81cf63
md"# Euler Bernoulli Beams"

# ╔═╡ 7171003c-f65f-11ea-0ed1-eb09bbaa9a32
begin
	
	md"""**Let us define the properties and loading characterstics of our beam**
	
	Analysis: $(@bind analysis Select(["static" => "Static" ,"dynamic" => "Dynamic"]))
	
	Type:  $(@bind t Select(["c" => "Cantilever" ,"ss" => "Simply Supported"]))
	
	Length of the beam (L): 1 $(@bind L Slider(1:1:10;default=1)) 10

	Number of nodes (N): 2 $(@bind N Slider(2:1:30;default=6)) 30
	
	Bending modulus (EI) = $(@bind ei NumberField(1:0.01:5;default=1))
	
	μ Mass/length = $(@bind μ NumberField(1:0.01:5;default=1))
	
	Uniformly distributed Load (q) = $(@bind ul NumberField(-10:0.1:10;default=0))"""
	
end

# ╔═╡ 4c1d7992-f668-11ea-01e2-6300503d267d
begin
if analysis == "dynamic"
	md"""**Dynamic pameters**
		
	Total time steps (Mtau) = $(@bind Mtau NumberField(1:1:500;default=200))
	
	stepping(tau) = $(@bind tau NumberField(0:0.001:1;default=0.1))
	
	β = $(@bind beta NumberField(0:0.001:1;default=0.25))
	
	γ = $(@bind gamma NumberField(0:0.001:1;default=0.5))"""
end
end

# ╔═╡ 65d5625e-f58b-11ea-1cb7-fb37ec26e2c3
begin
	if t=="c"
		md"""**Boundary Conditions for a cantilever beam**
	
			a = $(@bind a Slider(0:10)) 
	
			b = $(@bind b NumberField(-1:1;default=0))  
	
			End Load (Q) = $(@bind Q NumberField(-10:0.01:10;default=1)) 
	
			End Moment (M) = $(@bind ML NumberField(-1:0.01:1;default=0))"""
	elseif t=="ss"
		md"""**Boundary Conditions for a Simply Supported beam**
			
			a0 = $(@bind a html"<input type=range min=0 max=10 value=0>") 
	
			aL = $(@bind b html"<input type=number value=0>")  
	
			Moment at x=0 (M0) = $(@bind Q NumberField(-10:0.01:10;default=1))
	
			Moment at x=L (ML) = $(@bind ML NumberField(-10:0.01:10;default=1))"""
		end
end

# ╔═╡ f2b09b22-f03a-11ea-2baf-ff3e1721752d
md" And so the governing PDE for the above set of conditions will be <EDIT HERE>
.
We intend to solve this system with the Finite Element Method and analyze the results."

# ╔═╡ 00fbb85c-f508-11ea-21b0-1d20d4bdb05d
md"#### Generate domain"

# ╔═╡ 16101166-f4de-11ea-25f8-bdfa0e83d2a4
function gen_elements(L,N)
	xh=LinRange(0,L,N)
	ne=N-1
	conn=[[e,e+1] for e=1:ne]
	return xh,ne,conn
end

# ╔═╡ 548dde12-f508-11ea-2ef8-fffbcd320ec8
xh,ne,conn=gen_elements(L,N);

# ╔═╡ 2b5fe7a2-f04a-11ea-3f33-07ef382f7722
scatter(xh,10*ones(N),yaxis=nothing,ylims=(0,20),label="\$\\Omega\$",size=(700,200))

# ╔═╡ 2b84ae9c-f66d-11ea-38e4-c1e6e9e7b7b8
#Assumption of constant q,EI,mu throughout domain. Change for different distribution
function constant_dist(ul,μ,ei)
	qv = ul*ones(ne);
	EI = ei*ones(ne);
	mu = μ*ones(ne);
	return qv,mu,EI
end

# ╔═╡ 8cbd99bc-f66d-11ea-3c74-2b9b9999741d
qv,mu,EI=constant_dist(ul,μ,ei); #taking constant distribution of q, EI and μ on domain

# ╔═╡ 02a16210-f503-11ea-2552-7b72536ee06d
md"#### Generating our Ansatz space"

# ╔═╡ 466ce928-f03b-11ea-1b7e-a9e26cbccdcc
function basis1D(xh,N)
	
	# Define form function
	form(y) = [1 - 3*y^2 + 2*y^3, y*(y-1)^2,3*y^2 - 2*y^3, y^2*(y-1)]
	formdd(y) = [-6+12*y, -4+6*y, 6-12*y,-2+6*y] 
	
	h = xh[2]-xh[1] #uniform grid spacing assumed

	#basis function definition
	phi = Array{Function}(undef,2*N,1)
	phidd = Array{Function}(undef,2*N,1)
	
	phi[1] = x -> form(x./h)[1].*(x<=h)
	phi[2] = x -> h*form(x./h)[2].*(x<=h)
	for i = 2:(N-1)
		phi[2*i-1] = x ->
					form((x-xh[i-1])/h)[3].*
					((x >= xh[i-1]) & (x < xh[i]))+
					form((x-xh[i])/h)[1].*
					((x >= xh[i]) & (x <= xh[i+1]))

		phi[2*i] = x ->
					h*form((x-xh[i-1])/h)[4].*
					((x >= xh[i-1]) & (x < xh[i]))+
					h*form((x-xh[i])/h)[2].*
					((x >= xh[i]) & (x <= xh[i+1]))
	end
	phi[end-1] = x -> form((x - xh[end-1])/h)[3].*((x >= xh[end-1]) & (x <= xh[end]))
	phi[end] = x -> h*form((x - xh[end-1])/h)[4].*((x >= xh[end-1]) & (x <= xh[end]))
	
	#derivatives of basis functions
	#1
	phid = map(phi) do f
		x -> ForwardDiff.derivative(f, x);
	end
	
	#2
	phidd[1] = x -> formdd(x./h)[1].*(x<=h)
	phidd[2] = x -> h*formdd(x./h)[2].*(x<=h)
	for i = 2:(N-1)
		phidd[2*i-1] = x ->
					formdd((x-xh[i-1])/h)[3].*
					((x >= xh[i-1]) & (x < xh[i]))+
					formdd((x-xh[i])/h)[1].*
					((x >= xh[i]) & (x <= xh[i+1]))

		phidd[2*i] = x ->
					h*formdd((x-xh[i-1])/h)[4].*
					((x >= xh[i-1]) & (x < xh[i]))+
					h*formdd((x-xh[i])/h)[2].*
					((x >= xh[i]) & (x <= xh[i+1]))
	end
	phidd[end-1] = x -> formdd((x - xh[end-1])/h)[3].*((x >= xh[end-1]) & (x <= xh[end]))
	phidd[end] = x -> h*formdd((x - xh[end-1])/h)[4].*((x >= xh[end-1]) & (x <= xh[end]))
	
	return phi,phid,phidd
end

# ╔═╡ 2c02f7ba-f4ed-11ea-19c9-9bb69c5de35b
function local2global(ne)
	return [[2*e-1,2*e,2*e+1,2*e+2] for e=1:ne]
end

# ╔═╡ 3f7a5a08-f03e-11ea-1fcf-83e5b5e57d28
phi,phid,phidd= basis1D(xh,N);

# ╔═╡ 998838b6-f508-11ea-0949-8b40b202bab8
dofs=local2global(ne);

# ╔═╡ 098b7a0a-f2c1-11ea-31d6-91b4b485ac41
md"###### Let us check whether our basis functions can approximate a function"

# ╔═╡ dbb675c8-f2b4-11ea-3ffa-81e83739270d
function check(xh,f)
	n=length(xh)
	fd = xh -> ForwardDiff.derivative(f,xh)
	u = zeros(2*n)
	for i = 1:N
		u[2*i-1] = f.(xh[i])
		u[2*i] = fd.(xh[i])
	end
	p = sum((u[2*i-1]*phi[2*i-1].(xh) + u[2*i]*phi[2*i].(xh)) for i=1:n)
	plot(xh,f.(xh),label="actual")
	scatter!(xh,p,label="piecewise cubic")
end

# ╔═╡ c78ce07a-f2bc-11ea-28bd-a1df0c486ed1
begin
	f(x)=cos(x)*sin(x)
	check(xh,f::Function);
end

# ╔═╡ 04fb4452-f2c1-11ea-393c-4b9b88d2aaf8
md"### Calculation of element FE matrices and vectors"

# ╔═╡ 9aea9122-f5d0-11ea-1bc7-d5b5f6368033
function element_mass(mu_e,x_e,phi_e)
	dofs = length(phi_e)
	M_e = zeros(dofs,dofs)
	gp = [-sqrt((3/7) - (2/7)*sqrt(6/5)),sqrt((3/7) - (2/7)*sqrt(6/5)),-sqrt((3/7) + (2/7)*sqrt(6/5)),sqrt((3/7) + (2/7)*sqrt(6/5))]
	gw = [(18+sqrt(30))/36,(18+sqrt(30))/36,(18-sqrt(30))/36,(18-sqrt(30))/36]
	h=x_e[2]-x_e[1]
	
	for i=1:length(gp)
		p=(0.5*gp[i]+0.5)*h + x_e[1]
		for j=1:dofs
			for k=1:dofs	
				M_e[j,k] += (0.5*mu_e*gw[i]*phi_e[j](p)*phi_e[k](p))
			end
		end
	end
	return M_e
end

# ╔═╡ 6c4aa6b8-f39b-11ea-0d7e-6bf87aaa95be
function element_stiffness(EI_e,x_e,phidd_e)
	dofs = length(phidd_e)
	K_e = zeros(dofs,dofs)
	gp = [-1/sqrt(3),1/sqrt(3)]
	gw = [1,1]
	h=x_e[2]-x_e[1]
	for i=1:length(gp)
		p=(0.5*gp[i]+0.5)*h + x_e[1]
		for j=1:dofs
			for k=1:dofs	
				K_e[j,k] += (0.5*(EI_e/h^3)*gw[i]*phidd_e[j](p)*phidd_e[k](p))
			end
		end
	end
	return K_e
end

# ╔═╡ 8352e10c-f4ce-11ea-25a9-d75c910445fa
function element_forces(q_e,x_e,phi_e)
	dofs =length(phi_e)
	f_e = zeros(dofs)
	gp = [-sqrt(3/5),0,sqrt(3/5)]
	gw = [5/9,8/9,5/9]
	h=x_e[2]-x_e[1]
	for i=1:length(gp)
		p=(0.5*gp[i]+0.5)*h + x_e[1]
		for j=1:dofs
			f_e[j]+= (0.5*q_e*h*gw[i]*phi_e[j](p))
		end
	end
	return f_e
end

# ╔═╡ a94fb776-f66e-11ea-205e-21f70a68ef1d
md"### Global assembly"

# ╔═╡ 0c893646-f66f-11ea-1ca9-e32285b8c42b
md"##### Static case"

# ╔═╡ 680184dc-f422-11ea-3aa4-e7d07918b79b
function assemble(t,EI,Q,ML,qv,a,b,xh,ne,conn,phi,phidd,dofs)
	ndofs=length(phi);
	nphi= length(dofs[1]); 
	q = zeros(ndofs)
	ii = zeros(Int64, ne, nphi, nphi); # sparse i-index
	jj = zeros(Int64, ne, nphi, nphi); # sparse j-index
	aa = zeros(ne, nphi, nphi); # entry of Galerkin matrix
	
	for e=1:ne
		sloc= element_stiffness(EI[e],xh[conn[e]],phidd[dofs[e]])
		floc= element_forces(qv[e],xh[conn[e]],phi[dofs[e]])
		q[dofs[e]]+=floc[:]
		for j=1:nphi
			for k=1:nphi
				ii[e,j,k] = dofs[e][j]; # local-to-global
            	jj[e,j,k] = dofs[e][k]; # local-to-global
				aa[e,j,k] = sloc[j,k];
			end
		end
	end
	S = sparse(ii[:],jj[:],aa[:]);		
	
	#BCs
	e0 = [phi[i](0) for i=1:ndofs];
	eL = [phi[i](L) for i=1:ndofs];
	d0 = [phid[i](0) for i=1:ndofs];
	dL = [phid[i](L) for i=1:ndofs];
	 
	if t=="c"
		C = [e0 d0];
		q_e = [q + Q*eL + ML*dL; a; b];
	elseif t=="ss"
		C = [e0 -eL];
		q_e = [q - Q*d0 + ML*dL; a; -b];
	end
		
	#Extended matrix system with BCs
	S_e = [S C; C' zeros(2,2)];
	
	return S_e,q_e;
end	

# ╔═╡ 84c0f1b8-f5d0-11ea-14e6-c987ddf8da67
md"""##### Dynamic case"""

# ╔═╡ 986616bc-f5dd-11ea-3975-bb8d25062cce
function dyn_assemble(mu,EI,Q,ML,qv,xh,ne,conn,phi,phidd,dofs)
	ndofs=length(phi);
	nphi= length(dofs[1]); 
	q = zeros(ndofs)
	ii = zeros(Int64, ne, nphi, nphi); # sparse i-index
	jj = zeros(Int64, ne, nphi, nphi); # sparse j-index
	aa = zeros(ne, nphi, nphi); # entry of Galerkin matrix
	bb = zeros(ne, nphi, nphi); #entry of mass matrix
	for e=1:ne
		sloc= element_stiffness(EI[e],xh[conn[e]],phidd[dofs[e]])
		mloc= element_mass(mu[e],xh[conn[e]],phi[dofs[e]])
		floc= element_forces(qv[e],xh[conn[e]],phi[dofs[e]])
		q[dofs[e]]= q[dofs[e]]+floc[:]
		for j=1:nphi
			for k=1:nphi
				ii[e,j,k] = dofs[e][j]; # local-to-global
            	jj[e,j,k] = dofs[e][k]; # local-to-global
				aa[e,j,k] = sloc[j,k];
				bb[e,j,k] = mloc[j,k];
			end
		end
	end
	S = sparse(ii[:],jj[:],aa[:]);	
	M = sparse(ii[:],jj[:],bb[:]);
	
	#BCs
	e0 = [phi[i](0) for i=1:ndofs];
	eL = [phi[i](L) for i=1:ndofs];
	d0 = [phid[i](0) for i=1:ndofs];
	dL = [phid[i](L) for i=1:ndofs];
	 
	if t=="c"
		C = [e0 d0];
		q_e = [q + Q*eL + ML*dL; a; b];
	elseif t=="ss"
		C = [e0 -eL];
		q_e = [q - Q*d0 + ML*dL; a; -b];
	end
		
	#Extended matrix system with BCs
	S_e = [S C; C' zeros(2,2)];
	M_e = [S 0*C; 0*C' zeros(2,2)];
	
	return S_e,M_e,q_e;
end	

# ╔═╡ b59aacb2-f50d-11ea-34dd-67c479511398
md"### Solving the system"

# ╔═╡ b9cb119e-f5e8-11ea-302d-6db9fec26b2b
function newmark_step(Mtau,xh,w,tau,beta,gamma,M_e,S_e,phi)
	w_interp = zeros(length(xh),Mtau);
	udot = 0*w;
	udot[end-1:end] =w[end-1:end];
	udotdot =  0*w;
	udotdot[end-1:end] =w[end-1:end];
	for i = 1:Mtau
		ustar = w + udot*tau + (1/2 - beta).*udotdot.*tau^2;
		ustardot= udot + (1-gamma)*tau*udotdot;
		b= -S_e*ustar;
		udotdot = (M_e + (beta*tau^2)*S_e)\b;
		w = ustar + beta*tau^2*udotdot;
		udot = ustardot + gamma*tau*udotdot;
		w_interp[:,i] = sum((w[2*j-1]*phi[2*j-1].(xh) + w[2*j]*phi[2*j].(xh)) for j=1:length(xh));
	end
	return w_interp
end

# ╔═╡ 28da6d86-f4f7-11ea-010d-f7838116054e
begin
	if analysis=="static"
		S_e,q_e = assemble(t,EI,Q,ML,qv,a,b,xh,ne,conn,phi,phidd,dofs);
		w = S_e\q_e;
		w_interp = sum((w[2*j-1]*phi[2*j-1].(xh) + w[2*j]*phi[2*j].(xh)) for j=1:length(xh));
		
	elseif analysis=="dynamic"
		S_e,M_e,q_e = dyn_assemble(mu,EI,Q,ML,qv,xh,ne,conn,phi,phidd,dofs);
		w = S_e\q_e;
		w_interp = newmark_step(Mtau,xh,w,tau,beta,gamma,M_e,S_e,phi);
		
	end
	nothing
end

# ╔═╡ e2536ecc-f674-11ea-386d-274c40e5a846
md"### Results and Post Processing"

# ╔═╡ 1b42f33a-f510-11ea-128a-c5dadf723476
function check_analytical(t,a,b,Q,ML,q,L,ei,xh)
	w_exact=zeros(length(xh))
	if iszero(b)	
		if t=="c"
			if iszero(ML) && iszero(q)
				w_exact = (Q.*xh.*xh./(6*ei)).*(3*L.-xh)
				text="with end load"
			elseif iszero(ML) && iszero(Q)
				w_exact = (q*xh.*xh./(24*ei)).*(6*L^2 .- 4*L*xh .+ xh.^2)
				text="with distributed load"
			elseif iszero(q) && iszero(Q)
				w_exact = (ML*xh.^2)/(2*ei)
				text="with end moment"
			else 
				text="Analytical solution not coded in!"
			end
		elseif t=="ss"
			if iszero(ML) && iszero(Q)
				w_exact = (q*xh./(24*ei)).*(L^3-2*L*xh.*xh + xh.^3)
				text="with distibuted load"
			elseif iszero(q)
				w_exact = (ML*xh).*(L.-xh)/(2*ei)
				text="with end moments"
			else
				text="Analytical solution not coded in!"
			end
		end
	else
		text="Analytical solution not coded in!"
	end
	
return w_exact,text;
end

# ╔═╡ b2b57896-f696-11ea-07b0-2f01a202fbd0
function plotFEM(xh,wh,ylim=[-L/2,L/2])
	#assumption of linear interpolation in between nodes while plotting. 
	plot(xh,wh,
		xlims=(0,L+0.05),
		ylims=(ylim[1],ylim[2]),
		label="FEM",
		xlabel="Domain (Ω)",
		ylabel="Deflection (w)",
		marker=(:circle,:orange),title="Static Case")
end	

# ╔═╡ a32eccb4-f6ab-11ea-2b2c-1de81ae66813
md"#### Static"

# ╔═╡ f64ad968-f68f-11ea-2a36-7b7bdac2b87d
if analysis=="static"
		md"""**Check Analytical solution** $(@bind check_sol CheckBox())
			 **Plot Shear Force distribution** $(@bind shear_stress CheckBox())
			 **Plot Bending Moment distribution** $(@bind bending_moment CheckBox())
	"""
end

# ╔═╡ 68166060-f5b1-11ea-0e76-41f2220403c6
begin
	if analysis=="static"
	y_lim=[minimum(w_interp) maximum(w_interp)] 
	y_lim = y_lim+0.5*y_lim
		if check_sol==false
			plotFEM(xh,w_interp,y_lim)
		else 
			w_ex, ctext = check_analytical(t,a,b,Q,ML,ul,L,ei,xh)
			plot!(xh,w_ex,
				annotations=(L/3,y_lim[2]-0.5*y_lim[2],ctext),
				label="Analytical",
				line=(:dash)
				)
		end
		else
			"Select the static case to see these plots"
	end
end

# ╔═╡ 96f1baee-f6ab-11ea-2dd3-57f8e2e541d8
md"#### Dynamic"

# ╔═╡ a6cb3bd2-f6a1-11ea-09e0-198c23e9b2a8
if analysis=="dynamic"
	anim = @animate for T ∈ 1:Mtau
		plot(xh,w_interp[:,T],
			xlims=(0,L+0.05),
			ylims=(-1,1),
			label="FEM solution",
			xlabel="Domain (Ω)",
			ylabel="Deflection (w)",
			marker=(:circle,:orange),
			annotations=(L/3,L/1.15,"T=$(T)"))
		plot!(title="Dynamic case")
	end
	gif(anim)
else 
	"select the dynamic case to see this clip"
end

# ╔═╡ 16792298-f6a5-11ea-07a8-39d00521b9db
if analysis=="dynamic"
	md"#### Plotting first five modes"
end

# ╔═╡ f0bebeb0-f6a8-11ea-2eaf-a72972cd414e


# ╔═╡ Cell order:
# ╟─219fcfc0-ef96-11ea-3baa-85eb1e81cf63
# ╠═e125c7d6-ef97-11ea-39fd-b176103f48ba
# ╟─7171003c-f65f-11ea-0ed1-eb09bbaa9a32
# ╟─4c1d7992-f668-11ea-01e2-6300503d267d
# ╟─65d5625e-f58b-11ea-1cb7-fb37ec26e2c3
# ╟─f2b09b22-f03a-11ea-2baf-ff3e1721752d
# ╟─00fbb85c-f508-11ea-21b0-1d20d4bdb05d
# ╟─16101166-f4de-11ea-25f8-bdfa0e83d2a4
# ╠═548dde12-f508-11ea-2ef8-fffbcd320ec8
# ╠═2b5fe7a2-f04a-11ea-3f33-07ef382f7722
# ╟─2b84ae9c-f66d-11ea-38e4-c1e6e9e7b7b8
# ╠═8cbd99bc-f66d-11ea-3c74-2b9b9999741d
# ╟─02a16210-f503-11ea-2552-7b72536ee06d
# ╟─466ce928-f03b-11ea-1b7e-a9e26cbccdcc
# ╟─2c02f7ba-f4ed-11ea-19c9-9bb69c5de35b
# ╠═3f7a5a08-f03e-11ea-1fcf-83e5b5e57d28
# ╠═998838b6-f508-11ea-0949-8b40b202bab8
# ╟─098b7a0a-f2c1-11ea-31d6-91b4b485ac41
# ╟─dbb675c8-f2b4-11ea-3ffa-81e83739270d
# ╠═c78ce07a-f2bc-11ea-28bd-a1df0c486ed1
# ╟─04fb4452-f2c1-11ea-393c-4b9b88d2aaf8
# ╟─9aea9122-f5d0-11ea-1bc7-d5b5f6368033
# ╟─6c4aa6b8-f39b-11ea-0d7e-6bf87aaa95be
# ╟─8352e10c-f4ce-11ea-25a9-d75c910445fa
# ╟─a94fb776-f66e-11ea-205e-21f70a68ef1d
# ╟─0c893646-f66f-11ea-1ca9-e32285b8c42b
# ╟─680184dc-f422-11ea-3aa4-e7d07918b79b
# ╟─84c0f1b8-f5d0-11ea-14e6-c987ddf8da67
# ╟─986616bc-f5dd-11ea-3975-bb8d25062cce
# ╟─b59aacb2-f50d-11ea-34dd-67c479511398
# ╟─b9cb119e-f5e8-11ea-302d-6db9fec26b2b
# ╟─28da6d86-f4f7-11ea-010d-f7838116054e
# ╟─e2536ecc-f674-11ea-386d-274c40e5a846
# ╟─1b42f33a-f510-11ea-128a-c5dadf723476
# ╟─b2b57896-f696-11ea-07b0-2f01a202fbd0
# ╟─a32eccb4-f6ab-11ea-2b2c-1de81ae66813
# ╟─f64ad968-f68f-11ea-2a36-7b7bdac2b87d
# ╟─68166060-f5b1-11ea-0e76-41f2220403c6
# ╟─96f1baee-f6ab-11ea-2dd3-57f8e2e541d8
# ╟─a6cb3bd2-f6a1-11ea-09e0-198c23e9b2a8
# ╟─16792298-f6a5-11ea-07a8-39d00521b9db
# ╠═f0bebeb0-f6a8-11ea-2eaf-a72972cd414e
