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
	using PlutoUI
	using SparseArrays
	using Arpack
	using ForwardDiff
	using Plots	
	pyplot()
end

# ╔═╡ 219fcfc0-ef96-11ea-3baa-85eb1e81cf63
md"# Euler Bernoulli Beams"

# ╔═╡ c7eb8ae0-fffe-11ea-3f09-318f505b326e
md"""
**Numerical Analysis Project with Prof. Michael Karow** \
**TU Berlin** \
Summer semester 2020 \
**Group members**: Arvind Nayak, Obin Sturm, Fons van der Plans, Karolina Siemieniuk
"""

# ╔═╡ 6646b4f0-fdbe-11ea-38a6-9ff1504f8f5d
html"""
<style>

body.hide_all_inputs pluto-cell {
	min-height: 0.5px;
	margin-top: 6px;
}
</style>
<script>
document.body.classList.toggle("hide_all_inputs")
</script>
"""

# ╔═╡ 2ff6cc50-0002-11eb-1d90-e1f3bdada8fc
md"""
**Summary**\

Add summary
"""

# ╔═╡ d8e7f8c0-00d4-11eb-3aad-d5b4512646fa
md"""
#### **1. Introduction**
Structures are all around us. Constantly deflecting and deforming. To measure the performance of a structure deflections and deformations need to be investigated. There are many reason why those performance measures are important. One of them is safety, but the comfort of the users is a vital consideration as well. For example a floor needs to be stiff enough so that it does not deflect too much and users do not feel uncomfortable, even though the deflection might be safe. Another example would be buildings deflecting due to wind or bridges deflecting under the load of the cars going through it. There are various methods to calculate those parameters and this project is going to focus on the Euler Bernoulli beam theory. It is a model that was developed in the 18th century and is a method still used today to analyse the behaviour of bending elements and to determine how beams behave under axial forces and beding.



"""

# ╔═╡ 3840efd0-fffd-11ea-0445-4184c1500cf5
md"""

In this project the bending equation for an elastic Bernoulli beam is solved via the Finite Element Method (FEM). The aim of the project is to solve both, static and dynamic cases of the beam problem and generate:\
- pictures of bended beam in equilibrium (for the static case) and \
- video showing vibrating beam (for the dynamic case).\
The FEM was completed using piecewise differentiable functions instead of Sobolev Spaces that are most commonly used when working with FEM.
The project has been completed using Pluto editor for Julia."""

# ╔═╡ 7171003c-f65f-11ea-0ed1-eb09bbaa9a32
begin
	
	md"""
	#### 2. Parameters - user's choice
	
	**In this section, the properties and loading characterstics of the beam can be selected in order to find the solution of the case that the user is interested in. This will accordingly update the code and the Notebook will calculate and display the solution for the chosen case.**
	
	Analysis: $(@bind analysis Select(["static" => "Static" ,"dynamic" => "Dynamic"]))  Type:  $(@bind t Select(["c" => "Cantilever" ,"ss" => "Simply Supported"]))
	
	Length of the beam (L): 1 $(@bind L Slider(1:1:10;default=1)) 10

	Number of nodes (N): 2 $(@bind N Slider(2:1:100;default=21)) 100
	
	Bending modulus (EI) = $(@bind ei NumberField(1:0.01:10;default=1))
	
	μ Mass/length = $(@bind μ NumberField(1:0.01:5;default=1))"""
	
end

# ╔═╡ 65d5625e-f58b-11ea-1cb7-fb37ec26e2c3
begin
	if t=="c"
		md"""**Boundary Conditions for a cantilever beam**
	
			a = $(@bind a NumberField(0:2;default=0)) 
	
			b = $(@bind b NumberField(-1:1;default=0))  
	
			End Load (QL) = $(@bind Q NumberField(-10:0.01:10;default=1)) 
	
			End Moment (ML) = $(@bind ML NumberField(-1:0.01:1;default=0))"""
	elseif t=="ss"
		md"""**Boundary Conditions for a Simply Supported beam**
			
			a0 = $(@bind a NumberField(-10:0.01:10;default=0)) 
	
			aL = $(@bind b NumberField(-10:0.01:10;default=0))  
	
			Moment at x=0 (M0) = $(@bind Q NumberField(-10:0.01:10;default=-1))
	
			Moment at x=L (ML) = $(@bind ML NumberField(-10:0.01:10;default=1))"""
		end
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

# ╔═╡ 58ce39ee-fe93-11ea-3244-f5d7f3935485
md"""UDL (q) = $(@bind ul NumberField(-10:0.1:10;default=0))"""

# ╔═╡ f2b09b22-f03a-11ea-2baf-ff3e1721752d
md" Using the parameters above, we intend to solve the system with the `Finite Element Method`. "

# ╔═╡ 66e54c52-0002-11eb-27d1-45419b83cfa6
md"""
#### **3. Static bending equation**

The static bending equation is described below

$(EIw'')''(x) = q(x), x \in [0,L], q \in V,$ \
where the set of piecewise twice differentiable functions is defined as\
$V$ := $C^{2,p}(0,L) = \{ \phi : [0,L] \rightarrow \mathbb{R}| \phi \in C(0,L), \phi' \in C(0,L), \phi'' \in C^p(0,L) \}$ \
and\

$w(x)$ - height of the neutral axis (bending curve) at $x$ \
$E = E(x)$ Young's module \
$I = I(x)$ area moment of inertia: $I(x) = \int z^2 dydz$ \
$q(x)$ load (force density) at $x$\
Further notation:\

$M^x(w)$ = $EIw''(x)$ bending moment at $x$  \
$Q^x(w)$ = -$(EIw'')'(x)$ shear force at $x$ \
\
!!!!!!!!!!!Insert picture of the "cantilever" from the Graphics folder!!!!!\
\


"""

# ╔═╡ 00fbb85c-f508-11ea-21b0-1d20d4bdb05d
md"""#### 3. Domain generation
In this section, the domain is being generated based on the parameters selected by the user in Section 2.
"""


# ╔═╡ 16101166-f4de-11ea-25f8-bdfa0e83d2a4
function gen_elements(L,N)
	xh=LinRange(0,L,N)
	ne=N-1
	conn=[[e,e+1] for e=1:ne]
	return xh,ne,conn
end

# ╔═╡ 548dde12-f508-11ea-2ef8-fffbcd320ec8
xh,ne,conn=gen_elements(L,N);

# ╔═╡ c15ccb10-fdb7-11ea-2491-4557ac7265e5
#scatter(xh,5*ones(N),ylims=(0,10),yaxis=nothing,size=(700,150),label="Ω")

# ╔═╡ 2c02f7ba-f4ed-11ea-19c9-9bb69c5de35b
function local2global(ne)
	return [[2*e-1,2*e,2*e+1,2*e+2] for e=1:ne]
end

# ╔═╡ 998838b6-f508-11ea-0949-8b40b202bab8
dofs=local2global(ne);

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

# ╔═╡ 98d94f60-000f-11eb-1e15-6932b56cb721
md"""#### 4. Galerkin methods and Ansatz space
An approximate solution of w, wh needs to be computed. 
To calculate approximate solutions $w_h$ of the bending equation, a finite dimensional subspace $V_h$ (Ansatz space) of the space $V$ needs to be chosen with basis $\phi_1 \dots, \phi_N: [0,L] \rightarrow \mathbb{R}$, where
$w_h(x) = \sum_{k=1}^{N} w_k \phi_k(x), w_k \in R.$
\
Inserting this Ansatz into the weak formulation with Galerkin Bubnov scheme, where $\varphi = \phi_j, j=1,...,N$
\
$\int_{0}^L EIw_h'' \phi_j'' = \int_{0}^L q \phi_j + b(\phi j), j = 1,..., N$\
Then the equations can be written in the matrix-vector form:
S**w** = **q** + $Q_L$**e**$_L$ - $Q_0$**e**$_0$ + $M_L$**d**$_L$ - $M_0$**d**$_0$

**Choice of Ansatz space** \
In order that $\phi''$ is well defined, $\phi$ needs $\phi \in C^1$, i.e. both $\phi$ and $\phi'$ need to be continuous.\
A good choice can be:
$n\geq 2$. Let $h = L/(n-1)$, $x_i = h(i-1)$ and $i = 1,\dots,n$.\

$V_h = \{ \phi \in V | \phi|_{\left(x_i,x_{i+1}\right)} \text{ is a polynomial of degree} \leq 3, i=1,\dots,n-1\}$ \
$V_h$ was chosen to be the space of piecewise cubic polynomials.

Each $\phi \in V_h$ is uniquely determined by the values $u_{2i-1}:= \phi(x_i)$ and $u_{2i} = \phi'(x_i)$ and can be uniquely written as
$\phi = \sum_{k=1}^{2n} u_k \phi_k = \sum_{i=1}^{n} (u_{2i-1}\phi_{2i-1} + u_{2i}\phi_{2i})$
with the basis functions $\phi_i \in V_h$, for $i=1,\dots,2n$, define:\
$\phi_1(x) =\begin{cases}
        \bar{\phi}_1(\frac{x}{h}) & x\in[0,h]\\
        0 & \text{otherwise},
\end{cases}$
$\phi_2(x) = \begin{cases}
h\bar{\phi}_2(\frac{x}{h}) & x\in[0,h]\\
0 & \text{otherwise},
\end{cases}$\
for $i=2,\dots,n-1$:\
$\phi_{2i-1}(x) = \begin{cases}
        \bar{\phi}_3(\frac{x-x_{i-1}}{h}) & x\in[x_{i-1},x_i]\\
        \bar{\phi}_1(\frac{x-x_{i}}{h}) & x\in[x_{i},x_{i+1}]\\
        0 & \text{otherwise},
\end{cases}$
$\phi_{2i}(x) =\begin{cases}
h\bar{\phi}_4(\frac{x-x_{i-1}}{h}) & x\in[x_{i-1},x_i]\\
h\bar{\phi}_2(\frac{x-x_{i}}{h}) & x\in[x_{i},x_{i+1}]\\
0 & \text{otherwise},
\end{cases}$\
where,
\
$\bar{\phi}_1(\xi) = 1 - 3\xi^2 + 2\xi^3$\
$\bar{\phi}_2(\xi) = \xi (\xi - 1)^2$\
$\bar{\phi}_3(\xi) = 3\xi^2 - 2\xi^3$\
$\bar{\phi}_4(\xi) = \xi^2 (\xi - 1)$\
\
In the next section of the report, it is possible to check the approximation of $w_h(x)$ with basis functions chosen as above.


"""

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

# ╔═╡ 3f7a5a08-f03e-11ea-1fcf-83e5b5e57d28
phi,phid,phidd= basis1D(xh,N);

# ╔═╡ 098b7a0a-f2c1-11ea-31d6-91b4b485ac41
md"###### Do you want to check whether our basis functions can approximate a function?
Tick the box below"

# ╔═╡ d8fa3576-fea1-11ea-127e-01b24b313b6b
md"Check $(@bind approximate CheckBox()) "

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

# ╔═╡ 3f0f1d7a-fea2-11ea-1d7d-bffce1d81fec
f(x)=cos(x)*sin(x) ##Give a function

# ╔═╡ c78ce07a-f2bc-11ea-28bd-a1df0c486ed1
if approximate
	check(xh,f::Function);
end

# ╔═╡ 04fb4452-f2c1-11ea-393c-4b9b88d2aaf8
md"### 5. Calculation of element FE matrices and vectors"

# ╔═╡ 1f9d6fd0-002a-11eb-2fc2-87a8a7250e95
md"""In order to compute the approximate solution $w_h$, the finite dimensional subspace $V_h$ with basis functions is used. Inserting the Ansatz into the weak formulation with Galerkin Bubnov scheme, a system of $N$ linear equations with $N$ unknowns $w_1,....,w_N$ is obtained. And this linear system is as follows:
S**w** = **q** + $Q_L$**e**$_L$ - $Q_0$**e**$_0$ + $M_L$**d**$_L$ - $M_0$**d**$_0$

"""

# ╔═╡ bc0f5540-002a-11eb-32d1-558ea518bc1e
md"""where, 
$w = \begin{bmatrix}
            w_1\\
            \vdots\\
            w_N 
\end{bmatrix}$,
                 $\quad 
                 \hat{q}=
                 \begin{bmatrix}
                       \int q\phi_1\\ 
                        \vdots\\ 
                        \int q\phi_N\\
                  \end{bmatrix},
                  \quad$ 
                 $e_x=
                 \begin{bmatrix}
                       \phi_1(x)\\ 
                        \vdots\\
                        \phi_N(x)
                  \end{bmatrix},
                  \quad$
$d_x=
                 \begin{bmatrix}
                       \phi'_1(x)\\ 
                        \vdots\\
                        \phi'_N(x)
                  \end{bmatrix}$
"""

# ╔═╡ 682a04b0-002b-11eb-0c65-ff6186c4585b
md"""$S = \begin{bmatrix}
    s_{11}& \cdots& s_{1N}\\
    \vdots & \ddots& \vdots\\
    s_{N1}& \cdots& s_{NN}\\
\end{bmatrix} \quad \text{ and } \quad s_{jk} = \int^L_0 EI \phi''_k\phi''_j$

$S\in \mathbb{R}^{N \times N}$ is the stiffness matrix and is positive semi-definite for any $w \in \mathbb{R}^n$. \

For all x: $w_h(x)$ = **e**$_x^T$**w**, $w_h'(x)$ = **d**$_x^T$**w**

"""

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

# ╔═╡ 9e66cea0-0030-11eb-26e9-4d81e95c1cbc
md"### 6. Boundary conditions"

# ╔═╡ c59aa2f0-0033-11eb-0b4f-d53534e78f1a
md"""
The global stiffness matrix $S$ that is used to solve the system, is then obtained by an appropriate local to global transformation of the indices.
Note that $S \in \mathbb{R}^{2N\times 2N}$.

On a similar note, the mass matrix which is used to calculate the load vector(also used in the dynamic case as described further) is formulated as follows. 

$\begin{equation*}
M = \begin{bmatrix}
    m_{11}& \cdots& m_{1N}\\
    \vdots & \ddots& \vdots\\
    m_{N1}& \cdots& m_{NN}\\
    \end{bmatrix} \quad \text{ where, } \quad m_{ij} = \int^L_0 \phi_i\phi_j
\end{equation*}$
\
and the load vector is given by $\[\hat{q} = Mq\\\]\$

To get a unique solution of the beam equation, two essential and two natural boundary conditions need to be added.
$\begin{equation*}
    w(0)=a, \quad w'(0)=b, \quad Q^L(w)=Q_L, \quad M^L(w)=M_L
\end{equation*}$

where $a,b,Q_L,M_L \in \mathbb{R}$ are given. Incorporating these into the linear system results in the following extensions for the cantilever beam and for the simply supported beam, respectively. 

"""

# ╔═╡ 88668960-0035-11eb-25b3-c3aa828156fc
md"""
$\begin{bmatrix}
S& e_0& d_0\\
e^T_0& 0&  0\\
d^T_0& 0& 0\\
\end{bmatrix} \begin{bmatrix}
w \\
Q_0\\
M_0\\
                  \end{bmatrix} = \begin{bmatrix}
                                   q+Q_Le_L + M_Ld_L\\
                                   a\\
                                   b\end{bmatrix}$
"""

# ╔═╡ f35502b0-003a-11eb-1d90-e91d7447f1f2
md"""    
$\begin{bmatrix}
    S& e_0& -e_L\\
    e^T_0& 0&  0\\
    -e^T_L& 0& 0\\
    \end{bmatrix} \begin{bmatrix}
                   w \\
                   Q_0\\
                   Q_L\\
                  \end{bmatrix} = \begin{bmatrix}
                                   q-M_0d_0 + M_Ld_L\\
                                   -a_0\\
                                   -a_L
\end{bmatrix}$
"""

# ╔═╡ ed182740-0032-11eb-3715-d938b0f81ec2
md"""

To get a unique solution of the beam equation, two Dirichlet
boundary conditions need to be added. For a left sided clamped beam (cantilever) these conditions are $w(0) =a$; $w'(L) = b$, where a, b $\in R$ are given. The Galerkin approximation should satisfy these
conditions, too and incorporating these into the linear system of equations

"""

# ╔═╡ 957d80c0-003c-11eb-3846-850ad97d46c9
md"##### Dynamic case"

# ╔═╡ dce9f510-003c-11eb-24be-23294f7fac9f
md"""
In the dynamic case the bending curve as well as all forces and boundary conditions depend on time. In particular the bending curve at time t is given by $w(x, t), x \in [0,L]$. The bending equation is given by

$\mu \ddot{w}+\left(E I w^{\prime \prime}\right)^{\prime \prime}=q$

where $\mu = \mu(x)$ is the mass density (more precisely: the mass per unit length) and $\ddot{w}$ denotes the second derivative of $w$ with respect to $t$. During all the simulations, the mass density is assumed to be constant through the beam. One important fact about this differential equation is the absence of a disspation term,without it, the vibration movement of the beam will be non-damped over time. This not only mean that the movement of the beam will not stop, but also that all the vibration modes (especially the higher ones) will be present during the whole simulation. 

The same procedure as in the static case (multiplication with  and partial integration) yields the weak formulation: 

$\int_{0}^{L} \ddot{w} \psi+\int_{0}^{L} E I w^{\prime \prime} \psi^{\prime \prime}=\int_{0}^{L} q(\cdot, t) \psi+b(\psi, t),$\

where 

$b(\psi, t)=Q_{L}(t) \psi(L)-Q_{0}(t) \psi(0)+M_{L}(t) \psi^{\prime}(L)-M_{0}(t) \psi^{\prime}(0).$\

This domain discretization yields the following DAE that will be solved by the Newmark's method: \

$$\left[\begin{array}{ccc}
M & 0 & 0 \\
0 & 0 & 0 \\
0 & 0 & 0
\end{array}\right]\left[\begin{array}{c}
\underline{\ddot{w}} \\
\ddot{Q}_{0} \\
\ddot{M}_{0}
\end{array}\right]+\left[\begin{array}{ccc}
S & \underline{e}_{0} & \underline{d}_{0} \\
\underline{e}_{0}^{\top} & 0 & 0 \\
\underline{d}_{0}^{\top} & 0 & 0
\end{array}\right]\left[\begin{array}{c}
\underline{w} \\
Q_{0} \\
M_{0}
\end{array}\right]=\left[\begin{array}{c}
\underline{q}+Q_{L} \underline{e}_{L}+M_{L} \underline{d}_{L} \\
a \\
b
\end{array}\right]$$\
\
Here $M$ describes the mass matrix of the system, given by $M = [m_{ij}]_{i,j =1}^N$. As for the stiffness matrix $S$, that we derived above, we considered the local mass matrix $M_{loc}$ for the reference interval $[0,h]$. The global mass matrix $M$ will then again be calculated by a transformation from local to global coordinates.
It is $m_{ij} = \int_0^L \mu \phi_i \phi_j$
and 
$
M_{loc} = 
\frac{\mu\cdot A \cdot h}{420}
\begin{bmatrix}
156 & 22h & 54 & -13h\\
22h & 4h^2 & 13h & -3h^2\\
54 & 13h & 156 & -22h\\
-13h &-3h^2 & -22h & 4h^2\\
\end{bmatrix}.
$


"""

# ╔═╡ e1d4afc0-00cd-11eb-2cd3-6fe14acc0080
md"""
**The Newmark's algorithm** consists of : \
$\text { Compute } u_{j}^{*} \text { and } \dot{u}_{j}^{*} \text { using }:$

$$u_{j}^{*} = u_{j}+\dot{u}_{j} h_{j}+\left(\frac{1}{2}-\beta\right) \ddot{u}_{j} h_{j}^{2}$$

$$\dot{u}_{j}^{*} =\dot{u}_{j}+(1-\gamma) \ddot{u}_{j} h_{j}$$

Compute the solution $\ddot{u}_{j+1}$ with: 

$$f\left(\ddot{u}_{j+1}, \dot{u}_{j}^{*}+\gamma \ddot{u}_{j+1} h_{j}, u_{j}^{*}+\beta \ddot{u}_{j+1} h_{j}^{2}, t_{j+1}\right)=0$$\

and finally, compute $u_{j+1}$ and $\dot{u}_{j+1}$ using: \

$$u_{j+1}=u_{j}^{*}+\beta \ddot{u}_{j+1} h_{j}^{2}, \quad \dot{u}_{j+1}=\dot{u}_{j}^{*}+\gamma \ddot{u}_{j+1} h_{j}.$$

The parameter values were chosen to be $\beta = 1/4$ and $\gamma = 1/2$.

"""

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
	M_e = [M 0*C; 0*C' zeros(2,2)];
	
	return S_e,M_e,q_e;
end	

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

# ╔═╡ b59aacb2-f50d-11ea-34dd-67c479511398
md"### Solving the system"

# ╔═╡ 28da6d86-f4f7-11ea-010d-f7838116054e
begin
	if analysis=="static"
		S_e,q_e = assemble(t,EI,Q,ML,qv,a,b,xh,ne,conn,phi,phidd,dofs);
		w = S_e\q_e;
		w_interp = sum((w[2*j-1]*phi[2*j-1].(xh) + w[2*j]*phi[2*j].(xh)) for j=1:N);
		
	elseif analysis=="dynamic"
		S_e,M_e,q_e = dyn_assemble(mu,EI,Q,ML,qv,xh,ne,conn,phi,phidd,dofs);
		w = S_e\q_e;
		w_interp = newmark_step(Mtau,xh,w,tau,beta,gamma,M_e,S_e,phi);	
	end
	nothing
end

# ╔═╡ e2536ecc-f674-11ea-386d-274c40e5a846
md"### Results and Post Processing"

# ╔═╡ 1b080bca-fed5-11ea-169f-5930c4037901
md"###### Let us analyze the results!"

# ╔═╡ b2b57896-f696-11ea-07b0-2f01a202fbd0
function plotFEM(xh,wh,plt_title="Static Case",plt_label="FEM",ylim=[-L/2,L/2])
	#assumption of linear interpolation in between nodes while plotting. 
	plot(xh,wh,
		xlims=(0,L+0.05),
		ylims=(ylim[1],ylim[2]),
		label=plt_label,
		xlabel="Domain (Ω)",
		ylabel="Deflection (w)",
		marker=(:circle,:orange),title=plt_title)
end	

# ╔═╡ a32eccb4-f6ab-11ea-2b2c-1de81ae66813
if analysis=="static"
	md"#### The Static Case"
end

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
			elseif iszero(q)
				w_exact = (3*ML .+ 3*L*Q .- Q*xh).*xh.*xh/6*ei
				text="with end moment and load"
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
	w_exact=w_exact +  a*ones(length(xh))
return w_exact,text;
end

# ╔═╡ ae561da6-fe92-11ea-2d88-03ef3f13dff8
if analysis=="static"
md"make sure you uncheck the options before seeing new plots, otherwise the plot! will overlap on the same image. Julia `Plots()` package artifact :("
end

# ╔═╡ f64ad968-f68f-11ea-2a36-7b7bdac2b87d
if analysis=="static"
		md"""**Check Analytical solution** $(@bind check_sol CheckBox())
	"""
end

# ╔═╡ 68166060-f5b1-11ea-0e76-41f2220403c6
begin
	if analysis=="static"
	y_lim=[minimum(w_interp)-a maximum(w_interp)] 
	y_lim = y_lim+0.5*y_lim
		if check_sol==false
			plotFEM(xh,w_interp,"Static Case","FEM Solution",y_lim)
		else 
			w_ex, ctext = check_analytical(t,a,b,Q,ML,ul,L,ei,xh)
			plot!(xh,w_ex,
				annotations=(L/3,y_lim[2]-0.5*y_lim[2],ctext),
				label="Analytical",
				line=(:dash)
				)
		end
		else
			md"Select the static case to see the static plots!"
	end
end

# ╔═╡ 96f1baee-f6ab-11ea-2dd3-57f8e2e541d8
if analysis=="dynamic"
md"#### The Dynamic Case - Free and Undamped Vibrations"
end

# ╔═╡ 16792298-f6a5-11ea-07a8-39d00521b9db
if analysis=="dynamic"
	md"#### Plotting the first 5 modes"
end

# ╔═╡ 34eec6f6-fd9e-11ea-1563-e1a9c54012ad
if analysis=="dynamic"
	if t=="c"
		idx=3:2*N
	elseif t=="ss"
		part=[i for i in 3:(2*N-1)]
		idx=pushfirst!(part,1)
	end
	
	λ, ϕ = eigs(S_e[idx,idx],M_e[idx,idx]; nev=N,which=:SM);
	eigenpairs = [(λ[i], ϕ[:,i]) for i in 1:length(λ)];
	sort!(eigenpairs, by=pair -> abs(pair[1]))
	nothing
end

# ╔═╡ d2faeeb8-fdff-11ea-3fab-8d961cba8f25
begin
	if analysis=="dynamic"
		no_of_eigs=5
		w_interp_eig=zeros(N,no_of_eigs)
		for eig_i in 1:no_of_eigs
			w_mode = zeros(2*N)
			val, vect = eigenpairs[eig_i]
			vect = vect*sign(vect[5])
			w_mode[idx] = vect[:]
			w_interp_eig[:,eig_i] = sum((w_mode[2*j-1]*phi[2*j-1].(xh) + w_mode[2*j]*phi[2*j].(xh)) for j=1:N)
		end
		if t=="c"
			plotFEM(xh,w_interp_eig,"Mode Shapes-Cantilever",["Mode 1" "Mode 2" "Mode 3" "Mode 4" "Mode 5"])
		elseif t=="ss"
			plotFEM(xh,w_interp_eig,"Mode Shapes-Simply Supported",["Mode 1" "Mode 2" "Mode 3" "Mode 4" "Mode 5"])
		end
	end
end

# ╔═╡ 0af93070-fda3-11ea-06a3-6b34638a8696
if analysis=="dynamic"
	md"**Create Movie** $(@bind click CheckBox()) "
end

# ╔═╡ 5c0cb724-fec9-11ea-17de-817c5a58be19
if click
	if analysis=="dynamic"
		anim = @animate for T ∈ 1:Mtau
			plotFEM(xh,w_interp[:,T],
				"Free vibrations","FEM solution",[-1,1])
			plot!(annotations=(L/3,L/1.15,"T=$(T)"))
		end
		gif(anim)
	else 
		"select the dynamic case to see this clip"
	end
end

# ╔═╡ c94b321e-00d3-11eb-0c26-f1b5472f73cb
md"""
#### References
1. http://www.learnaboutstructures.com/Determinate-Deflections-Introduction
"""

# ╔═╡ Cell order:
# ╟─219fcfc0-ef96-11ea-3baa-85eb1e81cf63
# ╟─c7eb8ae0-fffe-11ea-3f09-318f505b326e
# ╟─e125c7d6-ef97-11ea-39fd-b176103f48ba
# ╟─6646b4f0-fdbe-11ea-38a6-9ff1504f8f5d
# ╟─2ff6cc50-0002-11eb-1d90-e1f3bdada8fc
# ╟─d8e7f8c0-00d4-11eb-3aad-d5b4512646fa
# ╟─3840efd0-fffd-11ea-0445-4184c1500cf5
# ╟─7171003c-f65f-11ea-0ed1-eb09bbaa9a32
# ╟─65d5625e-f58b-11ea-1cb7-fb37ec26e2c3
# ╟─4c1d7992-f668-11ea-01e2-6300503d267d
# ╟─58ce39ee-fe93-11ea-3244-f5d7f3935485
# ╟─f2b09b22-f03a-11ea-2baf-ff3e1721752d
# ╟─66e54c52-0002-11eb-27d1-45419b83cfa6
# ╟─00fbb85c-f508-11ea-21b0-1d20d4bdb05d
# ╟─16101166-f4de-11ea-25f8-bdfa0e83d2a4
# ╟─548dde12-f508-11ea-2ef8-fffbcd320ec8
# ╟─c15ccb10-fdb7-11ea-2491-4557ac7265e5
# ╟─2c02f7ba-f4ed-11ea-19c9-9bb69c5de35b
# ╟─998838b6-f508-11ea-0949-8b40b202bab8
# ╟─2b84ae9c-f66d-11ea-38e4-c1e6e9e7b7b8
# ╟─8cbd99bc-f66d-11ea-3c74-2b9b9999741d
# ╟─98d94f60-000f-11eb-1e15-6932b56cb721
# ╟─466ce928-f03b-11ea-1b7e-a9e26cbccdcc
# ╠═3f7a5a08-f03e-11ea-1fcf-83e5b5e57d28
# ╟─098b7a0a-f2c1-11ea-31d6-91b4b485ac41
# ╟─d8fa3576-fea1-11ea-127e-01b24b313b6b
# ╟─dbb675c8-f2b4-11ea-3ffa-81e83739270d
# ╠═3f0f1d7a-fea2-11ea-1d7d-bffce1d81fec
# ╟─c78ce07a-f2bc-11ea-28bd-a1df0c486ed1
# ╟─04fb4452-f2c1-11ea-393c-4b9b88d2aaf8
# ╟─1f9d6fd0-002a-11eb-2fc2-87a8a7250e95
# ╟─bc0f5540-002a-11eb-32d1-558ea518bc1e
# ╟─682a04b0-002b-11eb-0c65-ff6186c4585b
# ╟─9aea9122-f5d0-11ea-1bc7-d5b5f6368033
# ╟─6c4aa6b8-f39b-11ea-0d7e-6bf87aaa95be
# ╟─8352e10c-f4ce-11ea-25a9-d75c910445fa
# ╟─9e66cea0-0030-11eb-26e9-4d81e95c1cbc
# ╟─c59aa2f0-0033-11eb-0b4f-d53534e78f1a
# ╟─88668960-0035-11eb-25b3-c3aa828156fc
# ╟─f35502b0-003a-11eb-1d90-e91d7447f1f2
# ╟─ed182740-0032-11eb-3715-d938b0f81ec2
# ╟─957d80c0-003c-11eb-3846-850ad97d46c9
# ╟─dce9f510-003c-11eb-24be-23294f7fac9f
# ╟─e1d4afc0-00cd-11eb-2cd3-6fe14acc0080
# ╟─a94fb776-f66e-11ea-205e-21f70a68ef1d
# ╟─0c893646-f66f-11ea-1ca9-e32285b8c42b
# ╟─680184dc-f422-11ea-3aa4-e7d07918b79b
# ╟─84c0f1b8-f5d0-11ea-14e6-c987ddf8da67
# ╟─986616bc-f5dd-11ea-3975-bb8d25062cce
# ╟─b9cb119e-f5e8-11ea-302d-6db9fec26b2b
# ╟─b59aacb2-f50d-11ea-34dd-67c479511398
# ╟─28da6d86-f4f7-11ea-010d-f7838116054e
# ╟─e2536ecc-f674-11ea-386d-274c40e5a846
# ╟─1b080bca-fed5-11ea-169f-5930c4037901
# ╟─b2b57896-f696-11ea-07b0-2f01a202fbd0
# ╟─a32eccb4-f6ab-11ea-2b2c-1de81ae66813
# ╟─1b42f33a-f510-11ea-128a-c5dadf723476
# ╟─ae561da6-fe92-11ea-2d88-03ef3f13dff8
# ╟─f64ad968-f68f-11ea-2a36-7b7bdac2b87d
# ╟─68166060-f5b1-11ea-0e76-41f2220403c6
# ╟─96f1baee-f6ab-11ea-2dd3-57f8e2e541d8
# ╟─16792298-f6a5-11ea-07a8-39d00521b9db
# ╟─34eec6f6-fd9e-11ea-1563-e1a9c54012ad
# ╟─d2faeeb8-fdff-11ea-3fab-8d961cba8f25
# ╟─0af93070-fda3-11ea-06a3-6b34638a8696
# ╟─5c0cb724-fec9-11ea-17de-817c5a58be19
# ╟─c94b321e-00d3-11eb-0c26-f1b5472f73cb
