using SparseArrays
using Arpack
using ForwardDiff
using Plots	
pyplot()

include("domain.jl")
include("FE_utils.jl")

t="ss" #c -> cantilever, ss -> simply supported
eigen=true;
movie=false;

L=1; #Length
N=100; #Nodes
ul=1; #constant distributed load
μ=1; #mass density
ei=1; #bending stiffness

#boundary conditions
Q=1;
ML=0;
a=0;
b=0;

#dynamic parameters
Mtau=200
tau=1
beta=0.25
gamma=0.5

xh,ne,conn = gen_elements(L,N)
qv,mu,EI = const_dist(ul,μ,ei)
dofs = local2global(ne)

phi,phid,phidd = basis1D(xh,N)


S_e,M_e,q_e = dyn_assemble(t,mu,EI,Q,ML,qv,a,b,xh,ne,conn,phi,phidd,dofs);
w = S_e\q_e; #initial condition got from static_solve
w_interp = newmark_step(Mtau,xh,w,tau,beta,gamma,M_e,S_e,phi);

#create movie
if movie
    anim = @animate for T ∈ 1:Mtau
        plotFEM(xh,w_interp[:,T],
	        "Free vibrations","FEM solution",[-1,1])
        plot!(annotations=(L/3,L/1.15,"T=$(T)"))
    end
    gif(anim)
end

#plot eigenmodes and frequencies
if eigen
    if t=="c"
        idx=3:2*N
    elseif t=="ss"
        part=[i for i in 2:(2*N-2)]
        idx=push!(part,2*N)
    end

    λ, ϕ = eigs(S_e[idx,idx],M_e[idx,idx],nev=N,which=:SM);
    sort!(λ, by=val -> abs(val))
    eigenpairs = [(λ[i], ϕ[:,i]) for i in 1:length(λ)];
    sort!(eigenpairs, by=pair -> abs(pair[1]))
    
    no_of_eigs=5
    w_interp_eig=zeros(N,no_of_eigs)
    for eig_i in 1:no_of_eigs
        w_mode = zeros(2*N)
        val, vect = eigenpairs[eig_i]
        vect = vect*sign(vect[5])
        w_mode[idx] = vect[:]
        w_interp_eig[:,eig_i] = sum((w_mode[2*j-1]*phi[2*j-1].(xh) + w_mode[2*j]*phi[2*j].(xh)) for j=1:N)
    end
    plot()
    freq=[round(sqrt(Real(λ[i])),digits=2) for i=1:no_of_eigs]
    
    if t=="c"
	plotFEM(xh,w_interp_eig,"Mode Shapes-Cantilever",
                ["Mode 1, ω=$(freq[1])" "Mode 2, ω=$(freq[2])" "Mode 3, ω=$(freq[3])" "Mode 4, ω=$(freq[4])" "Mode 5 ω=$(freq[5])"])
    elseif t=="ss"
	plotFEM(xh,w_interp_eig,"Mode Shapes-Simply Supported",
                ["Mode 1, ω=$(freq[1])" "Mode 2, ω=$(freq[2])" "Mode 3, ω=$(freq[3])" "Mode 4, ω=$(freq[4])" "Mode 5, ω=$(freq[5])"])
    end
end
