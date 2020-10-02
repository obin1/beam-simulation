using SparseArrays
using Arpack
using ForwardDiff
using Plots	
pyplot()

include("domain.jl")
include("FE_utils.jl")
include("plotutils.jl")
include("analytical.jl")

t="ss" #c -> cantilever, ss -> simply supported

L=1; #Length
N=10; #Nodes
ul=0; #constant distributed load
ei=1; #bending stiffness

Q=1;#[-1,-2,-3,-4,-5]; #in case of simply supported this is moment
ML=0;
a=1;
b=0;

xh,ne,conn = gen_elements(L,N)

dofs = local2global(ne)

phi,phid,phidd = basis1D(xh,N)
plot()

qv,mu,EI = const_dist(ul,μ,ei)
S_e,q_e = assemble(t,EI,Q,ML,qv,a,b,xh,ne,conn,phi,phidd,dofs);
w = S_e\q_e;
w_interp = sum((w[2*j-1]*phi[2*j-1].(xh) + w[2*j]*phi[2*j].(xh)) for j=1:N);

w_exact,text=check_analytical(t,a,b,Q,ML,ul,L,ei,xh);

#Plotting
y_lim=[minimum(w_interp)-a maximum(w_interp)] 
y_lim = y_lim+0.5*y_lim
plotFEM(xh,w_interp,"Static Case","FEM Solution",y_lim)
#plotcompare(xh,w_interp,w_exact,"EndLoad")
