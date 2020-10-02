function gen_elements(L,N)
	xh=LinRange(0,L,N)
	ne=N-1
	conn=[[e,e+1] for e=1:ne]
	return xh,ne,conn
end

function local2global(ne)
    dofs = [[2*e-1,2*e,2*e+1,2*e+2] for e=1:ne]
    return dofs
end

#Assumption of constant q,EI,mu throughout domain. Change for different distribution
function const_dist(ul,μ,ei)
	qv = ul*ones(ne);
	EI = ei*ones(ne);
	mu = μ*ones(ne);
	return qv,mu,EI
end

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

	
