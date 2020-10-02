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

function dyn_assemble(t,mu,EI,Q,ML,qv,a,b,xh,ne,conn,phi,phidd,dofs)
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
