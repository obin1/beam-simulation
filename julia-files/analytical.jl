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
				w_exact = (q*xh./(24*ei)).*(L^3 .- 2*L*xh.*xh .+ xh.^3)
				text="with distibuted load"
			elseif iszero(q) && ML==Q
				w_exact = (-ML*xh).*(L.-xh)/(2*ei)
				text="with equal end moments"
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
