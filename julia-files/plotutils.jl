function plotFEM(xh,wh,plt_title="Static Case",plt_label="FEM",ylim=[-L/2,L/2])
	#assumption of linear interpolation in between nodes while plotting. 
	plot!(xh,wh,
		xlims=(0,L+0.05),
		ylims=(ylim[1],ylim[2]),
		label=plt_label,
		xlabel="Domain (Ω)",
		ylabel="Deflection (w)",
		marker=(:circle),title=plt_title)
end

function plotcompare!(xh,wh,w_exact,plt_title,ylim=[-L/2,L/2])
    plot(xh,wh,
	 xlims=(0,L+0.05),
	 ylims=(ylim[1],ylim[2]),
	 label= "FEM",
	 xlabel="Domain (Ω)",
	 ylabel="Deflection (w)",
	 marker=(:circle,:orange),title=plt_title)

    plot!(xh,w_exact,
      label="Analytical",
      line=(:dash)
          )
end
