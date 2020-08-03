import Pkg
Pkg.activate(@__DIR__)

import Pluto
println("From the Pluto main menu, open the file:")
println(joinpath(@__DIR__, "Beam.jl"))
println()
Pluto.run()