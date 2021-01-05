import Pkg
Pkg.activate(@__DIR__)

import PlutoUI
import SparseArrays
import LinearAlgebra
import Arpack
import ForwardDiff
import PyPlot

import Pluto
println("From the Pluto main menu, open the file:")
println(joinpath(@__DIR__, "Beam.jl"))
println()
Pluto.run()
