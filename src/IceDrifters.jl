module IceDrifters

const DATAPATH = "data/"

include("io.jl")
include("utils.jl")
include("buoy_calculations.jl")
include("dist2coast.jl")
include("deform.jl")
include("spectra.jl")
include("plots.jl")

end
