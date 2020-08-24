using RvSpectML
using Test

@testset "RvSpectML.jl" begin
    include("util.jl")
    include("alg/interp.jl")
    include("instruments/neid.jl")
end
