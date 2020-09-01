using RvSpectML
using Test

@testset "RvSpectML.jl" begin
    include("util.jl")
    include("instruments/instruments.jl")
    include("alg/alg.jl")
end
