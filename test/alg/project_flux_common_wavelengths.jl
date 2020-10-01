
using RvSpectMLBase
using Experimental
using Test

@testset "project_flux_common_wavelengths" begin
    len = 2000
    nspec = 100
    lgx = range(log(5000.0),stop=log(5100),length=len)
    x = exp.(lgx)
    P = 0.5*log(maximum(x)/minimum(x))
    y1 = sin.(2π.*lgx./P)
    #dx = vcat(x[2]-x[1], 0.5*(x[3:end].-x[1:end-2]), x[end]-x[end-1])
    v1 = 0.001
    yo1 = y1 .+ v1.*randn(length(y1))
    yo2 = mapreduce(z->y1 .+ v1.*randn(length(y1)),hcat,1:nspec)
    v2 = v1*ones(size(yo2))
    myo2 = Experimental.calc_mean_spectrum(yo2,v2)
    @assert sum(abs2.((myo2 .- y1)./v1)) < 3*length(y1)/sqrt(length(y1))

end
