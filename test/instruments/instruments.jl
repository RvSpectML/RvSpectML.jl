@testset "Instrument Traits" begin
include("expres.jl")
#include("harpsn.jl")
include("neid.jl")
include("theory.jl")
end

@testset "Code shared by Instruments" begin

    #include("common.jl")    # not much to test
    #include("io.jl")        # on hold, since would need to include test FITS file(s) in package
    include("linelists.jl")
    include("masks.jl")

end
