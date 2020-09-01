import RvSpectML
using Test

@testset "EXPRES" begin
    using RvSpectML.EXPRES
        @testset "1D Extracted Spectra Traits exist" begin
            n1 = EXPRES1D()
            @test typeof(n1) <: AbstractInstrument
        end
        @testset "2D Extracted Spectra Traits exist" begin
            n2 = EXPRES2D()
            @test typeof(n2) <: AbstractInstrument
            @test min_order(n2) == 1
            @test max_order(n2) > 1
            @test min_pixel_in_order(n2) == 1
            @test max_pixel_in_order(n2) > 1
            @test length(orders_to_use_default(n2)) >= 1
            @test min_col_default(n2,1) >=1
            @test max_col_default(n2,1) <= 10000
            @test length(metadata_symbols_default(n2)) >= 1
            @test length(metadata_strings_default(n2)) >= 1
        end
end
