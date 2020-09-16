using RvSpectML
using Test

@testset "Utils" begin

    @testset "Util" begin
        @test RvSpectML.calc_doppler_factor(10) ≈ 1.0000000333564094
        @test RvSpectML.predict_intrinsic_stellar_line_width(5780) ≈ 9883.42046054907
        @test RvSpectML.allequal([1,1,1,1])
        @test !RvSpectML.allequal([1,2,3,4])

        res = RvSpectML.findargminmax([10.0, 30.0, -20.0, 40.0])
        @test res.min == -20.0
        @test res.max ==  40.0
        @test res.argmin ==  3
        @test res.argmax == 4

        @test RvSpectML.searchsortednearest(-cos.(π*(1:16)./16), 0.0) == 8
    end

    @testset "Spectra" begin
        @test RvSpectML.calc_snr(8.0,4.0) ≈ 4
        @test RvSpectML.calc_snr(8*ones(100),4*ones(100)) ≈ 40
    end

    @testset "Chunks" begin
    end

end
