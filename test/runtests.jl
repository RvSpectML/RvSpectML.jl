using RvSpectML
using Test

@testset "RvSpectML.jl" begin
    @testset "Utils" begin
        @test RvSpectML.calc_doppler_factor(10) ≈ 1.0000000333564094
        @test RvSpectML.calc_snr(16*ones(100),4*ones(100)) ≈ 20
    end
    @testset "Masks" begin
        using DataFrames
        mask_df = read_mask_espresso("../data/masks/G2.espresso.mas")
        @test typeof(mask_df) <: DataFrame
        @test size(mask_df,1) >= 5551
        @test size(mask_df,2) == 2
    end

    include("neid.jl")
end
