using RvSpectML
using Test

@testset "Line lists" begin
    using DataFrames
    linelist_df = read_linelist_espresso(joinpath(pkgdir(RvSpectML),"data/masks","G2.espresso.mas"))
        @test typeof(linelist_df) <: DataFrame
        @test size(linelist_df,1) == 5484
        @test hasproperty(linelist_df,:lambda)
        @test hasproperty(linelist_df,:weight)
end


@testset "Wavelength conversions (air<->vacuumb)" begin
    λ = 4000.0
    @test λ ≈ λ_vac_to_air(λ_air_to_vac(λ))
    λ = 6000.0
    @test λ ≈ λ_vac_to_air(λ_air_to_vac(λ))
    λ = 8000.0
    @test λ ≈ λ_vac_to_air(λ_air_to_vac(λ))
end
