using RvSpectML
using Test

@testset "Masks" begin
using DataFrames
linelist_df = read_mask_espresso(joinpath(pkgdir(RvSpectML),"data/masks/G2.espresso.mas"))
    @test typeof(linelist_df) <: DataFrame
    @test size(linelist_df,1) == 5551
    #@test size(linelist_df,2) == 5
    @test hasproperty(linelist_df,:lambda)
    @test hasproperty(linelist_df,:lambda_lo)
    @test hasproperty(linelist_df,:lambda_hi)
    @test hasproperty(linelist_df,:weight)

linelist_file = joinpath(pkgdir(RvSpectML),"data/linelists","VALD_Fe1_DP_rejectTelluricSlope0.0_badLineFilterESPRESSO-strict-NEID-BIS_overlapcutoff6e-05_depthcutoff0.05_allowBlends0_wavesReiners_depthssolar_nbin1depth0.mas")
if isfile(linelist_file)
    linelist_df = read_mask_vald(linelist_file)
    @test size(linelist_df,1) == 280
    @test typeof(linelist_df) <: DataFrame
    @test hasproperty(linelist_df,:lambda)
    @test hasproperty(linelist_df,:lambda_lo)
    @test hasproperty(linelist_df,:lambda_hi)
    @test hasproperty(linelist_df,:weight)
end

end
