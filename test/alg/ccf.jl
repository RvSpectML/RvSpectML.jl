using RvSpectML
using Test

@testset "CCF" begin
    # Generate simuliated data & mask
    wave = range(5434.0, stop=5435.0, length=1000)
    spectrum = RvSpectML.absorption_line.(wave, mid=5434.5, depth=0.75, width=0.05)
    mask = [5434.49 5434.51 0.75]      # make a mask

    # measure the CCF
    ccfd = RvSpectML.CCF.calc_ccf_Δv_grid()
    ccf = RvSpectML.CCF.ccf_1D(wave, spec, mask)

    rv1 = RvSpectML.RVFromCCF.measure_rv_from_ccf(ccfd, ccf, mask, fit_type=:gaussian)[1]
    println(rv1)
    rv2 = RvSpectML.RVFromCCF.measure_rv_from_ccf(ccfd, ccf, mask, fit_type=:quadratic)[1]
    println(rv2)
    rv3 = RvSpectML.RVFromCCF.measure_rv_from_ccf(ccfd, ccf, mask, fit_type=:centroid)[1]
    println(rv3)
    @test abs(rv1) ≤ 0.04
    @test abs(rv2) ≤ 0.04
end
