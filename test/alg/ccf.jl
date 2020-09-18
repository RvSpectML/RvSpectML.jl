using RvSpectML
using Test

@testset "CCF" begin
    # Generate simuliated data & mask
    wave = range(5434.0, stop=5435.0, length=1000)
    spectrum = RvSpectML.absorption_line.(wave, mid=5434.5, depth=0.75, width=0.05)

    mask_shape = CCF.TopHatCCFMask(TheoreticalInstrument1D())
    line_list = CCF.BasicLineList([5434.5], [0.75])
    ccf_plan = CCF.BasicCCFPlan(mask_shape = mask_shape, line_list=line_list)
    # Get the velocities the CCF will be evaluated at
    v_grid = CCF.calc_ccf_v_grid(ccf_plan)

    # measure the CCF at each velocity
    ccf = RvSpectML.CCF.ccf_1D(wave, spectrum, ccf_plan)


    fit_gaussian_to_ccf = RVFromCCF.MeasureRvFromCCFGaussian()
    #rv1 = RvSpectML.RVFromCCF.measure_rv_from_ccf(ccfd, ccf, mask, fit_type=:gaussian)[1]
    rv1 = RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid, ccf, alg=fit_gaussian_to_ccf)[1]
    println(rv1)
    @test abs(rv1) ≤ 0.04

    fit_quadratic_to_ccf = RVFromCCF.MeasureRvFromCCFQuadratic()
    #=
    #rv2 = RvSpectML.RVFromCCF.measure_rv_from_ccf(ccfd, ccf, mask, fit_type=:quadratic)[1]
    rv2 = RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid, ccf, alg=fit_quadratic_to_ccf)[1]
    println(rv2)
    #@test abs(rv2) ≤ 0.04
    #rv3 = RvSpectML.RVFromCCF.measure_rv_from_ccf(ccfd, ccf, mask, fit_type=:centroid)[1]
    fit_centroid_to_ccf = RVFromCCF.MeasureRvFromCCFCentroid()
    rv3 = RvSpectML.RVFromCCF.measure_rv_from_ccf(v_grid, ccf, alg=fit_centroid_to_ccf)[1]
    println(rv3)
    =#
end
