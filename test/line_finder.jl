using RvSpectML
using Test

function gen_fake_normalized_spectra(λ_min::Float64, λ_max::Float64, lines::Vector{Float64})
   λ = collect(λ_min:0.01:λ_max)
   flux = ones(length(λ)) - rand(length(λ))*0.01 #continuum with a bit of noise
   for l in lines
      flux = flux .* (1.0.-0.98*rand()*exp.(-1.0 .* (λ.-l).^2 ./ 2 ./ (λ.*7000.0./3e8).^2)) #remove each line in lines
   end
   λ, flux, sqrt.(flux)
end

@testset "Line_Finder" begin

   lines1 = Float64[]
   λ1, flux1, var1 = gen_fake_normalized_spectra(4960.0,5000.0,lines1)
   chunk1 = RvSpectML.RvSpectMLBase.ChunkOfSpectrum(λ1,flux1,var1)
   lines2 = [5001,5005,5010,5010.3,5011,5011.5,5020,5020.2,5035]
   λ2, flux2, var2 = gen_fake_normalized_spectra(5000.0,5040.0,lines2)
   chunk2 = RvSpectML.RvSpectMLBase.ChunkOfSpectrum(λ2,flux2,var2)
   cl = RvSpectML.RvSpectMLBase.ChunkList([chunk1, chunk2],[1,2])
   @test_nowarn lines = RvSpectML.LineFinder.find_lines_in_chunklist(cl, plan=RvSpectML.LineFinder.LineFinderPlan(line_width=7000.0,min_deriv2=0.5, use_logλ=true, use_logflux=false), verbose=false)

end