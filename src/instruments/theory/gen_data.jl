
#import ..RvSpectML:   # Already in module
import ..RvSpectML: predict_intrinsic_stellar_line_width, speed_of_light_mps

function calc_λs(inst::AnyTheoreticalInstrument)
   Δlnλ = 1.0/inst.resolution
   n = ceil(Int,log(inst.λ_min)/log(inst.λ_max) / Δlnλ )
   r = range(log(inst.λ_min), stop=log(inst.λ_max), length=n)
   return exp.(r)
end

function calc_λ(pixel, inst::AnyTheoreticalInstrument)
   @assert min_pixel(inst) 1<= pixel <= max_pixel(inst)
   Δlnλ = 1.0/inst.resolution
   logλ_min = log(inst.λ_min)
   logλ_max = log(inst.λ_max)
   n = ceil(Int,(logλ_max-logλ_min) / Δlnλ )
   r = range(logλ_min, stop=logλ_max, n)
   return exp(r[pixel])
end


function generate_spectrum(line_list::DataFrame, inst::AnyTheoreticalInstrument; time::Real = 0.0, rv::Real = 0.0, ssbz::Real = 0.0,
               snr_per_pixel::Real = 1000.0, line_width::Real = predict_intrinsic_stellar_line_width(5780), add_noise::Bool = true ) # where { LLT<:AbstractLineList }
   @assert hasproperty(line_list,:lambda)
   @assert hasproperty(line_list,:weight)
   λ = calc_λs(inst)
   flux = ones(size(λ))
   #var = ones(size(λ))
   doppler_factor = calc_doppler_factor(rv) * (1+ssbz)  # TODO: Fix problem... Is ssbz in m/s or v/c?
   for i in 1:size(line_list,1)
      #width = doppler_factor* line_list.lambda[i]*line_width*1000/speed_of_light_mps   #TODO PUT PACK just for testing
      width = doppler_factor*line_list.lambda[i]*line_width*1000/speed_of_light_mps
      flux .*= RvSpectML.absorption_line.(λ, mid=line_list.lambda[i]*doppler_factor, depth=line_list.weight[i], width=width)
   end
   if add_noise
      flux .+= randn(size(flux)).*sqrt.(flux)./snr_per_pixel
      flux[flux .< 0.0] .= 0.0
   end
   var = flux ./ snr_per_pixel^2
   metadata = Dict{Symbol,Any}( :rv_true=>rv, :snr_per_pixel=>snr_per_pixel, :line_width=>line_width, :bjd=>time, :target=>"Simulation", :ssbz=>ssbz)
   Spectra1DBasic(λ,flux,var,inst, metadata)
end

function generate_spectra_timeseries(times::AbstractArray, line_list::DataFrame, inst::AnyTheoreticalInstrument, rvs::AbstractArray; ssbzs::AbstractArray = zeros(length(rvs)),
               snr_per_pixel::Real = 1000.0, line_width::Real = predict_intrinsic_stellar_line_width(5780) ) # where { LLT<:AbstractLineList }
  @assert length(times) == length(rvs) == length(ssbzs)
  @assert hasproperty(line_list,:lambda)
  @assert hasproperty(line_list,:weight)
  map(i->generate_spectrum(line_list, inst, snr_per_pixel=snr_per_pixel, line_width=line_width, time=times[i], rv=rvs[i], ssbz=ssbzs[i] ), 1:length(times) )
end

function generate_spectrum_line()
   inst = TheoreticalInstrument1D(λ_min=5495,λ_max=5505)
   ll = DataFrame(:lambda=>[5500.0], :weight=>0.6 )
   generate_spectrum(ll, inst)
end
