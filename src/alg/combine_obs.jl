"""
Code for combining spectra (currently just sums them)

Author: Eric Ford
Created: August 2020
Contact: https://github.com/eford/
"""




"""   bin_times_consecutive( times, n )
Computes mean times from conseuctive bins of n times (to go with bin_consecutive_spectra).
Returns floor(length(times)/n) elements.

WARNING:  Simply takes consecutive times, so some bins may be from spectra that weren't taken close together.
TODO:  Create version that pays attention to timestamps.
"""
function bin_times_consecutive(times::AT, n::Integer) where { T<:Real, AT<:AbstractVector{T} }
    local num_in = length(times)
    @assert num_in >= n
    local num_out = floor(Int,num_in//n)
    local times_out = Vector{eltype(times)}(undef,num_out)
    local time_idx_out  = Array{UnitRange,1}(undef,num_out)
    local idx_start = 1
    for i in 1:num_out
      idx_stop = min(idx_start+n-1,num_in)
      times_out[i] = mean(times[idx_start:idx_stop])
      time_idx_out[i] = idx_start:idx_stop
      idx_start = idx_stop + 1
    end
    return times_out
end

bin_rvs_consecutive(rvs::AT, n::Integer) where { T<:Real, AT<:AbstractVector{T} } = bin_times_consecutive(rvs,n)

"""   bin_spectra_nightly
Bins spectra from a SpectralTimeSeriesCommonWavelengths object with a maximum spacing between observation times
"""
function bin_spectra_nightly(spectra::AbstractSpectralTimeSeriesCommonWavelengths, times::AT) where { T<:Real, AT<:AbstractVector{T} }
    Δt_threshold = 0.5
    bin_spectra_max_Δt(spectra,times,Δt_threshold)
end

function bin_times_nightly(times::AT)  where { T<:Real, AT<:AbstractVector{T} }
    Δt_threshold = 0.5
    bin_times_max_Δt(times,Δt_threshold)
end

function bin_rvs_nightly(;times::AT2, rvs::AT1)  where { T1<:Real, AT1<:AbstractVector{T1}, T2<:Real, AT2<:AbstractVector{T2}  }
    Δt_threshold = 0.5
    bin_rvs_max_Δt(times=times,rvs=rvs,Δt_threshold=Δt_threshold)
end

function bin_times_and_rvs_nightly(;times::AT1, rvs::AT2 )  where { T1<:Real, AT1<:AbstractVector{T1}, T2<:Real, AT2<:AbstractVector{T2}  }
    Δt_threshold = 0.5
    bin_times_and_rvs_max_Δt(times=times, rvs=rvs,  Δt_threshold=Δt_threshold)
end

function bin_times_max_Δt(times::AT, Δt_threshold::Real = 0.5)  where { T<:Real, AT<:AbstractVector{T} }
    local bin_idx = make_bin_indices_for_binning_max_Δt(times,Δt_threshold=Δt_threshold)
    local bin_labels = unique(bin_idx)
    local binned_times = zeros(length(bin_labels))
    for (i,label) in enumerate(bin_labels)
        binned_times[i] = mean(times[findall(isequal(label),bin_idx)])
    end
    return binned_times
end

function bin_rvs_max_Δt(;times::AT1, rvs::AT2, Δt_threshold::Real = 0.5)  where { T1<:Real, AT1<:AbstractVector{T1}, T2<:Real, AT2<:AbstractVector{T2}  }
    local bin_idx = make_bin_indices_for_binning_max_Δt(times,Δt_threshold=Δt_threshold)
    local bin_labels = unique(bin_idx)
    local binned_rvs = zeros(length(bin_labels))
    for (i,label) in enumerate(bin_labels)
        binned_rvs[i] = mean(rvs[findall(isequal(label),bin_idx)])
    end
    return binned_rvs
end

function bin_times_and_rvs_max_Δt(;times::AT1, rvs::AT2, Δt_threshold::Real = 0.5)  where { T1<:Real, AT1<:AbstractVector{T1}, T2<:Real, AT2<:AbstractVector{T2}  }
    local bin_idx = make_bin_indices_for_binning_max_Δt(times,Δt_threshold=Δt_threshold)
    local bin_labels = unique(bin_idx)
    local binned_times = zeros(length(bin_labels))
    local binned_rvs = zeros(length(bin_labels))
    local num_obs_in_bin = zeros(length(bin_labels))
    for (i,label) in enumerate(bin_labels)
        local idx = findall(isequal(label),bin_idx)
        binned_times[i] = mean(times[idx])
        binned_rvs[i] = mean(rvs[idx])
        num_obs_in_bin[i] = length(idx)
    end
    return (times=binned_times, rvs=binned_rvs, num_obs_in_bin=num_obs_in_bin)
end

function make_bin_indices_for_binning_max_Δt(times::AT; Δt_threshold::Real = 0.5)  where { T<:Real, AT<:AbstractVector{T} }
    @assert length(times) >= 1
    n = length(times)
    bin_indices = ones(Int,n)
    bin = 1
    bin_start_time = times[1]
    for i in 2:n
        if times[i]-bin_start_time > Δt_threshold
            bin += 1
            bin_start_time = times[i]
        end
        bin_indices[i] = bin
    end
    return bin_indices
end

"""   bin_spectra_max_Δt
Bins spectra from a SpectralTimeSeriesCommonWavelengths object with a maximum spacing between observation times
"""
function bin_spectra_max_Δt(spectra::AbstractSpectralTimeSeriesCommonWavelengths, times::AT, Δt::Real) where { T<:Real, AT<:AbstractVector{T} }
  local num_in = size(spectra.flux,2)
  @assert length(times) >= 1
  @assert Δt>=zero(Δt)
  local bin_idx = make_bin_indices_for_binning_max_Δt(times,Δt_threshold=Δt)
  local bin_labels = unique(bin_idx)
  local num_out = length(bin_labels)
  @assert num_in >= num_out
  local binned_times = zeros(num_out)
  local flux_out = Array{eltype(spectra.flux),2}(undef,size(spectra.flux,1),num_out)
  local var_out  = Array{eltype(spectra.flux),2}(undef,size(spectra.var ,1),num_out)
  local time_ave = Array{eltype(times),1}(undef,num_out)
  local time_idx_out = Array{UnitRange,1}(undef,num_out)
  for (i,label) in enumerate(bin_labels)
      idx_to_bin = findall(isequal(label),bin_idx)
      time_ave[i] = mean(times[idx_to_bin])
      flux_out[:,i] .= vec(sum(spectra.flux[:,idx_to_bin],dims=2))
      var_out[:,i]  .= vec(sum((spectra.var[:,idx_to_bin]),dims=2))
      if length(idx_to_bin) > 1
          @assert all(idx_to_bin[2:end].-idx_to_bin[1:end-1] .== 1)
      end
      time_idx_out[i] = first(idx_to_bin):last(idx_to_bin)
  end
  metadata_out = MetadataT(:time_idx=>time_idx_out, :bjd=>time_ave)
  return typeof(spectra)(spectra.λ,flux_out,var_out,spectra.chunk_map,spectra.inst,metadata_out)
end

"""   bin_spectra_consecutive
Bins consecutive spectra from a SpectralTimeSeriesCommonWavelengths object

WARNING:  Simply takes consecutive spectra, so some bins may be from spectra that weren't taken close together.
TODO:  Create version that pays attention to timestamps.
"""
function bin_spectra_consecutive(spectra::AbstractSpectralTimeSeriesCommonWavelengths, n::Integer)
  local num_in = size(spectra.flux,2)
  @assert num_in >= n
  local num_out = floor(Int,num_in//n)
  local flux_out = Array{eltype(spectra.flux),2}(undef,size(spectra.flux,1),num_out)
  local var_out  = Array{eltype(spectra.flux),2}(undef,size(spectra.flux,1),num_out)
  local time_idx_out  = Array{UnitRange,1}(undef,num_out)
  #idx_binned = map(i->1+(i-1)*obs_per_bin:obs_per_bin*i,1:num_obs_binned)

  local idx_start = 1
  for i in 1:num_out
    idx_stop = min(idx_start+n-1,num_in)
    flux_out[:,i] .= vec(sum(spectra.flux[:,idx_start:idx_stop],dims=2))
    var_out[:,i]  .= vec(sum((spectra.var[:,idx_start:idx_stop]),dims=2))
    time_idx_out[i] = idx_start:idx_stop
    idx_start = idx_stop + 1
  end
  metadata_out = MetadataT(:time_idx=>time_idx_out)
  return typeof(spectra)(spectra.λ,flux_out,var_out,spectra.chunk_map,spectra.inst,metadata_out)
end

"""   bin_times()
Bins times from a SpectralTimeSeriesCommonWavelengths object using the groupings from time_idx in the metadata.
"""
function bin_times(spectra::AbstractSpectralTimeSeriesCommonWavelengths, times::AT, n::Integer) where { T<:Real, AT<:AbstractVector{T} }
    @assert haskey(spectra.metadata,:time_idx)
    map(idx->mean(times[idx]), spectra.metadata[:time_idx])
end

"""  rms_rv_within_night(times, rvs)
Return RMS of RVs taken within the same night
"""
function rms_rvs_within_night(;times::AbstractVector{T1}, rvs::AbstractVector{T2}) where { T1<:Real, T2<:Real }
    Δt_threshold = 0.5
    bin_idx = make_bin_indices_for_binning_max_Δt(times,Δt_threshold=Δt_threshold)
    bin_labels = unique(bin_idx)
    sum_rms = 0
    sum_weight = 0
    for (i,label) in enumerate(bin_labels)
        idx = findall(isequal(label),bin_idx)
        if length(idx) <= 1    continue    end
        weight = length(idx)-1
        sum_rms += std(rvs[idx])*weight
        sum_weight += weight
    end
    @assert sum_weight > 0
    within_night_rms = sum_rms/sum_weight
    return within_night_rms
end



""" make_template_spectra( chunk_list_timeseries, [ options ] )
Combine portions of spectra in a ChunkListTimeseries into a template stored as a SpectralTimeSeriesCommonWavelengths.
Optionally remove radial velocities before combining spectra based on the :rv_est field in metadata.
# Inputs:
- chunk_list_timeseries
# Optional Parametes:
- remove_rv_est (true)
- oversample_factor (1)
- smooth_factor (1)
- alg (:TemporalGP)
# Outputs
- matrix:  timeseries interpolated to each wavelength in the form of a SpectralTimeSeriesCommonWavelengths
- f_mean:  mean flux
- var_mean: variance of mean flux (WARNING: Not computed accurately yet)
- deriv:  dflux/dlnλ averaged over spectra
- deriv2: d²flux/dlnλ² averaged over spectra
"""
function make_template_spectra( chunk_list_timeseries::ACLT; remove_rv_est::Bool = true, oversample_factor::Real = 1, smooth_factor::Real = 1, alg::Symbol = :TemporalGP ) where { ACLT<:AbstractChunkListTimeseries }
    if remove_rv_est
        @assert all(map(md->haskey(md,:rv_est), chunk_list_timeseries.metadata))
    end
    chunk_λ_grids = map(c->RvSpectML.make_grid_for_chunk(chunk_list_timeseries,c,oversample_factor=oversample_factor, remove_rv_est=remove_rv_est), 1:num_chunks(chunk_list_timeseries) )
    local matrix, f_mean, var_mean, deriv, deriv2
    if remove_rv_est
        if alg == :TemporalGP
            ( matrix, f_mean, var_mean, deriv, deriv2 )  = RvSpectML.pack_shifted_chunk_list_timeseries_to_matrix(chunk_list_timeseries,chunk_λ_grids, alg=alg, smooth_factor=smooth_factor )
        else
            if smooth_factor != one(smooth_factor)   @warn "Interpolation algorithm selected does not currently support smoothing."   end
            ( matrix, f_mean, var_mean, deriv, deriv2 )  = RvSpectML.pack_shifted_chunk_list_timeseries_to_matrix(chunk_list_timeseries,chunk_λ_grids, alg=alg, smooth_factor=smooth_factor )
        end
    else
        if alg == :TemporalGP
            ( matrix, f_mean, var_mean, deriv, deriv2 )  = RvSpectML.pack_chunk_list_timeseries_to_matrix(chunk_list_timeseries,chunk_λ_grids, alg=alg, smooth_factor=smooth_factor )
        else
            if smooth_factor != one(smooth_factor)   @warn "Interpolation algorithm selected does not currently support smoothing."   end
            ( matrix, f_mean, var_mean, deriv, deriv2 )  = RvSpectML.pack_chunk_list_timeseries_to_matrix(chunk_list_timeseries,chunk_λ_grids, alg=alg, smooth_factor=smooth_factor )
        end
    end
    return ( matrix, f_mean, var_mean, deriv, deriv2, chunk_λ_grids )
end
