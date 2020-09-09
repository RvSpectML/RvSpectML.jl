module LineFinder

using ..RvSpectML
using DataFrames, Query
using Statistics
using LsqFit
#using Optim
#using Plots
#using ThreadedIterables
using ..TemporalGPInterpolation


default_min_deriv2 = 1.5e5
default_smooth_factor= 4.0
default_min_pixels_in_line = 5
default_use_logλ = true
default_use_logflux = true

struct LineFinderPlan
  min_deriv2::Float64
  smooth_factor::Float64
  min_pixels_in_line::Int
  use_logλ::Bool
  use_logflux::Bool
end

function LineFinderPlan(; min_deriv2::Real = default_min_deriv2, smooth_factor::Real = default_smooth_factor,
                          min_pixels_in_line::Int = default_min_pixels_in_line,
                          use_logλ::Bool = default_use_logλ, use_logflux::Bool = default_use_logflux )
  LineFinderPlan(min_deriv2, smooth_factor, min_pixels_in_line, use_logλ, use_logflux)
end

"""  find_lines_candidates_in_chunk( chunk, plan )
Convenience function to find links in one chunk of spectrum.
# Inputs:
- chunk
# Optional Arguments:
# Return:
- array of ranges of pixels wihtin chunk
"""
function find_line_candidates_in_chunk(chunk::AbstractChuckOfSpectrum, deriv2::AbstractVector{T}; plan::LineFinderPlan = LineFinderPlan() ) where { T<:Real }
  idx_d2_gt_0 = findall(d2->isless(zero(d2),d2), deriv2)
  line_candidates = UnitRange[]
  first_idx_in_line = last_idx_in_line = first(idx_d2_gt_0)
    for i in idx_d2_gt_0[2:end]
      if i-last_idx_in_line > 1   # No longer in last set of contiguous indices, process last line candidate
        if last_idx_in_line - first_idx_in_line + 1 < plan.min_pixels_in_line
          # noop continue
        elseif maximum(deriv2[first_idx_in_line:last_idx_in_line]) < plan.min_deriv2
          # noop continue
        else  # Add last range of pixels as line candidate
          push!(line_candidates,first_idx_in_line:last_idx_in_line)
        end
        first_idx_in_line = i
        last_idx_in_line = i
      else # ! i-last_idx_in_line>1   # Still in same line
        last_idx_in_line = i
      end
    end # for i in idx_d2_gt_0
    # See if we were in the middle of a line candidate at the end of the chunk
    if last_idx_in_line - first_idx_in_line + 1 < plan.min_pixels_in_line
      # noop continue
    elseif maximum(deriv2[first_idx_in_line:last_idx_in_line]) < plan.min_deriv2
      # noop continue
    else
      push!(line_candidates,first_idx_in_line:last_idx_in_line)
    end
  pixel_to_report = map(r-> r[argmax(deriv2[r])], line_candidates)

  return DataFrame(:pixels=>line_candidates, :pixel_max_d2fdlnλ2=>pixel_to_report )
  #return line_candidates
end

# Utility functions for fitting.  Don't assume will remain unchanged in future versions.
function pack(; a::T1, b::T2, λc::T3, depth::T4, σ²::T5) where { T1<:Real, T2<:Real, T3<:Real, T4<:Real, T5<:Real  }
  return [a, λc-6500, log(depth), log(σ²), b]
end
function unpack(θ::AbstractVector{T} ) where {T<:Real}
  return (a=θ[1], λc=θ[2]+6500, depth=exp(θ[3]), σ²=exp(θ[4]), b=θ[5]  )
end

function line_model(λ::Union{T1,AbstractVector{T1} }, θ::AbstractVector{T2}) where { T1<:Real, T2<:Real}
  a, λc, d, σ², b = unpack(θ)
  return @. a*(1.0+b*(λ-λc))*(1.0 - d * exp(-0.5*(λ-λc)^2 /σ²) )
end


function fit_line(λ::AbstractArray{T1,1}, flux::AbstractArray{T2,1}, var::AbstractArray{T3,1}, idx::UnitRange) where { T1<:Real, T2<:Real, T3<:Real }
  @assert length(λ) == length(flux) == length(var)
    λc0 = λ[idx[argmin(flux[idx])]]
  continuum0 = maximum(flux[idx])
  depth0 = 1-minimum(flux[idx])/continuum0
  width0 = 0.05^2   # might need to make more general to deal with different v sin i's
  slope0 = 0.0
  θ0 = pack(a=continuum0,b=slope0,λc=λc0,depth=depth0,σ²=width0)
  lower_bound = pack(a=0.0, b=-1.0, λc=minimum(λ[idx]), depth=0.0, σ²=0.01^2)
  upper_bound = pack(a=2.0, b= 1.0, λc=maximum(λ[idx]), depth=1.0, σ²=0.15^2)
  try
    res = curve_fit(line_model, λ[idx], flux[idx], 1.0 ./(var[idx]), θ0, lower=lower_bound, upper=upper_bound, show_trace=false, autodiff=:forwarddiff)
    covar = estimate_covar(res)
    param = unpack(res.param)
    chi2perdof = sum(abs2.(res.resid).*res.wt) / length(idx)
    return (param=param, covar=covar, χ²_per_dof=chi2perdof, is_converged = res.converged)
  catch
    param = (a=continuum0, λc=λc0, depth=depth0, σ²=width0, b=slope0  )
    covar = Array{Float64,2}(undef,0,0)
    return (param=param, covar=covar , χ²_per_dof=Inf, is_converged = false)
  end
end


""" fit_line( λ, flux, var )
Fits a basic Gaussian absorption line times a line (variable slope) to data.
Returns a tuple of (param, χ²_per_dof, is_converged)
Warning:  Has some hardcoded parameters that likely need to be generalized for different instruments.
"""
function fit_line(λ::AbstractArray{T1,1}, flux::AbstractArray{T2,1}, var::AbstractArray{T3,1} ) where { T1<:Real, T2<:Real, T3<:Real }
  @assert length(λ) == length(flux) == length(var)
  @warn "Need to check why this results in fits not converging!"
  λc0 = λ[argmin(flux)]
  continuum0 = maximum(flux)
  depth0 = 1-minimum(flux)/continuum0
  width0 = 0.05^2   # might need to make more general to deal with different v sin i's
  slope0 = 0.0
  θ0 = pack(a=continuum0,b=slope0,λc=λc0,depth=depth0,σ²=width0)
  lower_bound = pack(a=0.0, b=-1.0, λc=minimum(λ), depth=0.0, σ²=0.01^2)
  upper_bound = pack(a=2.0, b= 1.0, λc=maximum(λ), depth=1.0, σ²=0.15^2)
  try
    #res = curve_fit(line_model, λ, flux, 1.0 ./(var), θ0, lower=lower_bound, upper=upper_bound, show_trace=false, autodiff=:forwarddiff)
    λtmp = copy(λ)
    ftmp = copy(flux)
    invvar = 1.0 ./ var
    res = curve_fit(line_model, λtmp, ftmp,invvar, θ0, lower=lower_bound, upper=upper_bound, show_trace=false, autodiff=:forwarddiff)
    #res = curve_fit(line_model, λ, flux, 1.0 ./(var), θ0, lower=lower_bound, upper=upper_bound, show_trace=false, autodiff=:forwarddiff)
    param = unpack(res.param)
    covar = estimate_covar(fit)
    chi2perdof = sum(abs2.(res.resid).*res.wt) / length(λ)
    return (param=param, covar=covar, χ²_per_dof=chi2perdof, is_converged = res.converged)
  catch
    param = (a=continuum0, λc=λc0, depth=depth0, σ²=width0, b=slope0  )
    covar = Array{Float64,2}(undef,0,0)
    return (param=param, covar=covar, χ²_per_dof=Inf, is_converged = false)
  end
end

"""  find_lines_in_chunk( chunk, plan )
Convenience function to find links in one chunk of spectrum.
# Inputs:
- chunk
# Optional Arguments:
# Return:
- line_list
"""
function find_lines_in_chunk(chunk::AbstractChuckOfSpectrum; plan::LineFinderPlan = LineFinderPlan(), chunk_id::Integer = 0 #=, telluric_model::AbstractArray{T,1}=zeros(0)=# ) where {T <:Real}

  function fit_line(idx::UnitRange)
    #RvSpectML.LineFinder.fit_line(view(chunk.λ,idx), view(flux,idx), view(chunk.var, idx) )   # some bug in this version of function causes non-convergence
    RvSpectML.LineFinder.fit_line( chunk.λ, flux, chunk.var, idx)
  end

  # Fit flux with extra smoothing for measuring derivatives well
  (flux_smooth, deriv_smooth, deriv2_smooth) = predict_mean_and_derivs(chunk.λ, chunk.flux, chunk.λ, sigmasq_obs=chunk.var,
                                  smooth_factor=plan.smooth_factor, use_logx=plan.use_logλ, use_logy=plan.use_logflux )

  # Find pixel ranges that might be lines
  line_candidates = find_line_candidates_in_chunk(chunk, deriv2_smooth, plan=plan)
  if chunk_id != 0    line_candidates[!,:chunk_id] .= chunk_id  end
  #line_candidates[!,:pixels] = line_candidates[!,:pixels]
  line_candidates[!,:fit_min_λ] = chunk.λ[first.(line_candidates[!,:pixels])]
  line_candidates[!,:fit_max_λ] = chunk.λ[ last.(line_candidates[!,:pixels])]

  # Refit flux without extra smoothing, so depths and widths won't be biased by extra smoothing
  (flux, deriv, deriv2) = predict_mean_and_derivs(chunk.λ, chunk.flux, chunk.λ, sigmasq_obs=chunk.var,
                                  smooth_factor=1, use_logx=plan.use_logλ, use_logy=plan.use_logflux )

  # Compute model fit parameters and add results to dataframe
  result = map(fit_line, line_candidates.pixels)
  #line_candidates[!,:fit_param] = [result[i].param for i in 1:length(result) ]
  line_candidates[!,:fit_λc] = [result[i].param.λc for i in 1:length(result) ]
  line_candidates[!,:fit_depth] = [result[i].param.depth for i in 1:length(result) ]
  line_candidates[!,:fit_σ²] = [result[i].param.σ² for i in 1:length(result) ]
  line_candidates[!,:fit_a] = [result[i].param.a for i in 1:length(result) ]
  line_candidates[!,:fit_b] = [result[i].param.b for i in 1:length(result) ]
  line_candidates[!,:fit_covar] = [result[i].covar for i in 1:length(result) ]
  line_candidates[!,:fit_χ²_per_dof] = [result[i].χ²_per_dof for i in 1:length(result) ]
  line_candidates[!,:fit_converged] = [result[i].is_converged for i in 1:length(result) ]
  line_candidates[!,:λ_near_center] = map(i->chunk.λ[line_candidates[i,:pixel_max_d2fdlnλ2]], 1:size(line_candidates,1) )
  line_candidates[!,:f_near_center] = map(i->flux_smooth[line_candidates[i,:pixel_max_d2fdlnλ2]], 1:size(line_candidates,1) )
  line_candidates[!,:dfdlnλ_near_center] = map(i->deriv_smooth[line_candidates[i,:pixel_max_d2fdlnλ2]], 1:size(line_candidates,1) )
  line_candidates[!,:df2dlnλ2_near_center] = map(i->deriv2_smooth[line_candidates[i,:pixel_max_d2fdlnλ2]], 1:size(line_candidates,1) )
  line_candidates[!,:λ_min_flux] = line_candidates[!,:λ_near_center].-line_candidates[!,:dfdlnλ_near_center]./(line_candidates[!,:df2dlnλ2_near_center])
  line_candidates[!,:dfdlnλ_min_flux] = line_candidates[!,:dfdlnλ_near_center].+log.(line_candidates[!,:λ_min_flux]./line_candidates[!,:λ_near_center]).*line_candidates[!,:df2dlnλ2_near_center]

  lines = line_candidates |> @filter(_.fit_converged) |> @filter(_.fit_σ²<0.018) |> @filter(-0.5 <_.fit_b <0.5) |> @filter(_.fit_a < 1.8) |> DataFrame

  return lines
end

"""  find_lines_in_chunklist ( chunklist, line_list )
Convenience function to find lines in each chunk of a spectrum.
# Inputs:
- chunklist
# Optional Arguments:
# Return:
- line_list
"""
function find_lines_in_chunklist(chunk_list::AbstractChunkList ; plan::LineFinderPlan = LineFinderPlan() )
  mapreduce(i->find_lines_in_chunk(chunk_list[i], plan=plan, chunk_id=i), append!, 1:length(chunk_list) )
end

function find_lines_in_chunklist(chunk_list::AbstractChunkList, spectra::AS ; plan::LineFinderPlan = LineFinderPlan() ) where {AS<:AbstractSpectra}
  if haskey(spectra.metadata,:tellurics)
    return mapreduce(i->find_lines_in_chunk(chunk_list[i], plan=plan, chunk_id=i,
      telluric_model = view(spectra.metadata[:tellurics], chunk_list[i].λ.indices[1], chunk_list[i].λ.indices[2])  ),
        append!, 1:length(chunk_list) )
  else
    return find_lines_in_chunklist(chunk_list,plan=plan)
  end
end


"""  find_lines_in_chunklist_timeseries( chunklist_timeseries, line_list )
Convenience function to find lines in each chunk of each spectrum in a timeseries.
# Inputs:
- chunklist_timeseries
# Optional Arguments:
# Return:
- line_list

"""
function find_lines_in_chunklist_timeseries(clt::AbstractChunkListTimeseries ; plan::LineFinderPlan = LineFinderPlan() )
    # TODO? switch to creating one big dataframes?
    #@threaded
    map(cl->find_lines_in_chunklist(cl, plan=plan),clt.chunk_list)
end


function find_pixels_for_line_in_chunk( chunk::AbstractChuckOfSpectrum, λ_min::Real, λ_max::Real )# ; plan::LineFinderPlan = LineFinderPlan() )
  idx_lo = searchsortedfirst(chunk.λ, λ_min, by=x->x>=λ_min)
  idx_tmp = searchsortedlast(chunk.λ[idx_lo:end], λ_max, by=x->x<=λ_max, rev=true)
  idx_hi = idx_lo + idx_tmp - 1
  return idx_lo:idx_hi
end

function find_pixels_for_line_in_chunklist( chunk_list::AbstractChunkList, λ_min::Real, λ_max::Real; verbose::Bool = false )
  ch_idx_all = findall(c-> (λ_min <= minimum(chunk_list.data[c].λ)) && (maximum(chunk_list.data[c].λ) <= λ_max) ,1:length(chunk_list))
  #map(c->(chunk_idx=c, pixels=find_pixels_for_line_in_chunk(chunk_list.data[c], λ_min, λ_max) ), ch_idx)
  ch_idx = 0
  if length(ch_idx_all) > 1
    snr_of_chunks_with_line = map(c->RvSpectML.calc_snr(chunk_list.data[c].flux, chunk_list.data[c].var), ch_idx_all)
    ch_idx_to_keep = argmax(snr_of_chunks_with_line)
    ch_idx = ch_idx_all[ch_idx_to_keep]
    if verbose
      println(" Found λ=",λ_min,"-",λ_max," in chunks: ", ch_idx_all)
      println(" SNRs = ", snr_of_chunks_with_line)
      println(" Keeping chunk #",ch_idx)
    end
  elseif length(ch_idx_all) == 1
    ch_idx = first(ch_idx_all)
  end
  if ch_idx == 0
    error("Didn't find λ = " *string(λ_min)*" - " *string(λ_max)* " in chunklist.")

  end
  return (chunk_idx=ch_idx, pixels=find_pixels_for_line_in_chunk(chunk_list.data[ch_idx], λ_min, λ_max) )
end

function find_pixels_for_line_in_chunklist( chunk_list::AbstractChunkList, λ_min::Real, λ_max::Real, chunk_id::Integer)
  return (chunk_idx=chunk_id, pixels=find_pixels_for_line_in_chunk(chunk_list.data[chunk_id], λ_min, λ_max) )
end

function find_pixels_for_all_lines_in_chunklist( chunk_list::AbstractChunkList, lines::DataFrame)
  @assert hasproperty(lines,:fit_min_λ)
  @assert hasproperty(lines,:fit_max_λ)
  @assert hasproperty(lines,:chunk_id)
  map(l->find_pixels_for_line_in_chunklist(chunk_list[1], lines[l,:fit_min_λ],  lines[l,:fit_max_λ], lines[l,:chunk_id] ), 1:size(lines,1) )
end


#function fit_line_in_chunklist_timeseries(clt::AbstractChunkListTimeseries, lines::DataFrame, line_idx::Integer)
function fit_line_in_chunklist_timeseries(clt::AbstractChunkListTimeseries, λmin::Real, λmax::Real, chid::Integer)
  df = DataFrame()
  #λmin = lines[line_idx,:fit_min_λ]
  #λmax = lines[line_idx,:fit_max_λ]
  #chid = lines[line_idx,:chunk_id]
  nobs = length(clt.chunk_list)
  ParamT = NamedTuple{(:a, :λc, :depth, :σ², :b),NTuple{5,Float64}}  # maybe generalize by making part of plan?
  param = Vector{ParamT}(undef,nobs)
  fit_a = Vector{Float64}(undef,nobs)
  fit_λc = Vector{Float64}(undef,nobs)
  fit_depth = Vector{Float64}(undef,nobs)
  fit_σ² = Vector{Float64}(undef,nobs)
  fit_b = Vector{Float64}(undef,nobs)
  fit_covar = Vector{Array{Float64,2}}(undef,nobs)
  fit_converged = falses(nobs)
  pixels = Vector{UnitRange{Int64}}(undef,nobs)
  gof = zeros(nobs)
  #obs_idx = collect(1:nobs)
  for t in 1:nobs
    pixels[t] = find_pixels_for_line_in_chunklist(clt.chunk_list[t], λmin, λmax, chid).pixels
    #(param_tmp, fit_covar[t], gof[t], fit_converged[t] ) = fit_line(view(clt.chunk_list[t][chid].λ, pixels), view(clt.chunk_list[t][chid].flux,pixels), view(clt.chunk_list[t][chid].var, pixels) )  # some bug in this version of function causes non-convergence
    (param_tmp, fit_covar[t], gof[t], fit_converged[t] ) = fit_line(clt.chunk_list[t][chid].λ, clt.chunk_list[t][chid].flux, clt.chunk_list[t][chid].var, pixels[t] )
    fit_a[t] = param_tmp.a
    fit_λc[t] = param_tmp.λc
    fit_depth[t] = param_tmp.depth
    fit_σ²[t] = param_tmp.σ²
    fit_b[t] = param_tmp.b
  end

  return DataFrame(:fit_a=>fit_a, :fit_λc=>fit_λc, :fit_depth=>fit_depth, :fit_σ²=>fit_σ², :fit_b=>fit_b,
                    :fit_covar=>fit_covar, :χ²_per_dof=>gof, :fit_converged=>fit_converged, :obs_idx =>collect(1:nobs), :chunk_id=>chid, :pixels=>pixels  )
end


function fit_all_lines_in_chunklist_timeseries(clt::AbstractChunkListTimeseries, lines::DataFrame ; plan::LineFinderPlan = LineFinderPlan() )
  @assert size(lines,1) >= 2
  line_idx = 1
  λmin = lines[line_idx,:fit_min_λ]
  λmax = lines[line_idx,:fit_max_λ]
  chid = lines[line_idx,:chunk_id]
  df = fit_line_in_chunklist_timeseries(clt, λmin, λmax, chid)
  #df = fit_line_in_chunklist_timeseries(clt, lines, line_idx)
  df[!,:line_id] .= line_idx
  for line_idx in 2:size(lines,1)
    #df_tmp = fit_line_in_chunklist_timeseries(clt, lines, i)
    λmin = lines[line_idx,:fit_min_λ]
    λmax = lines[line_idx,:fit_max_λ]
    chid = lines[line_idx,:chunk_id]
    df_tmp = fit_line_in_chunklist_timeseries(clt, λmin, λmax, chid)
    df_tmp[!,:line_id] .= line_idx
    append!(df,df_tmp)
  end
  return df
end


function find_worst_telluric_in_each_line!( df::DataFrame, clt::AbstractChunkListTimeseries, data::AbstractArray{AS,1} )  where {AS<:AbstractSpectra}
  min_telluric_model_this_obs = ones(size(df,1))
  min_telluric_model_all_obs  = ones(size(df,1))
  for (i, row) in enumerate(eachrow(df))
    t_idx = row.obs_idx
    ch_idx = row.chunk_id
    pixels = row.pixels
    view_indices = clt.chunk_list[t_idx].data[ch_idx].λ.indices
    order = view_indices[1]
    cols = view_indices[2]
    min_telluric_model_this_obs[i] = minimum(view(data[t_idx].metadata[:tellurics], order, cols)[pixels])
  end
  df.min_telluric_model_this_obs = min_telluric_model_this_obs
  return df
  df_tmp = df |> @groupby(_.line_id) |> @map( { min_telluric_model_all_obs=minimum(_.min_telluric_model_this_obs), line_id=first(_.line_id) } ) |> DataFrame
  df.min_telluric_model_all_obs = df_tmp[telluric_info[!,:line_id],:min_telluric_model_all_obs]
  return df
end


default_Δλ_over_λ_threshold_check_if_line_match = 2.25e-5
function check_if_line_match( list::AbstractVector{T}, λ::Real ; threshold::Real = default_Δλ_over_λ_threshold_check_if_line_match ) where { T<:Real }
  idx = searchsortednearest(list, λ)
  Δ = abs(list[idx]-λ)/λ
  if Δ <= threshold
    return true
  else
    return false
  end
end

function find_which_line_fits_in_line_list(df::DataFrame, line_list::DataFrame ; threshold::Real = default_Δλ_over_λ_threshold_check_if_line_match )
  @assert hasproperty(line_list,:lambda)
  @assert hasproperty(df,:fit_λc)
  @assert issorted(line_list)
  perm = sortperm(df.fit_λc)
  idx = RvSpectML.searchsortednearest(line_list.lambda, df.fit_λc[perm])
  Δ = Array{Float64,1}(undef, size(df,1) )
  Δ[perm] .= (line_list.lambda[idx] .- df.fit_λc[perm])./df.fit_λc[perm]
  line_matches = Δ .<= threshold
  return line_matches
end



#=
function interp_linear(;x1::T1,x2::T1,y1::T2,y2::T2,xpred::T1) where { T1<:Real, T2<:Real }
  ypred = y1+(y2-y1)*((xpred-x1)/(x2-x1))
end
=#

end # module LineFinder
