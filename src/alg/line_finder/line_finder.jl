module LineFinder

using ..RvSpectML
using DataFrames, Query
using Statistics
using LsqFit
#using Optim
#using Plots
#using ThreadedIterables
using ..TemporalGPInterpolation


default_min_deriv2 = 3 # 1.5e5
default_smooth_factor= 4.0
default_min_pixels_in_line = 5
default_use_logλ = true
default_use_logflux = false

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
function find_line_candidates_in_chunk(chunk::AbstractChunkOfSpectrum, deriv2::AbstractVector{T}; plan::LineFinderPlan = LineFinderPlan() ) where { T<:Real }
  idx_d2_gt_0 = findall(d2->isless(zero(d2),d2), deriv2)
  line_candidates = UnitRange[]
  first_idx_in_line = last_idx_in_line = first(idx_d2_gt_0)
    for i in idx_d2_gt_0[2:end]
      if i-last_idx_in_line > 1   # No longer in last set of contiguous indices, process last line candidate
        if last_idx_in_line - first_idx_in_line + 1 < plan.min_pixels_in_line
          # noop continue
        elseif maximum(deriv2[first_idx_in_line:last_idx_in_line]) < plan.min_deriv2 #* norm_in_chunk
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
    elseif maximum(deriv2[first_idx_in_line:last_idx_in_line]) < plan.min_deriv2 #* norm_in_chunk
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


function fit_line(λ::AbstractArray{T1,1}, flux::AbstractArray{T2,1}, var::AbstractArray{T3,1}, idx::UnitRange; show_trace::Bool = false) where { T1<:Real, T2<:Real, T3<:Real }
  @assert length(λ) == length(flux) == length(var)
    λc0 = λ[idx[argmin(flux[idx])]]
  continuum0 = maximum(flux[idx])
  depth0 = 1-minimum(flux[idx])/continuum0
  width0 = 0.05^2   # might need to make more general to deal with different v sin i's
  slope0 = 0.0
  θ0 = pack(a=continuum0,b=slope0,λc=λc0,depth=depth0,σ²=width0)
  lower_bound = pack(a=0.0, b=-1.0, λc=minimum(λ[idx]), depth=0.0, σ²=0.01^2)
  upper_bound = pack(a=2.0, b= 1.0, λc=maximum(λ[idx]), depth=1.0, σ²=0.15^2)
  #upper_bound = pack(a=2.0, b= 1.0, λc=maximum(λ[idx]), depth=1.0, σ²=0.1)
  try
    res = curve_fit(line_model, λ[idx], flux[idx], 1.0 ./(var[idx]), θ0, lower=lower_bound, upper=upper_bound, show_trace=show_trace, autodiff=:forwarddiff)
    covar = estimate_covar(res)
    param = unpack(res.param)
    chi2perdof = sum(abs2.(res.resid).*res.wt) / length(idx)
    #=
    println("idx = ", idx)
    println("lambda = ",λ[idx])
    println("flux = ",flux[idx])
    println("var = ", var[idx])
    println("result = ", res)
    println("param = ", param)
    println("chi2perdof = ", chi2perdof)
    println("res.converged = ",  res.converged)
    @error("good luck!")
    =#
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
Convenience function to find lines in one chunk of spectrum.
# Inputs:
- chunk
# Optional Arguments:
# Return:
- line_fit_list
"""
function find_lines_in_chunk(chunk::AbstractChunkOfSpectrum; plan::LineFinderPlan = LineFinderPlan(),
                              chunk_id::Integer = 0, obs_id::Integer = 0, keep_bad_fits::Bool = false, verbose::Bool = false) where {T <:Real}

  function fit_line(idx::UnitRange)
    #RvSpectML.LineFinder.fit_line(view(chunk.λ,idx), view(flux,idx), view(chunk.var, idx) )   # some bug in this version of function causes non-convergence
    RvSpectML.LineFinder.fit_line( chunk.λ, flux, chunk.var, idx)
  end
  norm_in_chunk = mean(chunk.flux)
  @assert 0.9<=norm_in_chunk<=1.1  # Just making sure we're passing roughly normalized chunks
  # Fit flux with extra smoothing for measuring derivatives well
  #(flux_smooth, deriv_smooth, deriv2_smooth) = predict_mean_and_derivs(chunk.λ, chunk.flux./norm_in_chunk, chunk.λ, sigmasq_obs=chunk.var./norm_in_chunk^2,
  (flux_smooth, deriv_smooth, deriv2_smooth) = predict_mean_and_derivs(chunk.λ, chunk.flux, chunk.λ, sigmasq_obs=chunk.var,
                                  smooth_factor=plan.smooth_factor, use_logx=plan.use_logλ, use_logy=plan.use_logflux )

  # Find pixel ranges that might be lines
  line_candidates = find_line_candidates_in_chunk(chunk, deriv2_smooth, plan=plan)
  if verbose    println("# Found ", size(line_candidates,1), " line candidates.")   end
  if chunk_id != 0    line_candidates[!,:chunk_id] .= chunk_id  end
  if obs_id != 0      line_candidates[!,:obs_id] .= obs_id      end
  #line_candidates[!,:pixels] = line_candidates[!,:pixels]
  line_candidates[!,:fit_min_λ] = chunk.λ[first.(line_candidates[!,:pixels])]
  line_candidates[!,:fit_max_λ] = chunk.λ[ last.(line_candidates[!,:pixels])]


  # Refit flux without extra smoothing, so depths and widths won't be biased by extra smoothing
  #(flux, deriv, deriv2) = predict_mean_and_derivs(chunk.λ, chunk.flux./norm_in_chunk, chunk.λ, sigmasq_obs=chunk.var./norm_in_chunk^2,
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


  if ! keep_bad_fits
    # Avoid "lines" at edge of chunk, either due to contamination or GP fit using zero mean
    line_candidates = line_candidates |> @filter(minimum(_.pixels)>1) |> @filter(maximum(_.pixels)<length(flux)) |> #DataFrame
    #line_candidates = line_candidates |>
        @filter(_.fit_converged) |>
        #@filter(_.fit_σ²<0.018) |>
        @filter(-0.5 <_.fit_b <0.5) |> @filter(_.fit_a < 1.8) |>
        DataFrame
  end

  return line_candidates
end

"""  find_lines_in_chunklist ( chunklist, line_list )
Convenience function to find lines in each chunk of a spectrum.
# Inputs:
- chunklist
# Optional Arguments:
# Return:
- line_list
"""
function find_lines_in_chunklist(chunk_list::AbstractChunkList ; obs_id::Integer = 0, plan::LineFinderPlan = LineFinderPlan() )
  mapreduce(i->find_lines_in_chunk(chunk_list[i], plan=plan, chunk_id=i, obs_id=obs_id), append!, 1:length(chunk_list) )
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
    #@threaded
    map(cl->find_lines_in_chunklist(clt.chunk_list[cl], obsid=cl, plan=plan),1:num_obs(clt))
end


function find_pixels_for_all_lines_in_chunklist( chunk_list::AbstractChunkList, lines::DataFrame)
  @assert hasproperty(lines,:fit_min_λ)
  @assert hasproperty(lines,:fit_max_λ)
  @assert hasproperty(lines,:chunk_id)
  map(l->find_pixels_for_line_in_chunklist(chunk_list[1], lines[l,:fit_min_λ],  lines[l,:fit_max_λ], lines[l,:chunk_id] ), 1:size(lines,1) )
end

global fit_line_in_chunklist_timeseries_count_msgs = 0
function fit_line_in_chunklist_timeseries(clt::AbstractChunkListTimeseries, λmin::Real, λmax::Real, chid::Integer; show_trace::Bool = false)
  df = DataFrame()
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
  global count_msgs
  for t in 1:nobs
    pixels[t] = RvSpectML.find_pixels_for_line_in_chunklist(clt.chunk_list[t], λmin, λmax, chid).pixels
    mean_flux = mean(clt.chunk_list[t][chid].flux[pixels[t]])
    flux = clt.chunk_list[t][chid].flux ./ mean_flux
    var = clt.chunk_list[t][chid].var ./ mean_flux^2
    (param_tmp, fit_covar[t], gof[t], fit_converged[t] ) = fit_line(clt.chunk_list[t][chid].λ, flux, var, pixels[t], show_trace=show_trace )
    #=
    if fit_line_in_chunklist_timeseries_count_msgs<5
      (param_tmp, fit_covar[t], gof[t], fit_converged[t] ) = fit_line(clt.chunk_list[t][chid].λ, flux, var, pixels[t], show_trace=true )
        println("lambda = ",extrema(clt.chunk_list[t][chid].λ[pixels[t]]))
      println("flux = ",extrema(clt.chunk_list[t][chid].flux[pixels[t]]))
      println("var = ", extrema(clt.chunk_list[t][chid].var[pixels[t]]))
      println("flux = ",extrema(flux))
      println("var = ", extrema(var))
      println("pixels = ",extrema(pixels[t]))
      println("param = ", param_tmp)
      println("gof = ", gof[t])
      println("fit_converged = ", fit_converged[t])
      fit_line_in_chunklist_timeseries_count_msgs += 1
    else
      (param_tmp, fit_covar[t], gof[t], fit_converged[t] ) = fit_line(clt.chunk_list[t][chid].λ, flux, var, pixels[t] )
    end
    =#
    fit_a[t] = param_tmp.a
    fit_λc[t] = param_tmp.λc
    fit_depth[t] = param_tmp.depth
    fit_σ²[t] = param_tmp.σ²
    fit_b[t] = param_tmp.b
  end

  return DataFrame(:fit_a=>fit_a, :fit_λc=>fit_λc, :fit_depth=>fit_depth, :fit_σ²=>fit_σ², :fit_b=>fit_b,
                    :fit_covar=>fit_covar, :χ²_per_dof=>gof, :fit_converged=>fit_converged, :obs_idx =>collect(1:nobs), :chunk_id=>chid, :pixels=>pixels  )
end


function fit_all_lines_in_chunklist_timeseries(clt::AbstractChunkListTimeseries, lines::DataFrame ; plan::LineFinderPlan = LineFinderPlan(), show_trace::Bool = false )
  @assert size(lines,1) >= 2
  line_idx::Int64 = 1
  λmin = lines[line_idx,:fit_min_λ]
  λmax = lines[line_idx,:fit_max_λ]
  chid = lines[line_idx,:chunk_id]
  df = fit_line_in_chunklist_timeseries(clt, λmin, λmax, chid, show_trace=show_trace)
  df[!,:line_id] .= line_idx
  for line_idx in 2:size(lines,1)
    λmin = lines[line_idx,:fit_min_λ]
    λmax = lines[line_idx,:fit_max_λ]
    chid = lines[line_idx,:chunk_id]
    df_tmp = fit_line_in_chunklist_timeseries(clt, λmin, λmax, chid)
    df_tmp[!,:line_id] .= line_idx
    append!(df,df_tmp)
  end
  return df
end


# Patterned after find_worst_telluric_in_each_chunk in src/instruments/tellurics.jl
function find_worst_telluric_in_each_line_fit!( df::DataFrame, clt::AbstractChunkListTimeseries, data::AbstractArray{AS,1} )  where {AS<:AbstractSpectra}
  @assert hasproperty(df,:line_id)
  min_telluric_model_this_obs = ones(size(df,1))
  min_telluric_model_all_obs  = ones(size(df,1))
  for (i, row) in enumerate(eachrow(df))
    t_idx = row.obs_idx
    ch_idx = row.chunk_id
    pixels = row.pixels
    view_indices = clt.chunk_list[t_idx].data[ch_idx].λ.indices
    cols = view_indices[1][pixels]
    order = view_indices[2]
    min_telluric_model_this_obs[i] = minimum(view(data[t_idx].metadata[:tellurics], cols, order))
  end
  df.min_telluric_model_this_obs = min_telluric_model_this_obs
  return df
  df_tmp = df |> @groupby(_.line_id) |> @map( { min_telluric_model_all_obs=minimum(_.min_telluric_model_this_obs), line_id=first(_.line_id) } ) |> DataFrame
  df.min_telluric_model_all_obs = df_tmp[telluric_info[!,:line_id],:min_telluric_model_all_obs]
  return df
end

end # module LineFinder
