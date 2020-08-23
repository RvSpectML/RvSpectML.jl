using Pkg
Pkg.activate(".")
#Pkg.add("")
# Pkg.add(["CSV","FITSIO","Distributions","MultivariateStats","Plots","DataFrames","Query","PDMats"])
# Pkg.instantiate()

using DataFrames, Query
using CSV, FITSIO
using Interpolations
using Statistics
using MultivariateStats
#=
using PDMats
using Plots
using Optim
=#

using Stheno # , TemporalGPs

abstract type AbstractSpectra end
abstract type AbstractSpectra1D <: AbstractSpectra end
abstract type AbstractSpectra2D <: AbstractSpectra end

struct Spectra1DBasic{T1<:Real,T2<:Real,T3<:Real,#=T4<:Real,=# AA1<:AbstractArray{T1,1},AA2<:AbstractArray{T2,1},AA3<:AbstractArray{T3,1} } <: AbstractSpectra1D
    λ::AA1
    flux::AA2
    var::AA3
    #doppler_factor::T4
    metadata::Dict{Symbol,Any}
end

struct Spectra2DBasic{T1<:Real,T2<:Real,T3<:Real,#=T4<:Real,=# AA1<:AbstractArray{T1,2},AA2<:AbstractArray{T2,2},AA3<:AbstractArray{T3,2}} <: AbstractSpectra2D
    λ::AA1
    flux::AA2
    var::AA3
    #doppler_factor::T4
    metadata::Dict{Symbol,Any}
end


function Spectra1DBasic(λ::A1, flux::A2, var::A3; #= doppler_factor::T4=one(eltype(λ)), =# metadata::Dict{Symbol,Any} = Dict{Symbol,Any}() ) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,1}, A3<:AbstractArray{T3,1}} #, T4<:Real }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    min_pixels_in_spectra = 1
    max_pixels_in_spectra = 9128
    @assert min_pixels_in_spectra <= length(λ) <= max_pixels_in_spectra
    Spectra1DBasic{eltype(λ),eltype(flux),eltype(var),#=typeof(doppler_factor),=# typeof(λ),typeof(flux),typeof(var)}(λ,flux,var,#=doppler_factor,=# metadata)
end


function Spectra2DBasic(λ::A1, flux::A2, var::A3; #= doppler_factor::T4=one(eltype(λ)), =# metadata::Dict{Symbol,Any} = Dict{Symbol,Any}() ) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}} #, T4<:Real }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    min_pixels_in_spectra = 1
    max_pixels_in_spectra = 9128*128
    min_orders_in_spectra = 1
    max_orders_in_spectra = 128
    @assert min_pixels_in_spectra <= size(λ,1) <= max_pixels_in_spectra
    @assert min_orders_in_spectra <= size(λ,2) <= max_orders_in_spectra
    Spectra2DBasic{eltype(λ),eltype(flux),eltype(var),#=typeof(doppler_factor),=# typeof(λ),typeof(flux),typeof(var)}(λ,flux,var,#=doppler_factor,=# metadata)
end


abstract type AbstractChuckOfSpectra end

struct ChunkOfSpectra{T1<:Real,T2<:Real,T3<:Real,AA1<:AbstractArray{T1,1},AA2<:AbstractArray{T2,1},AA3<:AbstractArray{T3,1}} <: AbstractChuckOfSpectra
    λ::AA1
    flux::AA2
    var::AA3
end

function ChunkOfSpectra{T1,T2,T3}(λ::A1, flux::A2, var::A3) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,1}, A2<:AbstractArray{T2,1}, A3<:AbstractArray{T3,1} }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    min_pixels_in_chuck = 4
    max_pixels_in_chuck = 9128
    @assert min_pixels_in_chuck <= length(λ) <= max_pixels_in_chuck
    ChunkOfSpectra{eltype(λ),eltype(flux),eltype(var),typeof(λ),typeof(flux),typeof(var)}(λ,flux,var)
end

function ChunkOfSpectra(λ::A1, flux::A2, var::A3, order::Integer, pixels::AUR) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, AUR<:AbstractUnitRange }
    @assert size(λ) == size(flux)
    @assert size(λ) == size(var)
    @assert 1 <= order <= size(λ,2)
    @assert 1 <= first(pixels) < last(pixels) <= size(λ,1)
    ChunkOfSpectra{T1,T2,T3}(view(λ,pixels,order),view(flux,pixels,order),view(var,pixels,order))
end


function ChunkOfSpectra(λ::A1, flux::A2, var::A3, loc::NamedTuple{(:pixels, :order),Tuple{AUR,I1}}) where {  T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, AUR<:AbstractUnitRange, I1<:Integer }
    ChunkOfSpectra(λ,flux,var,loc.order,loc.pixels)
end

function ChunkOfSpectra(spectra::AS, order::Integer, pixels::AUR) where { AS<:AbstractSpectra, AUR<:AbstractUnitRange }
    ChunkOfSpectra(spectra.λ,spectra.flux,spectra.var,order,pixels)
end

function ChunkOfSpectra(spectra::AS, loc::NamedTuple{(:pixels, :order),Tuple{AUR,I1}}) where {  AS<:AbstractSpectra, AUR<:AbstractUnitRange, I1<:Integer }
    ChunkOfSpectra(spectra.λ,spectra.flux,spectra.var,loc.order,loc.pixels)
end

abstract type AbstractChunckList end
mutable struct ChunckList{CT<:AbstractChuckOfSpectra, AT<:AbstractArray{CT,1} } <: AbstractChunckList
      data::AT
end

function ChunckList(in::AT) where {CT<:AbstractChuckOfSpectra, AT<:AbstractArray{CT,1} }
    ChunckList{CT,AT}(in)
end
    

abstract type AbstractChunckListTimeseries end
mutable struct ChunckListTimeseries{CLT<:AbstractChunckList, ACLT<:AbstractArray{CLT,1}, TT<:Real, AT<:AbstractArray{TT,1} } <: AbstractChunckListTimeseries
    times::AT
    chuck_list::ACLT
end

function ChunckListTimeseries(t::AT, cl::ACLT) where {CLT<:AbstractChunckList, ACLT<:AbstractArray{CLT,1}, TT<:Real, AT<:AbstractArray{TT,1} }
    ChunckListTimeseries{CLT,ACLT,TT,AT}(t,cl)
end

import Base.length
length(cl::CLT) where {CLT<:AbstractChunckList} = length(cl.data)
length(cl::ACLT) where {ACLT<:AbstractChunckListTimeseries} = length(cl.chuck_list)
num_chunks(cl::ACLT) where {ACLT<:AbstractChunckListTimeseries} = length(first(cl.chuck_list))

const speed_of_light_mps = 299792458.0 # TODO: Update value
calc_doppler_factor(rv::Real) = one(rv) + rv/speed_of_light_mps

function read_metadata_from_fits(fn::String; fields::Array{Symbol,1}, fields_str::AbstractArray{AS,1} )  where { AS<:AbstractString }
    @assert length(fields) == length(fields_str)
    f = FITS(fn)
    hdr = read_header(f[1])
    for field in fields_str
        @assert findfirst(isequal(field),keys(hdr)) != nothing
    end
    values = map(s->hdr[s],fields_str)
    df = Dict{Symbol,Any}(zip(fields,values))    
    return df
end

function read_metadata_from_fits(fn::String, fields::Array{Symbol,1})
    fields_str=string.(fields) 
    read_metadata_from_fits(fn,fields=fields,fields_str=fields_str)
end    
    
function read_metadata_from_fits(fn::String, fields_str::AbstractArray{AS,1} )  where { AS<:AbstractString }
    fields = map(f->Symbol(f),fields_str)
    read_metadata_from_fits(fn,fields=fields,fields_str=fields_str)
end


const Δλoλ_fit_line_default = 5*(1.8*1000/speed_of_light_mps)
const Δλoλ_edge_pad_default = 0*(1.8*1000/speed_of_light_mps)

function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_center::Real; Δ::Real = Δλoλ_fit_line_default) where T<:Real
    findfirst(x->x>=line_center*(1-Δ),wavelengths):findlast(x->x<=line_center*(1+Δ),wavelengths)
end

function find_cols_to_fit(wavelengths::AbstractArray{T,1}, line_lo::Real, line_hi::Real; Δ::Real = Δλoλ_edge_pad_default) where T<:Real
    @assert line_lo < line_hi
    findfirst(x->x>=line_lo*(1-Δ),wavelengths):findlast(x->x<=line_hi*(1+Δ),wavelengths)
end

function find_orders_with_line(goal::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i->order_min(i)<=goal<=order_max(i), 1:size(lambda,2) )
end

function find_orders_with_line(goal_lo::Real,goal_hi::Real,lambda::AbstractArray{T,2}) where T<:Real
   order_min(i) = lambda[1,i]
   order_max(i) = lambda[end,i]
   findall(i->order_min(i)<=goal_lo && goal_hi<=order_max(i), 1:size(lambda,2) )
end

function findall_line(goal::Real,lambda::AbstractArray{T1,2},var::AbstractArray{T2,2}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real}
    @assert lambda[1,1] <= goal <= lambda[end,end]
    orders = find_orders_with_line(goal,lambda)
    @assert length(orders) >= 1
    locs = map(o->(pixels=find_cols_to_fit(lambda[:,o],goal,Δ=Δ),order=o), orders)
    locs_good_idx = findall(t->!any(isnan.(var[t[1],t[2]])),locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]    
    end
    return locs
end

function findall_line(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,2},var::AbstractArray{T2,2}; Δ::Real = Δλoλ_edge_pad_default) where {T1<:Real, T2<:Real}
    @assert lambda[1,1] <= goal_lo < goal_hi <= lambda[end,end]
    orders = find_orders_with_line(goal_lo,goal_hi,lambda)
    #if ! (length(orders) >= 1) return end
    @assert length(orders) >= 1
    locs = map(o->(pixels=find_cols_to_fit(lambda[:,o],goal_lo, goal_hi,Δ=Δ),order=o), orders)
    locs_good_idx = findall(t->!any(isnan.(var[t[1],t[2]])),locs)
    if length(locs) != length(locs_good_idx)
        locs = locs[locs_good_idx]    
    end
    return locs
end

function findall_line(goal::Real,spectra::AS; Δ::Real = Δλoλ_fit_line_default) where {AS<:AbstractSpectra} 
    findall_line(goal,spectra.λ,spectra.var, Δ=Δ)
end

function findall_line(goal_lo::Real,goal_hi::Real,spectra::AS; Δ::Real = Δλoλ_edge_pad_default) where {AS<:AbstractSpectra} 
    findall_line(goal_lo,goal_hi,spectra.λ,spectra.var, Δ=Δ)
end

function calc_snr(flux::AbstractArray{T1},var::AbstractArray{T2}) where {T1<:Real, T2<:Real}
    @assert size(flux) == size(var)
    sum(flux./var)/sum(1.0 ./ var)
end

function find_line_best(goal::Real,lambda::AbstractArray{T1,2},flux::AbstractArray{T2,2},var::AbstractArray{T3,2}; Δ::Real = Δλoλ_fit_line_default) where {T1<:Real, T2<:Real, T3<:Real}
    locs = findall_line(goal,lambda,var,Δ=Δ)
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs) 
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
end

function find_line_best(goal_lo::Real,goal_hi::Real, lambda::AbstractArray{T1,2},flux::AbstractArray{T2,2},var::AbstractArray{T3,2}; Δ::Real = Δλoλ_edge_pad_default) where {T1<:Real, T2<:Real, T3<:Real}
    locs = findall_line(goal_lo,goal_hi,lambda,var,Δ=Δ)
    #scores = map( t->sum( flux[t[1],t[2]] ./ var[t[1],t[2]])/sum( 1.0 ./ var[t[1],t[2]]), locs)
    scores = map( t->calc_snr(flux[t[1],t[2]],var[t[1],t[2]]), locs)
    idx_best = findmax(scores)
    locs[idx_best[2]]
end

function find_line_best(goal::Real,spectra::AS; Δ::Real = Δλoλ_fit_line_default) where {AS<:AbstractSpectra} 
    find_line_best(goal,spectra.λ,spectra.flux,spectra.var, Δ=Δ)
end

function find_line_best(goal_lo::Real,goal_hi::Real,spectra::AS; Δ::Real = Δλoλ_edge_pad_default) where {AS<:AbstractSpectra} 
    find_line_best(goal_lo,goal_hi,spectra.λ,spectra.flux,spectra.var, Δ=Δ)
end


function make_chunck_list(spectra::AS, line_list::AA; Δ::Real=Δλoλ_fit_line_default) where { AS<:AbstractSpectra, R<:Real, AA<:AbstractArray{R,1} }
   ChunckList(map(l->ChunkOfSpectra(spectra,find_line_best(l,spectra,Δ=Δ)), line_list) )
end

function make_chunck_list(spectra::AS, line_list::DataFrame; Δ::Real=Δλoλ_edge_pad_default) where { AS<:AbstractSpectra } 
    @assert haskey(line_list,:lambda_lo)
    @assert haskey(line_list,:lambda_hi)
    ChunckList(map(row->ChunkOfSpectra(spectra,find_line_best(row.lambda_lo,row.lambda_hi,spectra,Δ=Δ)), eachrow(line_list) ))
end

neid_orders_all = 1:90
neid_pixels_all = 1:9216
min_usable_pixels_in_order = 128
#neid_pixels_use = fill(451:neid_pixels_all[end],length(neid_orders_all))
#neid_orders_use = map(pr->length(neid_pixels_use[pr])>min_usable_pixels_in_order,neid_orders_all);

neid_metadata_symbols_default = Symbol[:bjd, :target, :ssbz]
neid_metadata_strings_default = ["OBSJD", "SKY-OBJ", "SSBZ000"]

function read_neid_solar_data(fn::String, metadata::Dict{Symbol,Any} = Dict{Symbol,Any}() )
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] == "Solar"
    if isempty(metadata)
        hdr = read_header(f[1])
        metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    end
    λ, flux, var  = read(f["SKYWAVE"]), read(f["Sky Flux"]), read(f["Sky Variance"])
    Spectra2DBasic(λ, flux, var, metadata=metadata)
end

function read_neid_solar_data(dfr::DataFrameRow{DataFrame,DataFrames.Index})
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_neid_solar_data(fn,metadata)
end
    
function read_metadata_neid_solar(fn::String)
    dict = read_metadata_from_fits(fn,fields=neid_metadata_symbols_default,fields_str=neid_metadata_strings_default)
    dict[:Filename] = fn 
    return dict
end

function add_metadata_from_files_with_neid_solar_data!(df::DataFrame, fn::String)
    metadata_keep = read_metadata_neid_solar(fn)
    if metadata_keep[:target] != "Solar" return df end
    push!(df, metadata_keep)
    return df
end

# 2458745.1296134139 to 2458745.4284185418. Out of those, I usually also require the JD to be less than 2458745.283 
data_path = "20190918"
num_spectra_to_use = 20 
bjd_first_good = 2458745.1296134139
bjd_last_good = 2458745.283

dir_filelist = readdir(data_path,join=true)
idx_spectra = map(fn->occursin(r"^neid\w+\.fits$", last(split(fn,'/')) ),dir_filelist)
spectra_filelist = dir_filelist[idx_spectra]

df_files = DataFrame(Filename = String[], target = String[], bjd = Float64[], ssbz=Float64[] )
map(fn->add_metadata_from_files_with_neid_solar_data!(df_files,fn),spectra_filelist)

df_files_use = df_files |>
  @filter(bjd_first_good <= _.bjd < bjd_last_good) |>
  #@take(num_spectra_to_use) |> 
  DataFrame

@time solar_data = map(read_neid_solar_data,eachrow(df_files_use))



if true
    drift_corrections = CSV.read("SolarRV20190918_JD_SciRV_CalRV.txt", header=["bjd", "sci_drift", "cal_drift"]);

    # TODO: WARNING: NEED TO CHECK SIGN OF DRIFT CORRECTION!!!
function apply_drift_correction!(spectra::S,time::Real,drift_cor::DataFrame; rv_field::Symbol=:cal_drift) where {S<:AbstractSpectra}
    @assert any(isequal(:bjd),propertynames(drift_cor))
    @assert any(isequal(rv_field),propertynames(drift_cor))
    #@assert haskey(spectra.metadata,:doppler_factor)
    idx = findmin(abs.(time .- drift_cor.bjd))[2]
    doppler_factor = calc_doppler_factor(drift_cor[rv_field][idx])
    # println("# t= ",time, " drift_cor= ",drift_cor[rv_field][idx], " df= ",doppler_factor)
    spectra.λ .*= doppler_factor
    if haskey(spectra.metadata,:doppler_factor)
        spectra.metadata[:doppler_factor] /= doppler_factor
    else
        spectra.metadata[:doppler_factor] = 1/doppler_factor    
    end
    return spectra
end

@time map(x->apply_drift_correction!(x[1],x[2],drift_corrections), zip(solar_data,df_files_use.bjd));
drift_corrections[:cal_drift].-mean(drift_corrections[:cal_drift])
end


if true
    ssb_corrections = CSV.read("SolarTelescope2019-09-18_inclGravRedshiftAirMassAltitude.csv", header=["bjd","rv_ssb"], select=[1,2], types=types=[Float64,Float64], datarow=2, silencewarnings=true);
@assert any(isequal(:bjd),propertynames(ssb_corrections))
@assert any(isequal(:rv_ssb),propertynames(ssb_corrections))
ssb_interp = LinearInterpolation(ssb_corrections.bjd, ssb_corrections.rv_ssb)

function apply_ssb_correction!(spectra::AS,time::Real,ssb_cor::AI) where {AS<:AbstractSpectra,AI<:AbstractInterpolation}    
    doppler_factor = calc_doppler_factor(ssb_cor(time))
    #println("# t= ",time, " ssb_cor= ",ssb_cor(time), " df= ",doppler_factor)
    spectra.λ .*= doppler_factor  
    if haskey(spectra.metadata,:doppler_factor)
        spectra.metadata[:doppler_factor] /= doppler_factor
    else
        spectra.metadata[:doppler_factor] = doppler_factor    
    end
    return spectra
end
@time map(x->apply_ssb_correction!(x[1],x[2],ssb_interp), zip(solar_data,df_files_use.bjd));
ssb_interp.(df_files_use.bjd)
end

lambda_range_with_data = (min = maximum(d->minimum(d.λ),solar_data), max = minimum(d->maximum(d.λ),solar_data) )

#vald_filename = "VALD_Fe1_DP_rejectTelluricSlope0.0_badLineFilterESPRESSO_overlapcutoff6e-05_depthcutoff0.05_allowBlends0_wavesoriginal_depthsoriginal_nbin1depth0.mas"
vald_filename = "VALD_Fe1_DP_rejectTelluricSlope0.0_badLineFilterESPRESSO-strict-NEID-BIS_overlapcutoff6e-05_depthcutoff0.05_allowBlends0_wavesReiners_depthssolar_nbin1depth0.mas"
vald_df = CSV.read(vald_filename,header=["lambda_lo","lambda_hi","depth_vald"])
line_list_df = vald_df |>
  @filter(lambda_range_with_data.min <= _.lambda_lo ) |>
  @filter( _.lambda_hi < lambda_range_with_data.max) |>
#  @filter( _.lambda_lo >6157 || _.lambda_hi < 6155  ) |>   # "Line" w/ large variability
  DataFrame;
size(line_list_df)

function predict_line_width(Teff::Real; v_rot::Real=0.0)
    @assert 3000 < Teff < 10000 # K 
    @assert 0 <= v_rot <=100 # km/s
    line_width_thermal = 13*sqrt(Teff/1e4) # km/s
    line_width = sqrt(v_rot^2+line_width_thermal^2) # km/s
end

#@time chuck_list = make_chuck_list(solar_data[1],line_list_df)
# Adjust width of chunks
line_list_df[:lambda_mid] = sqrt.(line_list_df.lambda_lo.*line_list_df.lambda_hi)
chunk_size_factor = 3       # TODO: Figure out what value to use
max_vald_line_offset = 0.0       # km/s
line_width = predict_line_width(5780,v_rot=1.8) # # km/s
Δλoλ_fit_line = (max_vald_line_offset+chunk_size_factor*line_width)*1000/speed_of_light_mps
println("# Δλ/λ = ",Δλoλ_fit_line)
line_list_df.lambda_hi .= line_list_df.lambda_mid*(1 + Δλoλ_fit_line)
line_list_df.lambda_lo .= line_list_df.lambda_mid/(1 + Δλoλ_fit_line)

chunk_list_df = line_list_df

# Check if chunks overlap
if any(line_list_df.lambda_hi[1:end-1] .>= line_list_df.lambda_lo[2:end])    
    idx_overlap = line_list_df.lambda_hi[1:end-1] .>= line_list_df.lambda_lo[2:end]
    println("# Overlapping chunks: ",length(findall(idx_overlap)))
end
minimum(line_list_df.lambda_lo), maximum(line_list_df.lambda_hi)

function merge_lines(line_list::DataFrame)
    chunk_list_df = DataFrame(:lambda_lo=>Float64[],:lambda_hi=>Float64[],
                                :line_λs=>Array{Float64,1}[],:line_depths=>Array{Float64,1}[])
    num_lines = size(line_list,1) 
    @assert num_lines >= 2
    lambda_lo_last = line_list[1,:lambda_lo]
    lambda_hi_last = line_list[1,:lambda_hi]
    line_λs = [line_list[1,:lambda_mid]]
    line_depths = [line_list[1,:depth_vald]]
    for i in 2:num_lines
        lambda_lo = line_list[i,:lambda_lo]
        lambda_hi = line_list[i,:lambda_hi]
        if lambda_lo>lambda_hi_last
            push!(chunk_list_df, (lambda_lo_last, lambda_hi_last, line_λs, line_depths))
            (lambda_lo_last, lambda_hi_last) = (lambda_lo, lambda_hi)
            line_λs = [line_list[i,:lambda_mid]]
            line_depths = [line_list[i,:depth_vald]]
        else
            lambda_hi_last = lambda_hi
            push!(line_λs,line_list[i,:lambda_mid])
            push!(line_depths,line_list[i,:depth_vald])
        end
    end
    if chunk_list_df[end,:lambda_hi] != lambda_hi_last
        #push!(chunk_list_df, (lambda_lo_last, lambda_hi_last))
        push!(chunk_list_df, (lambda_lo_last, lambda_hi_last, line_λs, line_depths))
    end
    return chunk_list_df
end

chunk_list_df = merge_lines(line_list_df)
size(chunk_list_df)
#chunk_list_df

solar_data[1].λ[:,1:90]

#@time time_series_of_chunk_lists = map(spec->make_chunck_list(spec,line_list_df),solar_data)
@time time_series_of_chunk_lists = map(spec->make_chunck_list(spec,chunk_list_df),solar_data)
chunk_list_timeseries = ChunckListTimeseries(df_files_use.bjd,time_series_of_chunk_lists)

chunk_list_timeseries.chuck_list[5].data[10].flux

#chunk_list_timeseries.chuck_list[2].data[10].flux
#solar_data[5].flux ./= 100

# If want to keep orders as one big chunk
min_col_neid_default = 451
max_col_neid_default = 9216 #- (min_col_neid_default-1)
neid_orders_to_use = 1:90
pixels_to_use_neid = fill(min_col_neid_default:max_col_neid_default,length(neid_orders_to_use))
if  73<maximum(neid_orders_to_use)   pixels_to_use_neid[73]=1940:max_col_neid_default    end

if true # If want to break orders into 1024 pixel chunks
neid_orders_to_use = 1:60
pixels_to_use_neid = repeat(map(i->560+i*1024:560+(i+1)*1024,0:7),length(neid_orders_to_use))
neid_orders_to_use = mapreduce(i->repeat([i],8),vcat,neid_orders_to_use)
end

function make_orders_into_chunks(spectra::AS; orders_to_use=1:size(spectra.flux,2), 
        min_col::Integer=min_col_neid_default, max_col::Integer=max_col_neid_default , 
        pixels_to_use=fill(min_col:max_col,length(orders_to_use)) ) where {AS<:AbstractSpectra, }
    ChunckList(map(order->ChunkOfSpectra(spectra,(pixels=pixels_to_use[order],order=order)), orders_to_use ))
end

@time time_series_of_order_lists = map( spec->make_orders_into_chunks(spec, 
                orders_to_use=neid_orders_to_use, pixels_to_use=pixels_to_use_neid) ,solar_data)
order_list_timeseries = ChunckListTimeseries(df_files_use.bjd,time_series_of_order_lists)

# Check that no NaN's included
map(t->any(map(o->any(isnan.(solar_data[t].λ[pixels_to_use_neid[o],o])),1:60)),length(solar_data))

bad_pixels = map(x->(x[1],x[2]),findall(isnan.(solar_data[5].flux)))
not_in_bad_collum_idx = findall(x->(x[1]<439 || x[1]>450),bad_pixels)
bad_pixels[not_in_bad_collum_idx]

function filter_bad_chunks(chunk_list_timeseries::ACLT, line_list::DataFrame; verbose::Union{Int,Bool} = false) where { ACLT<:AbstractChunckListTimeseries }
    @assert(haskey(line_list,:lambda_lo))
    @assert(haskey(line_list,:lambda_hi))
    @assert(length(chunk_list_timeseries)>=1)
    idx_keep = trues(num_chunks(chunk_list_timeseries))
    for t in 1:length(chunk_list_timeseries)
        idx_bad_λ = findall(c->any(isnan.(chunk_list_timeseries.chuck_list[t].data[c].λ)),1:num_chunks(chunk_list_timeseries))
        idx_bad_flux = findall(c->any(isnan.(chunk_list_timeseries.chuck_list[t].data[c].flux)),1:num_chunks(chunk_list_timeseries))
        idx_bad_var = findall(c->any(isnan.(chunk_list_timeseries.chuck_list[t].data[c].var)),1:num_chunks(chunk_list_timeseries))
        idx_keep[idx_bad_λ] .= false
        idx_keep[idx_bad_flux] .= false
        idx_keep[idx_bad_var] .= false
        if verbose && (length(idx_bad_λ)+length(idx_bad_flux)+length(idx_bad_var) > 0)
               println("# Removing chunks", vcat(idx_bad_λ,idx_bad_flux,idx_bad_var), " at time ", t, " due to NaNs (",length(idx_bad_λ),",",length(idx_bad_flux),",",length(idx_bad_var),").")
        end
        
    end
    chunks_to_remove = findall(.!idx_keep)
    if length(chunks_to_remove) == 0
        println("# No lines to remove.")
        return (chunk_timeseries=chunk_list_timeseries, line_list=line_list)
    else
        println("# Removing ", length(chunks_to_remove), " chunks due to NaNs.")
        map(c->println("# ",c,": ",line_list.lambda_lo[c]," - ",line_list.lambda_hi[c]),chunks_to_remove)
        new_line_list = line_list[findall(idx_keep),:]
        new_chunk_list_timeseries = [ChunckList(chunk_list_timeseries.chuck_list[t].data[idx_keep]) for t in 1:length(chunk_list_timeseries) ]
        return (chunk_timeseries=ChunckListTimeseries(chunk_list_timeseries.times,new_chunk_list_timeseries), line_list=new_line_list)
    end
end

println(size(chunk_list_df), " vs ", num_chunks(chunk_list_timeseries) )
(chunk_list_timeseries, chunk_list_df) = filter_bad_chunks(chunk_list_timeseries,chunk_list_df)
println(size(chunk_list_df), " vs ", num_chunks(chunk_list_timeseries) )

chunk_list_df[1,:]



function calc_normalization(chunk_list::ACL) where { ACL<:AbstractChunckList}
    total_flux = sum(sum(Float64.(chunk_list.data[c].flux))
                        for c in 1:length(chunk_list) )
    num_pixels = sum( length(chunk_list.data[c].flux) for c in 1:length(chunk_list) )
    scale_fac = num_pixels / total_flux
end

function normalize_spectrum!(spectrum::ST, scale_fac::Real) where { ST<:AbstractSpectra }
    @assert 0 < scale_fac < Inf
    @assert !isnan(scale_fac^2)
    spectrum.flux .*= scale_fac
    spectrum.var .*= scale_fac^2
    return spectrum
end

function normalize_spectra!(timeseries::ACLT, spectra::AS) where { ACLT<:AbstractChunckListTimeseries, ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    @assert length(timeseries) == length(spectra)
    for t in 1:length(timeseries)  
        scale_fac = calc_normalization(timeseries.chuck_list[t])
        println("# t= ",t, " scale_fac= ", scale_fac)
        normalize_spectrum!(spectra[t], scale_fac)
    end
    return timeseries
end

#normalize_spectra!(chunk_list_timeseries,solar_data)
normalize_spectra!(order_list_timeseries,solar_data);

solar_data[5].flux[4000:6000]
chunk_list_timeseries.chuck_list[5].data[100].flux

using Plots

#chunk_list_timeseries.chuck_list[1].data[1].λ

order_idx = 20:22
xmin = Inf # minimum(order_list_df.lambda_lo[chunk_idx])
xmax = 0   # maximum(order_list_df.lambda_hi[chunk_idx])
plt = plot(legend=:none)
for c in order_idx
    t = 1
    plot!(plt,order_list_timeseries.chuck_list[t].data[c].λ ,order_list_timeseries.chuck_list[t].data[c].flux)
    xmin = min(xmin,minimum(order_list_timeseries.chuck_list[t].data[c].λ))
    xmax = max(xmax,maximum(order_list_timeseries.chuck_list[t].data[c].λ))
end
#xlims!(xmin,xmax)
display(plt)


chunk_idx = 12:20
xmin = minimum(chunk_list_df.lambda_lo[chunk_idx])
xmax = maximum(chunk_list_df.lambda_hi[chunk_idx])
#plt = plot(legend=:none)
#xlims!(xmin,xmax)
for c in chunk_idx
    t = 1
    #if(sum(chunk_list_df.line_depths[c])<0.25) continue end
    λ_mid = 0# sqrt(chunk_list_df.lambda_hi[c]*chunk_list_df.lambda_lo[c])
    println("c= ",c , " λs= ",chunk_list_df.line_λs[c]," depths= ",chunk_list_df.line_depths[c])
    #println("  λlo= ",chunk_list_df.lambda_lo[c]," λhi= ",chunk_list_df.lambda_hi[c], " Δλ= ",chunk_list_df.lambda_hi[c]-chunk_list_df.lambda_lo[c])
    plot!(plt,chunk_list_timeseries.chuck_list[t].data[c].λ.-λ_mid,chunk_list_timeseries.chuck_list[t].data[c].flux)
end
#plot!(plt,solar_data[1].λ,solar_data[1].flux)
#xlims!(4560,4565)
display(plt)




function make_grid_for_chunck(timeseries::ACLT, c::Integer) where { ACLT<:AbstractChunckListTimeseries }
    num_obs = length(timeseries.chuck_list)
    λ_min = maximum(minimum(timeseries.chuck_list[t].data[c].λ) for t in 1:num_obs)
    λ_max = minimum(maximum(timeseries.chuck_list[t].data[c].λ) for t in 1:num_obs)
    Δλ_grid_obs = median((timeseries.chuck_list[t].data[c].λ[end]-
            timeseries.chuck_list[t].data[c].λ[1])/  
              (length(timeseries.chuck_list[t].data[c].λ)-1) for t in 1:num_obs)
    num_pixels_obs = mean(length(timeseries.chuck_list[t].data[c].λ) for t in 1:num_obs)
    oversample_factor = 1
    num_pixels_gen = (num_pixels_obs-1) * oversample_factor + 1
    Δλ_grid_gen = (λ_max-λ_min)/ (num_pixels_gen-1)
    range(λ_min,stop=λ_max,step=Δλ_grid_gen)  
end

chunk_grids = map(c->make_grid_for_chunck(chunk_list_timeseries,c), 1:num_chunks(chunk_list_timeseries) )

# Doesn't work because looking up lambda_min and lambda_max, but haven't been updated since splitting order into segments
#chunk_grids = map(c->make_grid_for_chunck(order_list_timeseries,c), 1:num_chunks(order_list_timeseries) )

function make_interpolator_linear_flux(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectra}
    LinearInterpolation(spectra.λ, spectra.flux)
end

function make_interpolator_linear_var(spectra::Union{AS,AC}) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectra}
    LinearInterpolation(spectra.λ, spectra.var)
end

using Stheno

function make_interpolator_gp(spectra::Union{AS,AC}; length_scale::Real = 0.1, σ_scale::Real = 1.0) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectra}
    xobs = spectra.λ
    yobs = copy(spectra.flux)
    obs_var = spectra.var
    med_y = mean(yobs)
    
    # Choose the length-scale and variance of the process.
    σ² = (σ_scale * med_y)^2
    # Construct a kernel with this variance and length scale.
    k = σ² * stretch(Matern52(), 1 / length_scale)

    # Specify a zero-mean GP with this kernel. Don't worry about the GPC object.
    f = GP(med_y, k, GPC())
    fx = f(xobs, obs_var)
    f_posterior = f | Obs(fx, yobs)
end

function interp_to_grid(spectra::Union{AS,AC}, grid::AR) where { AS<:AbstractSpectra, AC<:AbstractChuckOfSpectra, AR<:AbstractRange}
   #grid = chunk_grids[c]
   #make_interpolator_linear(spectra).(grid) 
   f_posterior = make_interpolator_gp(spectra,length_scale=6e-5*mean(spectra.λ))
   (mean=mean(f_posterior(grid)),  std=std.(marginals(gp_interp(grid))))
end
    

function pack_chunks_into_matrix(timeseries::ACLT, chunk_grids::AR) where { ACLT<:AbstractChunckListTimeseries, RT<:AbstractRange, AR<:AbstractArray{RT,1} }
   num_obs = length(timeseries)
   num_λ = sum(length.(chunk_grids))
   flux_matrix = Array{Float64,2}(undef,num_λ,num_obs)
   var_matrix = Array{Float64,2}(undef,num_λ,num_obs)
   λ_vec = Array{Float64,1}(undef,num_λ)
   chunk_map = Array{UnitRange{Int64}}(undef,length(chunk_grids))
    
   for t in 1:num_obs
        idx_start = 0
        for c in 1:length(chunk_grids)
            if length(chunk_grids[c]) <= 512
                gp_interp = make_interpolator_gp(timeseries.chuck_list[t].data[c],length_scale=1e-4*mean(chunk_grids[c]))
                idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
                flux_matrix[idx,t] .= mean(gp_interp(chunk_grids[c]))
                var_matrix[idx,t] .= var.(marginals(gp_interp(chunk_grids[c])))
                if t == 1 
                    λ_vec[idx] .= chunk_grids[c]
                    chunk_map[c] = idx
                    #= if 8500<idx_start <8700
                        println(c,": ",idx)
                    end
                    =#
                end
                idx_start += length(chunk_grids[c])
            else
                lin_interp_flux = make_interpolator_linear_flux(timeseries.chuck_list[t].data[c])
                lin_interp_var = make_interpolator_linear_var(timeseries.chuck_list[t].data[c])
                idx = (idx_start+1):(idx_start+length(chunk_grids[c]))
                flux_matrix[idx,t] .= lin_interp_flux.(chunk_grids[c])
                var_matrix[idx,t] .= lin_interp_var.(chunk_grids[c])
                if t == 1 
                    λ_vec[idx] .= chunk_grids[c]
                    chunk_map[c] = idx
                    #= if 8500<idx_start <8700
                        println(c,": ",idx)
                    end
                    =#
                end
                idx_start += length(chunk_grids[c])
            end
        end
    end
    
    #return flux_matrix
    return (flux=flux_matrix, var=var_matrix, λ=λ_vec, chunk_map=chunk_map)
end

mean(chunk_grids[1])

(fm, vm, λv, cm) = pack_chunks_into_matrix(chunk_list_timeseries,chunk_grids)

size(fm)
fm_mean = sum(fm./vm,dims=2)./sum(1.0./vm,dims=2)
#fm[cm[1]],λv[cm[1]]

function calc_deriv(flux::AbstractArray{T1,1}, λ::AbstractArray{T2,1}) where { T1<:Real, T2<:Real }
    @assert size(flux) == size(λ) 
    @assert length(flux) >= 3
    dfdlogλ = Array{T1,1}(undef,length(flux))
    dfdlogλ[1] = 0.5*(flux[2]-flux[1])/(λ[2]-λ[1])*(λ[2]+λ[1])
    dfdlogλ[2:end-1] .= 0.5*(flux[3:end].-flux[1:end-2])./(λ[3:end].-λ[1:end-2]).*(λ[3:end].+λ[end-2])
    dfdlogλ[end] = 0.5*(flux[end]-flux[end-1])/(λ[end]-λ[end-1])*(λ[end]+λ[end-1])
    return dfdlogλ
end

function calc_mean_deriv(flux::AbstractArray{T1,2}, var::AbstractArray{T1,2}, λ::AbstractArray{T3,1}, 
        chunk_map::AbstractArray{URT,1}) where
    { T1<:Real, T2<:Real, T3<:Real, URT<:AbstractUnitRange} #, RT<:AbstractRange }
    flux_mean = sum(flux./vm,dims=2)./sum(1.0./vm,dims=2)
    deriv = Array{T1,1}(undef,length(flux_mean))
    map(c->deriv[c] .= calc_deriv(flux_mean[c],λv[c]),chunk_map )
    return deriv
end

#fm
#fm_mean = mean(fm,dims=2)
deriv_old = copy(deriv)
deriv = calc_mean_deriv(fm,vm,λv,cm)
idx_plt = cm[13]
rvs = vec(sum((fm[idx_plt,:].-fm_mean[idx_plt]).*deriv[idx_plt],dims=1)./sum(abs2.(deriv[idx_plt]))).*speed_of_light_mps
#rvs = vec(sum((fm[idx_plt,:].-fm_mean[idx_plt]).*deriv[idx_plt],dims=1)./sum(abs2.(deriv[idx_plt]))).*speed_of_light_mps
rvs .-= mean(rvs)
println(rvs)
plot(λv[idx_plt] ,fm_mean[idx_plt].-1.0, label="mean")
plot!(λv[idx_plt], deriv[idx_plt]./std(deriv[idx_plt]), label="deriv")
plot!(λv[idx_plt], deriv_old[idx_plt]./std(deriv_old[idx_plt]), label="deriv old")

# Compute RVs
zs = vec(sum((fm.-fm_mean).*deriv,dims=1)./sum(abs2.(deriv)))
rvs = (zs.-mean(zs)).*speed_of_light_mps

σ_rvs = sqrt.(vec(sum(vm.*abs.(deriv),dims=1)./sum(abs2.(deriv))))./length(deriv).*speed_of_light_mps

std(rvs[idx_good])

idx_good = ( (chunk_list_timeseries.times.-minimum(chunk_list_timeseries.times))*hours_per_day .>= 1.0 )
                  # .&  (-3 .<= pca_out[1,:] .<= 3 )

hours_per_day = 24
plot((chunk_list_timeseries.times[idx_good].-minimum(chunk_list_timeseries.times))*hours_per_day,rvs[idx_good],
        yerr=σ_rvs[idx_good],xlabel="Time (hr)",ylabel="RV (m/s)", legend=:none)

obs_per_bin = 5
num_obs_binned = floor(Int,length(findall(idx_good))//obs_per_bin)
idx_binned = map(i->1+(i-1)*obs_per_bin:obs_per_bin*i,1:num_obs_binned)
rvs_binned = map(i->mean(rvs[findall(idx_good)][i]),idx_binned)
times_binned = (map(i->mean(chunk_list_timeseries.times[findall(idx_good)][i]),idx_binned).-minimum(chunk_list_timeseries.times))*hours_per_day
scatter!(times_binned,rvs_binned)

rms_rvs_binned = std(rvs_binned)

fm_perp = fm .- zs'.* deriv
M = fit(PCA, fm_perp[:,idx_good]; maxoutdim=12)
pca_out = MultivariateStats.transform(M,fm_perp[:,idx_good])

scatter(1.0.-cumsum(principalvars(M))./tvar(M), xlabel="Number of PCs", ylabel="Frac Variance Unexplained")

idx_plt = cm[13]
plt0 = plot((fm_mean[idx_plt].-1.0)./std(fm_mean[idx_plt]),legend=:none)
plt0 = plot!(deriv[idx_plt]./std(deriv[idx_plt]),legend=:none)
plt0 = plot!(M.proj[idx_plt,1]./std(M.proj[idx_plt,1]),legend=:none)
plt0 = plot!(M.proj[idx_plt,2]./std(M.proj[idx_plt,2]),legend=:none)
#plt1 = plot(M.proj[idx_plt,1],legend=:none)
plt1 = plot(M.proj[idx_plt,1]./std(M.proj[idx_plt,1]),legend=:none)
plt1 = plot!(M.proj[idx_plt,2]./std(M.proj[idx_plt,2]),legend=:none)
plt2 = plot(M.proj[idx_plt,2],legend=:none)
plt3 = plot(M.proj[idx_plt,3],legend=:none)
plt4 = plot(M.proj[idx_plt,4],legend=:none)
plot(plt0,plt1,plt2,plt3,plt4, layout = (5,1) )

plot((chunk_list_timeseries.times[idx_good].-minimum(chunk_list_timeseries.times))*hours_per_day,rvs[idx_good],
        yerr=σ_rvs[idx_good],xlabel="Time (hr)",ylabel="RV (m/s)", legend=:none)
plot!((chunk_list_timeseries.times[idx_good].-minimum(chunk_list_timeseries.times))*hours_per_day,vec(pca_out[1,:]),
        xlabel="Time (hr)",ylabel="PC1", legend=:none)
plot!((chunk_list_timeseries.times[idx_good].-minimum(chunk_list_timeseries.times))*hours_per_day,vec(pca_out[2,:]),
        xlabel="Time (hr)",ylabel="PC2", legend=:none)

plt1 = scatter(rvs[idx_good],vec(pca_out[1,:]),xlabel="RV",ylabel="PC1",legend=:none)
plt2 = scatter(rvs[idx_good],vec(pca_out[2,:]),xlabel="RV",ylabel="PC2",legend=:none)
plt3 = scatter(rvs[idx_good],vec(pca_out[3,:]),xlabel="RV",ylabel="PC3",legend=:none)
plt4 = scatter(vec(pca_out[1,:]),vec(pca_out[2,:]),xlabel="PC1",ylabel="PC2",legend=:none)
plot(plt1,plt2,plt3,plt4,layout=(2,2))

M2 = fit(PCA, pca_out; maxoutdim=1)

pca_out2 = vec(MultivariateStats.transform(M2, pca_out ))


scatter(rvs[idx_good],pca_out2,xlabel="RV",ylabel="PC1",legend=:none)
plot!(rvs[idx_good],pca_out2,xlabel="RV",ylabel="PC1",legend=:none)



if false
   ts = chunk_list_timeseries
   t = 1
   c = 1
   interp = LinearInterpolation(ts.chuck_list[t].data[c].λ, ts.chuck_list[t].data[c].flux)
   
    xobs = (ts.chuck_list[t].data[c].λ)
    yobs = copy(ts.chuck_list[t].data[c].flux)
    obs_var = (ts.chuck_list[t].data[c].var)
    med_y = median(yobs)
    
    # Choose the length-scale and variance of the process.
    l = 0.1
    σ² = 1.0 * med_y^2
    # Construct a kernel with this variance and length scale.
    k = σ² * stretch(Matern52(), 1 / l)

    # Specify a zero-mean GP with this kernel. Don't worry about the GPC object.
    
    f = GP(med_y, k, GPC())
    fx = f(xobs, obs_var)
    f_posterior = f | Obs(fx, yobs)
    #mean(f_posterior(chunk_grids[c]))
    
   grid = chunk_grids[c]
   plt1 = plot(grid,mean(f_posterior(grid)))
   plot!(plt1,grid,interp.(grid))
   #plot!(grid,mean(f_posterior(grid)))
   scatter!(plt1,ts.chuck_list[t].data[c].λ, ts.chuck_list[t].data[c].flux)
   plt2 = scatter(ts.chuck_list[t].data[c].λ, mean(f_posterior(ts.chuck_list[t].data[c].λ)).-ts.chuck_list[t].data[c].flux)
   #scatter!(plt2,ts.chuck_list[t].data[c].λ, interp.(ts.chuck_list[t].data[c].λ).-ts.chuck_list[t].data[c].flux)
   plot(plt1,plt2,layout=(2,1)) #display(plot)
    
end

@inline function allequal(x::AbstractArray{T,1}) where {T<:Real}
    length(x) < 2 && return true
    e1 = x[1]
    i = 2
    @inbounds for i=2:length(x)
        x[i] == e1 || return false
    end
    return true
end

if false
    
    map(c->allequal(map(t->length(chunk_list_timeseries.chuck_list[t].data[c].λ),1:length(chunk_list_timeseries))),1:10)

    c=2
    println(map(t->length(chunk_list_timeseries.chuck_list[t].data[c].λ),1:length(chunk_list_timeseries)))
end
    

if false
    θ_init = [-0.02, 0.1]
xobs = sol_wave[idx_cols,order] .- line_center 
yobs = convert(Array{Float64,1},sol_flux[idx_cols,order] )
sigmaobs = convert(Array{Float64,1},sol_var[idx_cols,order] )
(ypred, coeff) = fit_and_predict(θ_init, x=xobs, yobs=yobs, sigmaobs=sigmaobs, degree_poly=1)
println("RMS = ",calc_rms_error(ypred,yobs), " χ^2 = ", calc_chi_sq(ypred,yobs,sigmaobs), " dof = ", length(xobs))
println("coeff = ", coeff)
scatter(xobs,yobs,legend=:none)
plot!(xobs,ypred,legend=:none)
end
