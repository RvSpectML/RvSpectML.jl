"""
Convenience code to store info about what should be (re)calculated, saved to disk and/or plotted,
and what can be skipped.

Author: Eric Ford
Created: Sept 2020
"""

"""
The PipelinePlan stores what work needs to be done, what plots should be made, and what data/plots should be saved to disk.
Users will querty it via  need_to(plan, symbol), make_plot(plan, symbol), save_plot(plan, symbol), and save_data(plan, symbol).
Users can update the plan with make_plot!(plan, symbol),  dont_make_plot!(plan, symbol), etc.

In principle, we could cache data here, too.  But I'm not sure if we'll actulaly use that.  So consider the cache as experimental.
"""
module Pipeline

export PipelinePlan
export make_plot, save_plot, save_data, need_to, has_cache

export make_plot!, dont_make_plot!, make_all_plots!,  make_no_plots!
export save_plot!, dont_save_plot!, save_all_plots!,  save_no_plots!
export save_data!, dont_save_data!, save_all_data!,   save_no_data!
export need_to!,   dont_need_to!,   reset_all_needs!, reset_no_needs!
export has_cache,  read_cache, set_cache!, reset_all_cache!

# Default steps for the PipelinePlan to track
default_pipeline_input_collect(keys)  = [ :read_spectra, :read_line_list ]
default_pipeline_work_collect(keys)   = [ :extract_orders, :clean_line_list_tellurics, :fit_lines, :clean_line_list_blends  ]
default_pipeline_output_collect(keys) = [ :ccf_pass0, :ccf_total, :rvs_ccf_total, :ccf_orders, :rvs_ccf_orders, :scalpels, :template, :dcpca ]
default_pipeline_need_to_collect(keys) = vcat( default_pipeline_input_collect(keys),  default_pipeline_work_collect(keys), default_pipeline_output_collect(keys) )

DictSB = Dict{Symbol,Bool}
CacheT = Dict{Symbol,Any}

struct PipelinePlan
    make_plot::DictSB
    save_plot::DictSB
    save_data::DictSB
    need_to::DictSB
    cache::CacheT
end


make_all_plots_dict = DictSB( zip(default_pipeline_output_collect(keys),  trues(length(default_pipeline_output_collect(keys)))  ) )
save_all_plots_dict = DictSB( zip(default_pipeline_output_collect(keys),  trues(length(default_pipeline_output_collect(keys)))  ) )
save_all_files_dict = DictSB( zip(default_pipeline_output_collect(keys),  trues(length(default_pipeline_output_collect(keys)))  ) )
need_to_do_everything_dict = DictSB( zip(default_pipeline_need_to_collect(keys), trues(length(default_pipeline_need_to_collect(keys)))  ) )

make_no_plots_dict = DictSB( zip(default_pipeline_output_collect(keys),  falses(length(default_pipeline_output_collect(keys)))  ) )
save_no_plots_dict = DictSB( zip(default_pipeline_output_collect(keys),  falses(length(default_pipeline_output_collect(keys)))  ) )
save_no_files_dict = DictSB( zip(default_pipeline_output_collect(keys),  falses(length(default_pipeline_output_collect(keys)))  ) )
need_to_do_nothing_dict = DictSB( zip(default_pipeline_need_to_collect(keys), falses(length(default_pipeline_need_to_collect(keys)))  ) )

default_pipeline_make_plot = make_all_plots_dict
default_pipeline_save_plot = save_no_plots_dict
default_pipeline_save_file = save_no_files_dict
default_pipeline_need_to   = need_to_do_everything_dict
default_pipeline_cached_data = CacheT()

function PipelinePlan(; make_plot::DictSB = default_pipeline_make_plot,
                        save_plot::DictSB = default_pipeline_save_plot,
                        save_file::DictSB = default_pipeline_save_file,
                        need_to::DictSB   = default_pipeline_need_to,
                        cached_data::Dict{Symbol,Any} = default_pipeline_cached_data )
    PipelinePlan(make_plot, save_plot, save_file, need_to, cached_data )
end

make_plot(p::PipelinePlan, s::Symbol) = haskey(p.make_plot,s) && p.make_plot[s]
make_plot!(p::PipelinePlan, s::Symbol) = p.make_plot[s] = true
dont_make_plot!(p::PipelinePlan, s::Symbol) = p.make_plot[s] = false
make_all_plots!(p::PipelinePlan, s::Symbol) = map(k -> p.make_plot[k]=true, collect(keys(p.make_plot)) )
make_no_plots!(p::PipelinePlan, s::Symbol) = map(k -> p.make_plot[k]=false, collect(keys(p.make_plot)) )

save_plot(p::PipelinePlan, s::Symbol) = haskey(p.save_plot,s) && p.save_plot[s]
save_plot!(p::PipelinePlan, s::Symbol) = p.save_plot[s] = true
dont_save_plot!(p::PipelinePlan, s::Symbol) = p.save_plot[s] = false
save_all_plots!(p::PipelinePlan, s::Symbol) = map(k -> p.save_plot[k]=true, collect(keys(p.save_plot)) )
save_no_plots!(p::PipelinePlan, s::Symbol) = map(k -> p.save_plot[k]=false, collect(keys(p.save_plot)) )

save_data(p::PipelinePlan, s::Symbol) = haskey(p.save_data,s) && p.save_data[s]
save_data!(p::PipelinePlan, s::Symbol) = p.save_data[s] = true
dont_save_data!(p::PipelinePlan, s::Symbol) = p.save_data[s] = false
save_all_data!(p::PipelinePlan, s::Symbol) = map(k -> p.save_data[k]=true, collect(keys(p.save_data)) )
save_no_data!(p::PipelinePlan, s::Symbol) = map(k -> p.save_data[k]=false, collect(keys(p.save_data)) )

need_to(p::PipelinePlan, s::Symbol) = haskey(p.need_to,s) && p.need_to[s]
function need_to!(p::PipelinePlan, s::Symbol)
    p.need_to[s] = true
    invalidate!(p,s)
end
dont_need_to!(p::PipelinePlan, s::Symbol) = p.need_to[s] = false
reset_all_needs!(p::PipelinePlan) = map(k -> p.need_to[k]=true, collect(keys(p.need_to)) )
reset_no_needs!(p::PipelinePlan) = map(k -> p.need_to[k]=false, collect(keys(p.need_to)) )

has_cache(p::PipelinePlan, s::Symbol ) = haskey(p.cache,s)
read_cache(p::PipelinePlan, s::Symbol ) = p.cache[s]
set_cache!(p::PipelinePlan, s::Symbol, data ) = p.cache[s] = data
reset_all_cache!(p::PipelinePlan ) = p.cache = needs_to_do_everything

# TODO: Eventually make invalidate! more intelligent, so it automatically invalidates other steps/cached data that depend on this step.
invalidate!(p::PipelinePlan, s::Symbol ) = delete!(p.cache, s)

end # module Pipeline
