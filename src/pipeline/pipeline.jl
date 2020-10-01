module Pipeline

using ..RvSpectML

using RvSpectMLBase

using CSV, DataFrames #, Query
using Dates # for timing

include("extract_orders.jl")
export extract_orders

using ..EchelleCCFs
using ..EchelleInstruments
include("prep_line_list.jl")
export prepare_line_list

CCFs = EchelleCCFs
include("ccf_total.jl")
export ccf_total

include("ccf_orders.jl")
export ccf_orders

using ..EchelleCCFs.RVFromCCF
using Statistics
include("calc_rvs_from_ccf_total.jl")
export calc_rvs_from_ccf_total

end # module Pipeline
