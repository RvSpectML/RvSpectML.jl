```@meta
CurrentModule = RvSpectML
```
# RvSpectML

## Getting Started

- [Install Julia 1.5](https://julialang.org/downloads/).  On Penn State's ICS-ACI, it is avaliable at  `/gpfs/group/ebf11/default/julia/bin/julia`.
- Install the [RvSpectML package](https://github.com/eford/RvSpectML.jl) and it's dependencies.  From julia
```julia
import Pkg
Pkg.add("https://github.com/eford/RvSpectML.jl")
Pkg.instantiate()
```
- Create a file `examples/data_paths.jl` specifying what directories on your system contain the relevant input data files.  For the two example scripts, you'd set  `solar_data_path` and/or `expres_data_path` like.
```
solar_data_path = "/gpfs/group/ebf11/default/ebf11/neid_solar/data"
expres_data_path = "/gpfs/group/ebf11/default/ebf11/expres/inputs/"
```
which are the paths to the required files for the examples ICS-ACI.

- Start julia in the RvSpectML directory and activate the associated Project.
```
> julia --project=.
```
- Run an example script or two.  E.g.,
```julia
include("examples/neid_pipeline_1.jl")
include("examples/expres_pipeline_1.jl")
```
- Tinker with some of the parameters in `examples/param.jl` or the example scripts.
- Let us know as you encounter any issues.
