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
- Create a file `examples/data_paths.jl` specifying what directories on your system contain the relevant input data files.  For some of the first example scripts, you'd set  `expres_data_path` or `solar_data_path` or `ancilary_solar_data_path` like:
```
expres_data_path = "/gpfs/group/ebf11/default/ebf11/expres/inputs/"
solar_data_path = "/gpfs/group/ebf11/default/ebf11/neid_solar/data"
ancilary_solar_data_path = "/gpfs/group/ebf11/default/ebf11/neid_solar/data"
output_dir = joinpath(homedir(),"examples/output")
```
which are the paths to the required files for the examples ICS-ACI.  If you're saving outputs, then you'll likely want to set `output_dir`, too.

- Start julia in the RvSpectML directory and activate the associated Project.
```
> julia --project=.
```
- Run an example script or two.  E.g.,
```julia
include("examples/calc_rvs_ccf_std.jl")
```
- Tinker with some of the parameters in `examples/param.jl` or the example scripts.
- Let us know as you encounter any issues.
- If you intend to contribute to the RvSpectML package, then please fork the main repository, use Julia's package manager to add _your_ repo and set it into develop mode.  For a simpel bug fix, a simple pull request is probably ok.  For feature additions or non-trivial changes, please create a branch of your repo to use for the pull request.

## Other packages in the [RvSpectML ecosystem](https://github.com/RvSpectML):
- [RvSpectMLBase](https://rvspectml.github.io/RvSpectMLBase.jl/stable/)
- [RvSpectML](https://github.com/eford/RvSpectML.jl)
- [EchelleInstruments](https://rvspectml.github.io/EchelleInstruments.jl/stable/)
- [EchelleCCFs](https://rvspectml.github.io/EchelleCCFs.jl/stable)
- [Scalpels](https://rvspectml.github.io/Scalpels.jl/stable/)
- [RvSpectMLPlots](https://rvspectml.github.io/RvSpectMLPlots.jl/stable/)^[unreg]


 [^unreg]: This package is not yet registerd in Julia's general registry yet.
