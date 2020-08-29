# RvSpectML [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://eford.github.io/RvSpectML.jl/stable) [![Build Status](https://github.com/eford/RvSpectML.jl/workflows/CI/badge.svg)](https://github.com/eford/RvSpectML.jl/actions) 
<!--- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://eford.github.io/RvSpectML.jl/dev) --->  
<!--- [![Coverage](https://codecov.io/gh/eford/RvSpectML.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/eford/RvSpectML.jl) --->

[RvSpectML.jl](https://github.com/eford/RvSpectML.jl) is a package to facilitate the analysis of stellar spectroscopic times series.
The primary goal is to measure extremely precise radial velocities (EPRVs).  
To support that goal, it will also include tools to deal with intrinsic stellar variability and telluric variability.  

## Scope
[RvSpectML.jl](https://github.com/eford/RvSpectML.jl)
It will eventually include tools to:
- read datafiles from multiple spectrographs into a common format,
- perform basic pre-processing,
- measure RVs, and
- predict contamination due to stellar variability.

## Contributing
For now, please start by contributing code that you anticipate is likely useful for collaborators or other researchers.
Please  keep code where you are actively experimenting with new approaches in separate github repositories.  Once you have a basic working example of how to apply your methods, then please create an example demonstrating that.  
Once it is reasonably mature, then please contact Eric to discuss whether to merge your code into this repo or to leave it as an example showing how to use your method in its separate repository.  At some point we may break [RvSpectML.jl](https://github.com/eford/RvSpectML.jl) into multiple packages.
