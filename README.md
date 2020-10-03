# RvSpectML 
[![GitHub tag](https://img.shields.io/github/tag/RvSpectML/RvSpectML.jl.svg)](https://GitHub.com/RvSpectML/RvSpectML.jl/tags/)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://eford.github.io/RvSpectML.doc/) [![Build Status](https://github.com/RvSpectML/RvSpectML.jl/workflows/CI/badge.svg)](https://github.com/RvSpectML/RvSpectML.jl/actions)
<!--- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://RvSpectML.github.io/RvSpectML.jl/stable)  --->
<!--- [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://RvSpectML.github.io/RvSpectML.jl/dev) --->  
[![Coverage](https://codecov.io/gh/RvSpectML/RvSpectML.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/RvSpectML/RvSpectML.jl)


[RvSpectML.jl](https://github.com/RvSpectML/RvSpectML.jl) is a package to facilitate the analysis of stellar spectroscopic times series.
The primary goal is to measure extremely precise radial velocities (EPRVs).  
To support that goal, it will also include tools to deal with intrinsic stellar variability and telluric variability.  RvSpectML works with several other related packages.  

## Scope
[RvSpectML.jl](https://github.com/RvSpectML/RvSpectML.jl) is currently able to:
- Call [EchelleInstruments.jl](https://github.com/RvSpectML/EchelleInstruments.jl) to:
  - create a manifest of files (as a DataFrame) to be ingested from a directory (custom filtering via Query.jl)
  - read datafiles from NEID and EXPRES into a common set of data structures,
  - perform basic pre-processing (filtering out some orders, pixels within an order, chunks of spectra with NaNs, normalize spectra,...)
  - read a line list or cross-correlation function (CCF) mask file based on ESPRESSO or VALD,
- Call [EchelleCCFs.jl](https://github.com/RvSpectML/EchelleCCFs.jl):
  - compute cross-correlation function (CCF) of spectra relative to multiple CCF mask shapes efficiently,
  - measure RVs based on either the CCF or a direct Taylor expansion of the flux,
- interpolate spectra to a new set of wavelengths using linear, sinc, or Gaussian process regression algorithms,
- combine many files into a template spectra, interpolating them to a common wavelength grid and applying Doppler shift by estimated RV,
- perform Doppler-constrained PCA analysis.
- Call [RvSpectMLPlots.jl](https://github.com/RvSpectML/RvSpectMLPlots.jl) to:
   - make some common plots

[RvSpectML.jl and/or its companion packages](https://github.com/RvSpectML/RvSpectML-Overview) will eventually include tools to:
- read datafiles from additional spectrographs into a common set of data structures,
- perform additional pre-processing steps as needed,
- measure RVs using additional methods,
- calculated additional stellar activity indicators, and
- predict contamination due to stellar variability.

## Contributing
For now, please start by contributing code that you anticipate is likely useful for collaborators or other researchers.
Please keep code where you are actively experimenting with new approaches in separate github repositories.  Once you have a basic working example of how to apply your methods, then please create an example demonstrating that.  
Once it is reasonably mature, then please contact Eric to discuss whether to merge your code into this repo, one of the other associated repo or to keep it as an example showing how to use your method in its separate repository.  

## Related Packages & Repos
- [RvSpectMLBase](https://github.com/RvSpectML/RvSpectMLBase.jl): Types, common small utilities.  Minimal dependancies.  
- [EchelleInstruments.jl](https://github.com/RvSpectML/EchelleInstruments.jl): Code specific to each instrument
- [EchelleCCFs.jl](https://github.com/RvSpectML/EchelleCCFs.jl):  Computes CCFs with an anlytic mask
- RVSpectML (this package) holds larger algorithms and code that interconnects the component packages.  (Any plotting should be outside of src and not in the Project.toml.)
- [RvSpectMLPlots.jl](https://github.com/RvSpectML/RvSpectMLPlots.jl):  Plotting functions/scripts/notebooks, so other packages don't get bogged down by Plots.jl
- [Scalpels.jl](https://github.com/RvSpectML/Scalpels.jl):  Provides Scalpels algorithm for analyzing an ensemble of CCFs and estimating RVs and contamination from stellar variability.  
- [GPLinearODEMaker](https://github.com/christiangil/GPLinearODEMaker.jl):  Implements a multi-variate GP time-series likelihood and optimization functions.
