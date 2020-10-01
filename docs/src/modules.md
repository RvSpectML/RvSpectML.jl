```@meta
CurrentModule = RvSpectML
```
# RvSpectML Modules

```@contents
Pages = ["modules.md"]
Depth = 3
```
## RV-Related Algorithms
```@autodocs
Modules = [ RvSpectML.DCPCA, RvSpectML.LineFinder ] #, RvSpectML.PPCA ]
Order = [:module]
```

## Interpolation Algorithms
```@autodocs
Modules = [RvSpectML.LinearInterpolation, RvSpectML.SincInterpolation, RvSpectML.TemporalGPInterpolation ]  # RvSpectML.GPInterpolation,
Order = [:module]
```

## Other Modules
```@autodocs
Modules = [Pipeline  ]
Order = [:module]
```
