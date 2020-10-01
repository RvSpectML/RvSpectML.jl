```@meta
CurrentModule = RvSpectML
```
# Types Exported by RvSpectML

```@contents
Pages = ["types.md"]
Depth = 3
```

## General purpose
```@autodocs
Modules = [RvSpectML ]
Private = false
Order = [:type]
```
## RV-Related Algorithms
```@autodocs
Modules = [RvSpectML.DCPCA, RvSpectML.LineFinder ] #, RvSpectML.PPCA ]
Private = false
Order = [:type]
```

## Interpolation Algorithms
```@autodocs
Modules = [RvSpectML.LinearInterpolation, RvSpectML.SincInterpolation, RvSpectML.TemporalGPInterpolation ]  # RvSpectML.GPInterpolation,
Private = false
Order = [:type]
```

## Other
```@autodocs
Modules = [Pipeline  ]
Private = false
Order = [:type]
```
