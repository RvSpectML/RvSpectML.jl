```@meta
CurrentModule = RvSpectML
```
# RvSpectML Internals

As a heads up, these functions and types are more likely to change going forward than functions and types that are exported.  

```@contents
Pages = ["internals.md"]
Depth = 3
```
## Functions

### General purpose
```@autodocs
Modules = [RvSpectML ]
Public = false
Order = [ :function ]
```

### Radial Velocity Related
```@autodocs
Modules = [RvSpectML.DCPCA, RvSpectML.LineFinder ] #, RvSpectML.PPCA ]
Public = false
Order = [ :function ]
```

### Interpolation
```@autodocs
Modules = [RvSpectML.LinearInterpolation, RvSpectML.SincInterpolation, RvSpectML.TemporalGPInterpolation ]  # RvSpectML.GPInterpolation,
Public = false
Order = [ :function]
```

## Other
```@autodocs
Modules = [Pipeline  ]
Public = false
Order = [:function]
```

## Types

### General purpose
```@autodocs
Modules = [RvSpectML ]
Public = false
Order = [:type ]
```

### Radial Velocity Related
```@autodocs
Modules = [ RvSpectML.DCPCA, RvSpectML.LineFinder ] #, RvSpectML.PPCA ]
Public = false
Order = [:type ]
```

### Interpolation
```@autodocs
Modules = [RvSpectML.LinearInterpolation, RvSpectML.SincInterpolation, RvSpectML.TemporalGPInterpolation ]  # RvSpectML.GPInterpolation,
Public = false
Order = [:type ]
```

## Other
```@autodocs
Modules = [Pipeline  ]
Public = false
Order = [:type]
```
