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
Modules = [RvSpectML.CCF, RvSpectML.RVFromCCF, RvSpectML.DCPCA, RvSpectML.Scalpels, RvSpectML.LineFinder ] #, RvSpectML.PPCA ]
Order = [:module]
```

## Interpolation Algorithms
```@autodocs
Modules = [RvSpectML.LinearInterpolation, RvSpectML.SincInterpolation, RvSpectML.TemporalGPInterpolation ]  # RvSpectML.GPInterpolation,
Order = [:module]
```

## Instrument-specific Modules
```@autodocs
Modules = [EXPRES, HARPSN, NEID, TheoreticalInstrument  ]
Order = [:module]
```

## Other Modules
```@autodocs
Modules = [Pipeline  ]
Order = [:module]
```
