```@meta
CurrentModule = RvSpectML
```
# RvSpectML
# Internals Types & Functions

```@contents
```

## General purpose
```@autodocs
Modules = [RvSpectML ]
Public = false
Order = [:type, :function ]
```

## Radial Velocity Related
```@autodocs
Modules = [ RvSpectML.CCF, RvSpectML.RVFromCCF, RvSpectML.DCPCA ] #, RvSpectML.PPCA ]
Public = false
Order = [:type, :function ]
```

## Interpolation
```@autodocs
Modules = [RvSpectML.LinearInterpolation, RvSpectML.SincInterpolation, RvSpectML.TemporalGPInterpolation ]  # RvSpectML.GPInterpolation,
Public = false
Order = [:type, :function]
```

## Instrument specific
```@autodocs
Modules = [RvSpectML.EXPRES, RvSpectML.HARPSN, RvSpectML.NEID, RvSpectML.TheoreticalInstrument  ]
Public = false
Order = [:type, :function]
```

