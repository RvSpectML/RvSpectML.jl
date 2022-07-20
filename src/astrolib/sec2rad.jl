# This file is a part of AstroLib.jl. License is MIT "Expat".
# Copyright (C) 2016 Mosè Giordano.

"""
    sec2rad(sec) -> radians

### Purpose ###

Convert from seconds to radians.

### Argument ###

* `sec`: number of seconds.

### Output ###

The number of radians corresponding to `sec`.

### Example ###

### Notes ###

Use `rad2sec` to convert radians to seconds.
"""
sec2rad(sec::Real) = deg2rad(sec / 3600)
