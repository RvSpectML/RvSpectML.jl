# AstroLib.jl[https://github.com/JuliaAstro/AstroLib.jl] was breaking compatability with other packages
# So this is a hack to include the just needed functions and avoid extra dependancies

module AstroLib

using Dates
const JULIANCENTURY = 36_525
const ct2lst_c  = (280.46061837, 360.98564736629, 0.000387933, 38710000.0)
include("jdcnv.jl")
const J2000 = Int(jdcnv(2000, 01, 01, 12)) # 2451545

include("sec2rad.jl")
include("ct2lst.jl")
export ct2lst
include("hadec2altaz.jl")
export hadec2altaz
include("sunpos.jl")
export sunpos

end
