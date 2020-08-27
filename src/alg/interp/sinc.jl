"""
   Sinc interpolation functionality
From
https://github.com/christiangil/GP-Research/blob/master/julia/src/interpolation_functions.jl
https://github.com/christiangil/GP-Research/blob/master/julia/src/general_functions.jl
"""

# interpolation_functions.jl

using DSP
using Interpolations

"""
Original author: Joe Ninan
Converted to Julia and optimized by Christian Gilbertson
Returns a cubit interpolator for windowed sinc Filter curve.
no_of_points: number of intepolation points to use in cubic inteprolator
"""
function create_filter_curve(no_of_points::Int; filter_size::Int=23, kaiserB::Int=13)
	x = linspace(-Int(floor(filter_size / 2)), Int(floor(filter_size / 2)), no_of_points)
	Sinc = sinc.(x)
	Window = kaiser(length(x), kaiserB / π)
	FilterResponse = Window .* Sinc
	x_step = (x[end] - x[1]) / (no_of_points - 1)
	# append 0 to both ends far at the next node for preventing cubic spline
	# from extrapolating spurious values
	x = linspace(x[1] - x_step, x[end] + x_step, no_of_points + 2)
	return CubicSplineInterpolation(x, multiple_append!([0.], FilterResponse, [0.]); extrapolation_bc=Line())
end


"""
Original author: Joe Ninan
Converted to Julia and optimized by Christian Gilbertson
"""
function spectra_interpolate(
	newX::Vector{T},
	oldX::Vector{T},
	oldY::Vector{T};
	PeriodicBoundary::Bool=false,
	filter_size::Int=23,
	kaiserB::Int=13) where T<:Real

	@assert issorted(newX)
	@assert issorted(oldX)

	Filter = create_filter_curve(filter_size*21; filter_size=filter_size, kaiserB=kaiserB)
	pixarray = collect(-Int(floor(filter_size/2)):Int(floor(filter_size/2)))

	oXsize = length(oldX)
    # First generate a 2D array of difference in pixel values
    # OldXminusNewX = np.array(oldX)[:,np.newaxis] - np.array(newX)
    # Find the minimum position to find nearest pixel for each each newX
	minargs = searchsortednearest(oldX, newX)
    # Pickout the those minumum values from 2D array
	minvalues = oldX[minargs] - newX
	sign = minvalues .< 0  # True means new X is infront of nearest old X
	# coordinate of the next adjacent bracketing point
	Nminargs = minargs .+sign -.~sign
	Nminargs = Nminargs .% oXsize  # Periodic boundary
	# In terms of pixel coordinates the shift values will be
	shiftvalues = minvalues./abs.(oldX[minargs]-oldX[Nminargs])
	# Coordinates to calculate the Filter values
	FilterCoords = shiftvalues .+ pixarray'
	FilterValues = Filter.(FilterCoords)
	# Coordinates to pick the values to be multiplied with Filter and summed
	OldYCoords = minargs .+ pixarray'
	if PeriodicBoundary
		OldYCoords = OldYCoords .% oXsize  # Periodic boundary
	else  # Extrapolate the last value till end..
		OldYCoords[OldYCoords .>= oXsize + 1] .= oXsize
		OldYCoords[OldYCoords .< 1] .= 1
	end
	OldYSlices = oldY[OldYCoords]
	return vec(sum(OldYSlices .* FilterValues, dims=2))
end
