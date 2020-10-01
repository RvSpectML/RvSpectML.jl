"""
Original author: Joe Ninan
Converted to Julia and optimized by Christian Gilbertson & Eric Ford
Returns a cubit interpolator for windowed sinc Filter curve.
no_of_points: number of intepolation points to use in cubic inteprolator
"""
function create_filter_curve_orig(no_of_points::Int; filter_size::Int=23, kaiserB::Int=13)
	x = range(-Int(floor(filter_size / 2)), stop=Int(floor(filter_size / 2)), length=no_of_points)
	Sinc = sinc.(x)
	Window = kaiser(length(x), kaiserB / π)
	FilterResponse = Window .* Sinc
	x_step = (x[end-1] - x[2]) / (no_of_points - 1)
	# append 0 to both ends far at the next node for preventing cubic spline
	# from extrapolating spurious values
	x = range(x[1] - x_step, stop=x[end] + x_step, length=no_of_points+2)
	#return CubicSplineInterpolation(x, multiple_append!([0.], FilterResponse, [0.]); extrapolation_bc=Line())
	return CubicSplineInterpolation(x, vcat(0., FilterResponse, 0.); extrapolation_bc=Line())
end

function create_filter_curve(no_of_points::Int; filter_size::Int=23, kaiserB::Int=13)
	ret_t = Float64
	x = range(-floor(Int,filter_size // 2), stop=floor(Int,filter_size // 2), length=no_of_points)
	# append 0 to both ends far at the next node for preventing cubic spline
	# from extrapolating spurious values
	Sinc = Array{ret_t,1}(undef,no_of_points+2)
	Window = Array{ret_t,1}(undef,no_of_points+2)
	Sinc[1] = zero(ret_t)
	Sinc[2:end-1] .= sinc.(x)
	Sinc[end] = zero(ret_t)
	Window[1] = zero(ret_t)
	Window[2:end-1] .= kaiser(length(x), kaiserB / π)
	Window[end] = zero(ret_t)
	FilterResponse = Window .* Sinc
	x_step = (x[end-1] - x[2]) / (no_of_points - 1)
	x = range(x[1] - x_step, stop=x[end] + x_step, length=no_of_points+2)
	return CubicSplineInterpolation(x, FilterResponse; extrapolation_bc=Line())
end

function next_min_arg(minarg::Integer, minvalue::Real, len::Integer)
	local sign = minvalue < 0
	local next_min_arg = minarg +sign -~sign
	if next_min_arg < 1       next_min_arg = len-next_min_arg end
	if next_min_arg > len     next_min_arg = next_min_arg-len end
	return next_min_arg
end

"""
Original author: Joe Ninan
Converted to Julia and optimized by Christian Gilbertson
Further adapted/optimized by Eric Ford
Additional optimizations possible by preallocating arrasy for minargs, Nminargs, minvalues, FilterValues and OldYCoords.
"""
function spectra_interpolate(
	newX::V1,
	oldX::V2,
	oldY::V3;
	PeriodicBoundary::Bool=false,
	filter_size::Int=23,
	kaiserB::Int=13,
	Filter::V4 = create_filter_curve(filter_size*21; filter_size=filter_size, kaiserB=kaiserB),
	 assume_sorted::Bool=false, verbose::Bool = false) where {
	 	T1<:Real, V1<:AbstractVector{T1}, T2<:Real, V2<:AbstractVector{T2}, T3<:Real, V3<:AbstractVector{T3}, T4<:Real, V4<:AbstractVector{T4}  }
	if !assume_sorted
		@assert issorted(newX)
		@assert issorted(oldX)
	end

	#Filter = create_filter_curve(filter_size*21; filter_size=filter_size, kaiserB=kaiserB)
	#pixarray = collect(-floor(Int,filter_size//2):floor(Int,filter_size//2))
	pixarray = -floor(Int,filter_size//2):floor(Int,filter_size//2)

	# oXsize = length(oldX)
    # First generate a 2D array of difference in pixel values
    # OldXminusNewX = np.array(oldX)[:,np.newaxis] - np.array(newX)
    # Find the minimum position to find nearest pixel for each each newX
	minargs = RvSpectMLBase.searchsortednearest(oldX, newX)
    # Pickout the those minumum values from 2D array
	minvalues = oldX[minargs] - newX
	Nminargs = next_min_arg.(minargs,minvalues,length(oldX))
	# EBF: used above to replace next few due to 1-based arrays, simplification and reduced temporaries
	#=
	sign = minvalues .< 0  # True means new X is infront of nearest old X
	# coordinate of the next adjacent bracketing point
	Nminargs = minargs .+sign -.~sign
	#Nminargs = Nminargs .% oXsize  # Periodic boundary  # WARN: Not right for 0-based arrays
	# In terms of pixel coordinates the shift values will be
	=#
	# EBF: consolidated to help reduce temporaries.  I'm not sure if if mattered.
	#=
	shiftvalues = minvalues./abs.(oldX[minargs]-oldX[Nminargs])
	# Coordinates to calculate the Filter values
	FilterCoords = shiftvalues .+ pixarray'
	FilterValues = Filter.(FilterCoords)
	=#
	FilterValues = Filter.(minvalues./abs.(oldX[minargs]-oldX[Nminargs]) .+ pixarray')
	# Coordinates to pick the values to be multiplied with Filter and summed
	OldYCoords = minargs .+ pixarray'
	if PeriodicBoundary
		# OldYCoords = OldYCoords .% oXsize  # Periodic boundary  # EBF: changed for 1-based arrays, TODO: WARN: need to test line below
		for i in 1:length(OldYCoords)
			if OldYCoords[i] < 1               		 OldYCoords[i] = length(oldX)-OldYCoords[i]   end
			if OldYCoords[i] > length(OldYCoords)+1  OldYCoords[i] = OldYCoords[i]-length(oldX)   end
		end

	else  # Extrapolate the last value till end..
		OldYCoords[OldYCoords .>= length(oldX) + 1] .= length(oldX)
		OldYCoords[OldYCoords .< 1] .= 1
	end
	# EBF consolidated to reduce temporaries.
	#OldYSlices = oldY[OldYCoords]
	#return vec(sum(OldYSlices .* FilterValues, dims=2))
	return vec(sum(oldY[OldYCoords] .* FilterValues, dims=2))
end
