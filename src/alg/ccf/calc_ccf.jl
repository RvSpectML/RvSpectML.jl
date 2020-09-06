"""
    Code to compute CCFs
Author: Michael Palumbo
Created: December 2019
Contact: mlp95@psu.edu
Based on code by Alex Wise (aw@psu.edu)
Refactors and optimized by Eric Ford
"""



"""
    ccf_1D!(ccf_out, λs, fluxes, line_list; mask_shape, plan )

Compute the cross correlation function of a spectrum with a mask.
# Inputs:
- ccf_out: 1-d array of size length(calc_ccf_v_grid(plan)) to store output
- λs: 1-d array of wavelengths
- fluxes:  1-d array of fluxes
- line_linst:  Each mask entry consists of the left and right end of tophats and a weight.
# Optional Arguments:
- mask_shape: shape of mask to use (currently only works with TopHatCCFMask())
- plan:  parameters for computing ccf (BasicCCFPlan())
"""
function ccf_1D!(ccf_out::A1, λ::A2, flux::A3,
    #, line_list::ALL, #mask_shape1::A3
                #; mask_shape::ACMS = TopHatCCFMask(),
                plan::PlanT = BasicCCFPlan() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1},
                PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
                # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape, AP<:AbstractCCFPlan }
    # make sure wavelengths and spectrum are 1D arrays
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    v_grid = calc_ccf_v_grid(plan)
    #ccf_out = zeros(size(v_grid))
    @assert size(ccf_out,1) == length(v_grid)
    # TODO: Opt could move this assertion into the ccf_plan
    @assert all( map(i->λ_min(plan.mask_shape,plan.line_list.λ[i+1]) .> λ_max(plan.mask_shape,plan.line_list.λ[i]) , 1:(length(plan.line_list.λ)-1) ) )

    mask_projections = zeros(length(λ),1)
    mask_workspace = zeros(length(λ)+2)

    # loop through each velocity, project the mask and compute ccf at given v
    for i in 1:size(ccf_out,1)
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask_opt!(mask_projections, λ, plan, shift_factor=doppler_factor, workspace=mask_workspace) #line_list)

        # compute the ccf value at the current velocity shift
        ccf_out[i] = sum(flux .* mask_projections)
    end
    return ccf_out
end


"""
    ccf_1D_expr!(ccf_out, λs, fluxes, line_list; mask_shape, plan )

Compute the cross correlation function of a spectrum with a mask.
    Experimental version that should work with different mask shapes.
    Need to understand why difference before merging this in.
# Inputs:
- ccf_out: 1-d array of size length(calc_ccf_v_grid(plan)) to store output
- λs: 1-d array of wavelengths
- fluxes:  1-d array of fluxes
- line_linst:  Each mask entry consists of the left and right end of tophats and a weight.
# Optional Arguments:
- mask_shape: shape of mask to use (TopHatCCFMask())
- plan:  parameters for computing ccf (BasicCCFPlan())
"""
function ccf_1D_expr!(ccf_out::A1, λ::A2, flux::A3,
    #, line_list::ALL, #mask_shape1::A3
                #; mask_shape::ACMS = TopHatCCFMask(),
                plan::PlanT = BasicCCFPlan() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1},
                PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
                # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape, AP<:AbstractCCFPlan }
    # make sure wavelengths and spectrum are 1D arrays
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    v_grid = calc_ccf_v_grid(plan)
    #ccf_out = zeros(size(v_grid))
    @assert size(ccf_out,1) == length(v_grid)
    # TODO: Opt could move this assertion into the ccf_plan
    @assert all( map(i->λ_min(plan.mask_shape,plan.line_list.λ[i+1]) .> λ_max(plan.mask_shape,plan.line_list.λ[i]) , 1:(length(plan.line_list.λ)-1) ) )

    mask_projections = zeros(length(λ),1)
    mask_workspace = zeros(length(λ)+2)

    # loop through each velocity, project the mask and compute ccf at given v
    for i in 1:size(ccf_out,1)
        # project shifted mask on wavelength domain
        doppler_factor = calc_doppler_factor(v_grid[i])
        project_mask_expr!(mask_projections, λ, plan, shift_factor=doppler_factor, workspace=mask_workspace) #line_list)

        # compute the ccf value at the current velocity shift
        ccf_out[i] = sum(flux .* mask_projections)
    end
    return ccf_out
end


"""
    ccf_1D!(ccf_out, λs, fluxes, line_list; mask_shape, plan )

Compute the cross correlation function of a spectrum with a mask.
# Inputs:
- λs: 1-d array of wavelengths
- fluxes:  1-d array of fluxes
- line_linst:  Each mask entry consists of the left and right end of tophats and a weight.
# Optional Arguments:
- mask_shape: shape of mask to use (TopHatCCFMask())
- plan:  parameters for computing ccf (BasicCCFPlan())
# Returns:
- 1-d array of size length(calc_ccf_v_grid(plan))
"""
function ccf_1D(λ::A2, flux::A3, #line_list::ALL, #mask_shape1::A3
                plan::PlanT = BasicCCFPlan() ) where {
                #; mask_shape::ACMS = TopHatCCFMask(), plan::PlanT = BasicCCFPlan() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1},
                #ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape, AP<:AbstractCCFPlan }
                PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    # TODO: Opt could move this assertion into the ccf_plan
    @assert all( map(i->λ_min(plan.mask_shape,plan.line_list.λ[i+1]) .> λ_max(plan.mask_shape,plan.line_list.λ[i]) , 1:(length(plan.line_list.λ)-1) ) )

    v_grid = calc_ccf_v_grid(plan)
    ccf_out = zeros(size(v_grid))
    ccf_1D!(ccf_out, λ, flux, plan )
    return ccf_out
end

function ccf_1D_expr(λ::A2, flux::A3, #line_list::ALL, #mask_shape1::A3
                plan::PlanT = BasicCCFPlan() ) where {
                #; mask_shape::ACMS = TopHatCCFMask(), plan::PlanT = BasicCCFPlan() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,1}, T3<:Real, A3<:AbstractArray{T3,1},
                #ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape, AP<:AbstractCCFPlan }
                PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
    @assert ndims(λ) == 1
    @assert ndims(flux) == 1
    @assert length(λ) == length(flux)
    @assert all( map(i->λ_min(plan.mask_shape,plan.line_list.λ[i+1]) .> λ_max(plan.mask_shape,plan.line_list.λ[i]) , 1:(length(plan.line_list.λ)-1) ) )

    v_grid = calc_ccf_v_grid(plan)
    ccf_out = zeros(size(v_grid))
    ccf_1D_expr!(ccf_out, λ, flux, plan )
    return ccf_out
end


"""
    project_mask!( output, λs, line_list; shift_factor, mask_shape )

Compute the projection of the mask onto the 1D array of wavelengths (λs) at a given shift factor (default: 1).
The mask is computed from a line_list and mask_shape (default: tophat).

TODO: Implement/test other mask_shapes.
"""
function project_mask_opt!(projection::A2, λs::A1, plan::PlanT ; shift_factor::Real=one(T1),
                            workspace::A3 = Array{eltype(λ),1}(undef,length(λ)+2) ) where {
    #     }, line_list::ALL
    #        mask_shape::ACMS = TopHatCCFMask() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T1,1},
                PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }

    # find bin edges
    #λsle = find_bin_edges_opt(λs)  # Read as λ's left edges
    @assert length(workspace) == length(λs)+2
    find_bin_edges_opt!(workspace,λs)  # Read as λ's left edges
    λsle = workspace
    # TODO: OPT: Can get remove this allocation and compute these as needed, at least for tophat mask.  But might be useful to keep for more non-const masks.

    # allocate memory for mask projection
    nm = length(plan.line_list) # size(mask, 1)
    nx = size(λsle, 1)-2
    ny = size(λsle, 2)
    @assert ny == 1                   # TODO:  Ask Alex what purpose of this is.  Maybe running CCF with multiple masks from the same lines at once?
    @assert size(projection,1) == nx
    @assert size(projection,2) == ny
    #projection = zeros(nx, ny)
    projection .= zero(eltype(projection))  # Important to zero out projection before accumulating there!

    # set indexing variables
    p = 1
    λsle_cur = λsle[p]   # Current pixel's left edge
    λsre_cur = λsle[p+1] # Current pixel's right edge
    m = 1
    on_mask = false
    #mask_lo = line_list.λ_lo[m] * shift_factor   # current line's lower limit
    #mask_hi = line_list.λ_hi[m] * shift_factor   # current line's upper limit
    mask_lo = λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor   # current line's lower limit
    mask_hi = λ_max(plan.mask_shape,plan.line_list.λ[m]) * shift_factor   # current line's upper limit

    mask_weight = plan.line_list.weight[m]            # current liine's weight
    # loop through lines in mask, weight mask by amount in each wavelength bin
    while m <= nm
        if !on_mask
            if λsre_cur > mask_lo
                if λsre_cur > mask_hi   # Pixel overshot this mask entry, so weight based on mask filling only a portion of the pixel,
                    projection[p] += mask_weight * (mask_hi - mask_lo) / (λsre_cur - λsle_cur)
                    m += 1
                    if m<=length(plan.line_list)      # Move to next line
                        mask_lo = λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                        mask_hi = λ_max(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                        mask_weight = plan.line_list.weight[m]
                    else         # We're out of lines, so exit
                        break
                    end
                else                       # Right edge of current pixel has entered the mask for this line, but left edge hasn't
                    projection[p] += mask_weight * (λsre_cur - mask_lo) / (λsre_cur - λsle_cur)
                    on_mask = true         # Indicate next pixel will hav something to contribute based on the current line
                    p+=1                   # Move to next pixel
                    λsle_cur = λsle[p]
                    λsre_cur = λsle[p+1]
                end
            else   # ALl of pixel is still to left of beginning of mask for this line.
                p+=1
                λsle_cur = λsle[p]
                λsre_cur = λsle[p+1]
            end
        else
            if λsre_cur > mask_hi               # Right edge of this pixel moved past the edge of the current line's mask
                projection[p] += mask_weight * (mask_hi - λsle_cur) / (λsre_cur - λsle_cur)
                on_mask = false                 # Indicate that we're done with this line
                m += 1                          # Move to next line
                if m<=length(plan.line_list)
                    mask_lo = λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                    mask_hi = λ_max(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                    mask_weight = plan.line_list.weight[m]
                else                            # We're done with all lines, can return early
                    break
                end
            else                                # Mask window is entirely within this pixel
                projection[p] += mask_weight    # Add the full mask weight
                p += 1                          # move to next pixel
                λsle_cur = λsle[p]
                λsre_cur = λsle[p+1]
            end
        end
        if p>length(projection) break end       # We're done with all pixels in this chunk.
    end
    return projection
end


""" project_mask_exper
TODO: Experimental code generalizing mask shape.
"""
function project_mask_expr!(projection::A2, λs::A1, plan::PlanT ; shift_factor::Real=one(T1),
                            workspace::A3 = Array{eltype(λ),1}(undef,length(λ)+2) ) where {
    #     }, line_list::ALL
    #        mask_shape::ACMS = TopHatCCFMask() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T1,1},
                PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }

    # find bin edges
    #λsle = find_bin_edges_opt(λs)  # Read as λ's left edges
    @assert length(workspace) == length(λs)+2
    find_bin_edges_opt!(workspace,λs)  # Read as λ's left edges
    λsle = workspace
    # TODO: OPT: Can get remove this allocation and compute these as needed, at least for tophat mask.  But might be useful to keep for more non-const masks.

    # allocate memory for mask projection
    nm = length(plan.line_list) # size(mask, 1)
    nx = size(λsle, 1)-2
    ny = size(λsle, 2)
    @assert ny == 1                   # TODO:  Ask Alex what purpose of this is.  Maybe running CCF with multiple masks from the same lines at once?
    @assert size(projection,1) == nx
    @assert size(projection,2) == ny
    #projection = zeros(nx, ny)
    projection .= zero(eltype(projection))  # Important to zero out projection before accumulating there!

    # set indexing variables
    p = 1
    λsle_cur = λsle[p]   # Current pixel's left edge
    λsre_cur = λsle[p+1] # Current pixel's right edge
    m = 1
    on_mask = false
    #mask_lo = line_list.λ_lo[m] * shift_factor   # current line's lower limit
    #mask_hi = line_list.λ_hi[m] * shift_factor   # current line's upper limit
    mask_mid = plan.line_list.λ[m] * shift_factor
    mask_lo = λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor   # current line's lower limit
    mask_hi = λ_max(plan.mask_shape,plan.line_list.λ[m]) * shift_factor   # current line's upper limit
    c_mps = RvSpectML.speed_of_light_mps
    mask_weight = plan.line_list.weight[m]
    mask_val_at_zero = plan.mask_shape(0.0)            # current liine's weight
    #mask_at_last_mask_lo = plan.mask_shape(mask_lo)
    # loop through lines in mask, weight mask by amount in each wavelength bin
    # assumes no overlap in mask entries
    while m <= nm
        if !on_mask
            if λsre_cur > mask_lo
                if λsre_cur > mask_hi   # Pixel overshot this mask entry, so weight based on mask filling only a portion of the pixel,
                    # TODO check this is where should evaluate plan.mask_shape().
                    integrand = 0.5*(0.5*(plan.mask_shape((mask_lo-mask_mid)/mask_mid*c_mps)+
                                          plan.mask_shape((mask_hi-mask_mid)/mask_mid*c_mps))+
                                          #plan.mask_shape((0.5*(mask_lo+mask_hi)-mask_mid)/mask_mid*c_mps))
                                          mask_val_at_zero )

                    projection[p] += (integrand * c_mps/mask_mid) * mask_weight * (mask_hi - mask_lo) # / (λsre_cur - λsle_cur)
                    m += 1
                    if m<=length(plan.line_list)      # Move to next line
                        mask_mid = plan.line_list.λ[m] * shift_factor
                        mask_lo = λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                        mask_hi = λ_max(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                        mask_weight = plan.line_list.weight[m]
                    else         # We're out of lines, so exit
                        break
                    end
                else                       # Right edge of current pixel has entered the mask for this line, but left edge hasn't
                     # TODO check this is where should evaluate plan.mask_shape().
                     # TODO warning, trapezoid rule not good for masks that go to zero at boundary.  Make plan specify how to integrate mask
                    integrand =  0.5*(0.5*(plan.mask_shape((mask_lo-mask_mid)/mask_mid*c_mps)+
                                           plan.mask_shape((λsre_cur-mask_mid)/mask_mid*c_mps))+
                                           plan.mask_shape((0.5*(mask_lo+λsre_cur)-mask_mid)/mask_mid*c_mps))
#                                           plan.mask_shape((sqrt(mask_lo*λsre_cur)-mask_mid)/mask_mid*c_mps))
                    projection[p] +=   (integrand * c_mps/mask_mid) * mask_weight * (λsre_cur - mask_lo)  # / (λsre_cur - λsle_cur)
                    on_mask = true         # Indicate next pixel will hav something to contribute based on the current line
                    p+=1                   # Move to next pixel
                    λsle_cur = λsle[p]
                    λsre_cur = λsle[p+1]
                end
            else   # ALl of pixel is still to left of beginning of mask for this line.
                p+=1
                λsle_cur = λsle[p]
                λsre_cur = λsle[p+1]
            end
        else
            if λsre_cur > mask_hi               # Right edge of this pixel moved past the edge of the current line's mask
                # TODO check this is where should evaluate plan.mask_shape().
                integrand =  0.5*(0.5*(plan.mask_shape((mask_hi-mask_mid)/mask_mid*c_mps)+
                                       plan.mask_shape((λsle_cur-mask_mid)/mask_mid*c_mps))+
                                       plan.mask_shape((0.5*(mask_hi+λsle_cur)-mask_mid)/mask_mid*c_mps))
#                                       plan.mask_shape((sqrt(mask_hi*λsle_cur)-mask_mid)/mask_mid*c_mps))
                projection[p] += (integrand * c_mps/mask_mid) * mask_weight * (mask_hi - λsle_cur) #/ (λsre_cur - λsle_cur)
                on_mask = false                 # Indicate that we're done with this line
                m += 1                          # Move to next line
                if m<=length(plan.line_list)
                    mask_mid = plan.line_list.λ[m] * shift_factor
                    mask_lo = λ_min(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                    mask_hi = λ_max(plan.mask_shape,plan.line_list.λ[m]) * shift_factor
                    mask_weight = plan.line_list.weight[m]
                else                            # We're done with all lines, can return early
                    break
                end
            else                                # Mask window is entirely within this pixel
                # assumes plan.mask_shape is normalized to integrate to unity and flux is constant within pixel
                #projection[p] += mask_weight    # Add the full mask weight
                @assert false
                projection[p] += mask_weight    # Add the full mask weight
                p += 1                          # move to next pixel
                λsle_cur = λsle[p]
                λsre_cur = λsle[p+1]
            end
        end
        if p>length(projection) break end       # We're done with all pixels in this chunk.
    end
    return projection
end




"""
    calc_doppler_factor(vel)

Compute the longitudinal relativistic doppler factor given a velocity
in meters per second.
"""
function calc_doppler_factor(vel::Real)
    one(vel) + vel/c_ms
end

#=
 Eric replaced this version with above.  Should double check that won't cause anyeone suprises.
function calc_doppler_factor_old(vel::Real)
    num = one(vel) + vel/c_ms
    den = one(vel) - vel/c_ms
    return sqrt(num/den)
end
=#


"""
   find_bin_edges( pixel_centers )

Internal function used by project_mask!.
TODO: OPT may be able to eliminate need for this memory allocation.  Or at least preallocate it.
"""
function find_bin_edges_opt(fws::A) where { T<:Real, A<:AbstractArray{T,1} }
    fwse = Array{eltype(fws),1}(undef,length(fws)+2)
    fwse[2:end-2] .= (fws[2:end] .+ fws[1:end-1]) ./ 2.0
    fwse[1] = 2.0 * fws[1] - fwse[2]
    fwse[end-1] = 2.0 * fws[end] - fwse[end-2]
    fwse[end] = zero(eltype(A))
    return fwse
end

function find_bin_edges_opt!(output::A1, fws::A2) where { T<:Real, A1<:AbstractArray{T,1}, A2<:AbstractArray{T,1} }
    #fwse = Array{eltype(fws),1}(undef,length(fws)+2)
    output[2:end-2] .= (fws[2:end] .+ fws[1:end-1]) ./ 2.0
    output[1] = 2.0 * fws[1] - output[2]
    output[end-1] = 2.0 * fws[end] - output[end-2]
    output[end] = zero(eltype(output))
    return output
end

#=
# Less optimized version of code in case need to debug.  Probably doesn't work with plans yet.
function find_bin_edges_compare(fws::A) where { T<:Real, A<:AbstractArray{T,1} }
    fwse = (fws[2:end] + fws[1:end-1]) ./ 2.0
    le = 2.0 * fws[1] - fwse[1]
    re = 2.0 * fws[end] - fwse[end]
    fwsle = vcat(le, fwse)
    fwsre = vcat(fwse, re)
    return fwsle, fwsre
end

function project_mask_compare!(projection::A2, λs::A1, plan::PlanT = BasicCCFPlan() ) where {
    #mask_in::ALL; shift_factor::Real=one(T1),
            #mask_shape::ACMS = TopHatCCFMask() ) where {
                T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2},
                PlanT<:AbstractCCFPlan } # ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }
                #plan::PlanT = BasicCCFPlan() ) where {
                #ALL<:AbstractLineList, ACMS<:AbstractCCFMaskShape }

#function project_mask_compare(fws::A1, mask::A2; shift_factor::Real=one(T1)) where { T1<:Real, A1<:AbstractArray{T1,1}, T2<:Real, A2<:AbstractArray{T2,2} }
    # find bin edges
    fwsle, fwsre = find_bin_edges_compare(λs)

    # allocate memory for line_list projection
    nm = length(mask_in)
    nx = size(fwsle, 1)
    ny = size(fwsle, 2)
    projection = zeros(nx, ny)
    #flush(stdout);  println("size(projection) = ", size(projection), " size(fws) = ", size(fws), " size(mask) = ", size(mask), " shift_factr = ", shift_factor)

    # shift the mask
    mask_shifted = zeros(length(mask_in),3)
    mask_shifted[:,1] = plan.mask_shape.λ_lo  .* shift_factor  # view(mask,:,1) .* shift_factor
    mask_shifted[:,2] = plan.mask_shape.λ_hi  .* shift_factor  # view(mask,:,2) .* shift_factor
    mask_shifted[:,3] = plan.mask_shape.weight                 # view(mask,:,3)
    local mask = mask_shifted

    # set indexing variables
    p = 1
    m = 1
    on_mask = false

    # loop through lines in mask, weight mask by amount in each wavelength bin
    while m <= nm
        if !on_mask
            if fwsre[p] > mask[m,1]
                if fwsre[p] > mask[m,2]
                    projection[p] += mask[m,3] * (mask[m,2] - mask[m,1]) / (fwsre[p] - fwsle[p])
                    m += 1
                else
                    projection[p] += mask[m,3] * (fwsre[p] - mask[m,1]) / (fwsre[p] - fwsle[p])
                    on_mask = true
                    p+=1
                end
            else
                p+=1
            end
        else
            if fwsre[p] > mask[m,2]
                projection[p] += mask[m,3] * (mask[m,2] - fwsle[p]) / (fwsre[p] - fwsle[p])
                on_mask = false
                m += 1
            else
                projection[p] += mask[m,3]
                p += 1
            end
        end
        if p>length(projection) break end
    end
    return projection
end
=#
