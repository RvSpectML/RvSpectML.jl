
# Declare functions that should be specialized for each instrument here, so they can be imported into their instrument's module.
function min_order end
function max_order end
function min_pixel_in_order end
function max_pixel_in_order end
function min_pixel end
function max_pixel end

function orders_to_use_default end
function min_col_default end
function max_col_default end

function metadata_symbols_default end
function metadata_strings_default end

# Trait-based functions that provide defaults (can be overwritten by instrument-specific versions)
orders_all(inst::AbstractInstrument2D) = min_order(inst):max_order(inst)
pixels_all(inst::AbstractInstrument2D) = min_pixels_in_order(inst):max_pixel_in_order(inst)
pixels_all(inst::AbstractInstrument1D) = min_pixel(inst):max_pixel(inst)
max_pixels_in_spectra(inst::AbstractInstrument1D) = length(pixels_all(inst))
max_pixels_in_spectra(inst::AbstractInstrument2D) = (max_order(inst)-min_order(inst)+1) * (max_pixel_in_order(inst)-min_pixel_in_order(inst)+1)
min_pixels_in_chunk(inst::AbstractInstrument1D) = 6
