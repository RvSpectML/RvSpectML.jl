abstract type AbstractInstrument end
abstract type AbstractInstrument1D <: AbstractInstrument end
abstract type AbstractInstrument2D <: AbstractInstrument end

# TODO: Decide if these are really a good idea or if we shouldn't allow them.
""" Trait for genric 1D Extracted spectra """
struct Generic1D <: AbstractInstrument1D end

""" Trait for generic 2D Extracted spectra  """
struct Generic2D <: AbstractInstrument2D end
