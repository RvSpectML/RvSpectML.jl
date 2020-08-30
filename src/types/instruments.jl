"""
Code declaring abstract types for instruments' trait system.

Author: Eric Ford
Created: August 2020
"""


""" Abstract Base type for Instruments """
abstract type AbstractInstrument end

""" Abstract Base type for 1-D extracted spectra
    Should be specialized for specific instruments. """
abstract type AbstractInstrument1D <: AbstractInstrument end

""" Abstract Base type for 2-D extracted spectra.
    Should be specialized for each instrument """
abstract type AbstractInstrument2D <: AbstractInstrument end

# TODO: Decide if these are really a good idea or if we shouldn't allow them.
""" Trait for genric 1D Extracted spectra """
struct Generic1D <: AbstractInstrument1D end

""" Trait for generic 2D Extracted spectra  """
struct Generic2D <: AbstractInstrument2D end
