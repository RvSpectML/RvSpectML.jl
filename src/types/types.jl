include("instruments.jl")
export AbstractInstrument, AbstractInstrument1D, AbstractInstrument2D

include("spectra.jl")
export AbstractSpectra, AbstractSpectra1D, AbstractSpectra2D
export Spectra1DBasic, Spectra2DBasic
export apply_doppler_boost!

include("chunks.jl")
export ChuckOfSpectrum, ChunkList, ChunkListTimeseries
export length, num_chunks
