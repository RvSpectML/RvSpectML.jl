include("spectra.jl")
export AbstractSpectra, AbstractSpectra1D, AbstractSpectra2D
export Spectra1DBasic, Spectra2DBasic

include("chunks.jl")
export ChunkOfSpectra, ChunkList, ChunkListTimeseries
export length, num_chunks
