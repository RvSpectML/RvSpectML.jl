""" Delegates loading of code specifying types essential to the package.  """

include("instruments.jl")
export AbstractInstrument, AbstractInstrument1D, AbstractInstrument2D
export Generic1D, Generic2D   # WARNING: Might remove these in future

include("spectra.jl")
export AbstractSpectra, AbstractSpectra1D, AbstractSpectra2D
export Spectra1DBasic, Spectra2DBasic
export AbstractSpectralTimeSeriesCommonWavelengths
export SpectralTimeSeriesCommonWavelengths
export make_vec_metadata_from_spectral_timeseries

include("chunks.jl")
export AbstractChuckOfSpectrum, ChuckOfSpectrum
export AbstractChunkList, ChunkList
export AbstractChunkListTimeseries, ChunkListTimeseries
export length, num_chunks
