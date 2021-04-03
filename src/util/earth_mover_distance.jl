function wasserstein_distance(u_values::AbstractVector{T1}, v_values::AbstractVector{T2}, p::Integer) where { T1<:Real, T2<:Real}
  # code adapted from SciPy.stats._cdf_distance
  u_sorter = issorted(u_values) ? (1:length(u_values)) : sortperm(u_values)
  v_sorter = issorted(v_values) ? (1:length(v_values)) : sortperm(v_values)
  all_values = sort(vcat(u_values,v_values),alg=MergeSort)
  deltas = all_values[2:end] .- all_values[1:end-1]
  u_cdf = map(x->(searchsortedlast(u_values[u_sorter],x))/ length(u_values),all_values[1:end-1])
  v_cdf = map(x->(searchsortedlast(v_values[v_sorter],x))/ length(v_values),all_values[1:end-1])
  #println("deltas: ",deltas)
  #println("u_cdf:  ",u_cdf)
  #println("v_cdf:  ",u_cdf)
  #return nothing
  if p==1
     return sum(abs.(u_cdf.-v_cdf).*deltas)
  elseif p==2
      return sqrt(sum((u_cdf.-v_cdf).^2 .* deltas))
  else
      return pow(sum(abs.(u_cdf.-v_cdf).^p .* deltas),1/p)
  end
end

function wasserstein_distance(u_values::AbstractVector{T1}, v_values::AbstractVector{T2},
                              u_weights::AbstractVector{T3}, v_weights::AbstractVector{T4}, p::Integer
            ) where { T1<:Real, T2<:Real, T3<:Real, T4<:Real }
  # code adapted from SciPy.stats._cdf_distance
  @assert length(u_values) == length(u_weights)
  @assert length(v_values) == length(v_weights)
  u_sorter = issorted(u_values) ? (1:length(u_values)) : sortperm(u_values)
  v_sorter = issorted(v_values) ? (1:length(v_values)) : sortperm(v_values)
  all_values = sort(vcat(u_values,v_values),alg=MergeSort)
  deltas = all_values[2:end] .- all_values[1:end-1]
  u_cdf_indices = map(x->searchsortedlast(u_values[u_sorter],x),all_values[1:end-1]) .+1
  v_cdf_indices = map(x->searchsortedlast(v_values[v_sorter],x),all_values[1:end-1]) .+1
  u_sorted_weights = vcat(0.,cumsum(u_weights[u_sorter]))
  v_sorted_weights = vcat(0.,cumsum(v_weights[v_sorter]))
  u_cdf_max = u_sorted_weights[end]
  v_cdf_max = v_sorted_weights[end]
  u_cdf = u_sorted_weights[u_cdf_indices] / u_cdf_max
  v_cdf = v_sorted_weights[v_cdf_indices] / v_cdf_max
  if p==1
     return sum(abs.(u_cdf.-v_cdf).*deltas)
    elseif p==2
      return sqrt(sum((u_cdf.-v_cdf).^2 .* deltas))
    else
      return pow(sum(abs.(u_cdf.-v_cdf).^p .* deltas),1/p)
    end
end

function wasserstein_distance_presorted(u_values::AbstractVector{T1}, v_values::AbstractVector{T2},
                              u_weights::AbstractVector{T3}, v_weights::AbstractVector{T4}, p::Integer
            ) where { T1<:Real, T2<:Real, T3<:Real, T4<:Real }
  # code adapted from SciPy.stats._cdf_distance
  @assert length(u_values) == length(u_weights)
  @assert length(v_values) == length(v_weights)
  @assert issorted(u_values) && issorted(v_values)
  all_values = sort(vcat(u_values,v_values),alg=MergeSort)
  deltas = all_values[2:end] .- all_values[1:end-1]
  u_cdf_indices = map(x->searchsortedlast(u_values,x),all_values[1:end-1]) .+1
  v_cdf_indices = map(x->searchsortedlast(v_values,x),all_values[1:end-1]) .+1
  #return u_cdf_indices, v_cdf_indices, deltas, all_values
  u_sorted_weights = vcat(0.,cumsum(u_weights))
  v_sorted_weights = vcat(0.,cumsum(v_weights))
  u_cdf_max = u_sorted_weights[end]
  v_cdf_max = v_sorted_weights[end]
  u_cdf = u_sorted_weights[u_cdf_indices] / u_cdf_max
  v_cdf = v_sorted_weights[v_cdf_indices] / v_cdf_max
  if p==1
     return sum(abs.(u_cdf.-v_cdf).*deltas)
    elseif p==2
      return sqrt(sum((u_cdf.-v_cdf).^2 .* deltas))
    else
      return pow(sum(abs.(u_cdf.-v_cdf).^p .* deltas),1/p)
    end
end

function wasserstein_distance_presorted_common_x(u_values::AbstractVector{T1}, #v_values::AbstractVector{T2},
                              u_weights::AbstractVector{T3}, v_weights::AbstractVector{T4}, p::Integer
            ) where { T1<:Real, T2<:Real, T3<:Real, T4<:Real }
  # code adapted from SciPy.stats._cdf_distance
  @assert length(u_values) == length(u_weights)
  @assert length(u_values) == length(v_weights)
  @assert issorted(u_values) #&& issorted(v_values)
  all_values = zeros(T1,2*length(u_values))
  all_values[1:2:end-1] .= u_values
  all_values[2:2:end] .= u_values
  #deltas = all_values[2:end] .- all_values[1:end-1]
  deltas = zeros(T1,2*length(u_values)-1)
  deltas[1:2:end] .= zero(T1)
  deltas[2:2:end-1] .= u_values[2:end].-u_values[1:end-1]
  u_cdf_indices = map(x->searchsortedlast(u_values,x),all_values[1:end-1]) .+1
  #v_cdf_indices = map(x->searchsortedlast(v_values,x),all_values[1:end-1]) .+1
  u_sorted_weights = vcat(0.,cumsum(u_weights))
  v_sorted_weights = vcat(0.,cumsum(v_weights))
  u_cdf_max = u_sorted_weights[end]
  v_cdf_max = v_sorted_weights[end]
  u_cdf = u_sorted_weights[u_cdf_indices] / u_cdf_max
  v_cdf = v_sorted_weights[u_cdf_indices] / v_cdf_max
  if p==1
     return sum(abs.(u_cdf.-v_cdf).*deltas)
    elseif p==2
      return sqrt(sum((u_cdf.-v_cdf).^2 .* deltas))
    else
      return pow(sum(abs.(u_cdf.-v_cdf).^p .* deltas),1/p)
    end
end

#earth_mover_distance(u_values::AbstractVector{T1}, v_values::AbstractVector{T2}) where { T1<:Real, T2<:Real} = wasserstein_distance(u_values,v_values,1)
earth_mover_distance(u_values::AbstractVector{T1}, v_values::AbstractVector{T2}, u_weights::AbstractVector{T3}, v_weights::AbstractVector{T4}) where { T1<:Real, T2<:Real, T3<:Real, T4<:Real } = wasserstein_distance(u_values,v_values,u_weights,v_weights,1)
earth_mover_distance(u_values::AbstractVector{T1}, u_weights::AbstractVector{T3}, v_weights::AbstractVector{T4} ) where { T1<:Real, T3<:Real, T4<:Real } = wasserstein_distance_presorted_common_x(u_values,u_weights,v_weights,1)
