include("match.jl")

"""
Fallback MIP formulation for only 2-cycles and 2-3-cycles.
"""
function cycle_match(weights::Array{Float64,2};
                           graph::DiGraph=DiGraph(),
                           verbose::Int64 = 0,
                           max_cycle_length::Int64 = 2,
                           max_chain_length::Float64 = Inf64,
                           solver::String = "GLPK",
                           method::String = "fallback", #
                           ndds::Array{Int64, 1} = zeros(Int64, 0))

  if nv(graph) == 0
      graph = generate_graph(weights)
  end
  m, x, in = basic_MIP(weights, graph, ndds)
  n = nv(graph)
  if max_cycle_length == 2
      @constraint(m,
                  c5[i=1:n, j=intersect(outneighbors(graph, i), inneighbors(graph, i))],
                  x[Edge((i,j))] == x[Edge((j,i))])
      @constraint(m,
                  c6[i=1:n, j=setdiff(outneighbors(graph, i), inneighbors(graph, i))],
                  x[Edge((i,j))] == 0
                  )
  elseif max_cycle_length == 3
      @constraint(m,
                  c5[i=1:n,
                     j=outneighbors(graph, i),
                     k=intersect(outneighbors(graph, j), inneighbors(graph, i))],
                  x[Edge((k,i))] >= x[Edge((i,j))] + x[Edge((j,k))] - 1
                 )
     @constraint(m,
                 c6[i=1:n, j=outneighbors(graph, i), k=setdiff(outneighbors(graph, j), inneighbors(graph, i))],
                 Int(k == i) >= x[Edge((i,j))] + x[Edge((j,k))] - 1
                )
  else
      error("No fallback implemented for this setting")
  end
  status = optimize!(m)
  match_edges = JuMP.value.(x)
  match_vertices = JuMP.value.(in)
  value = JuMP.objective_value(m)

  return  (match_vertices, match_edges, value)
end
