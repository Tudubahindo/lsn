using TravelingSalesmanExact, SCIP
using CSV, DataFrames
using DelimitedFiles
map = DataFrame(CSV.File("map_square.dat"; delim='\t', header=false))
cities = Tuple.(eachrow(map))
set_default_optimizer!(SCIP.Optimizer)
n = 50
tour, cost = get_optimal_tour(cities; verbose = true)
plot_cities(cities[tour])
print(tour, '\n')
writedlm("output_bruteforce.dat", tour, '\t')
