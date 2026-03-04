using FoodwebPlots
using SpeciesInteractionNetworks
using AnnotatedHypergraphs
using Distributions
import WGLMakie

web = structuralmodel(NicheModel, 30, 0.15)

node_weights = Dict(species(web)[1:3] .=> rand(Uniform(1, 10), 3))
edge_weights = Dict(interactions(web)[1:50] .=> rand(LogNormal(0, 10), 50))

node_colours = Dict([sp => :black for sp in species(web)])

node_colours[:node_4] = :red

foodwebplot(web; draw_loops = false, 
    node_weights = node_weights,
    edge_weights = edge_weights,
    node_colours = node_colours,
)
    
rand(length(interactions(web)))
interactions(web)
typeof(web)