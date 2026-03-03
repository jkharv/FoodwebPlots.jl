using FoodwebPlots
using SpeciesInteractionNetworks
using AnnotatedHypergraphs
import WGLMakie

web = structuralmodel(NicheModel, 30, 0.15)

foodwebplot(web; draw_loops = false, 
    node_weights = Dict([sp => rand() for sp in species(web)]),
    edge_weights = Dict([(x[1], x[2]) => rand() for x in interactions(web)])
    )

interactions(web)


typeof(web)