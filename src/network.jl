@recipe(FoodwebPlot, foodweb) do scene

    Attributes(
        node_weights = 1.0, # Scalar or Dict(String => Float64)
        edge_weights = 1.0, # Scalar or Dict(Tuple(2, String) => Float64)
        draw_loops = true
    )
end

function Makie.plot!(fp::FoodwebPlot)

    web_observable = fp.foodweb
    S = richness(fp.foodweb[])

    g = SimpleDiGraph(S)

    for i in HigherOrderFoodwebs.interactions(fp.foodweb[])

        src_index = findfirst(x -> x == object(i), species(fp.foodweb[]))
        dst_index = findfirst(x -> x == subject(i), species(fp.foodweb[]))

        if isloop(i) & fp.draw_loops[]

            add_edge!(g, src_index, dst_index)
        else

            continue # Break this iteration and move onto the next inteaction.
        end

        add_edge!(g, src_index, dst_index)
    end

    # Y position is set to the trophic level.
    init_pos = trophicpositions(web_observable)
    # Pin nodes y position but allow them to move along x.
    pin = Dict(1:richness(web_observable[]) .=> Ref((false, true)));

    # Setup network layout alg with our new pinning requirements.
    foodweb_layout = Spring(;initialpos = init_pos, pin = pin, seed = 1);
   
    spp_biomass = node_weight_vector(web_observable[], fp.node_weights[])
    biomass_flux = edge_weight_vector(web_observable[], g, fp.edge_weights[])

    edge_colors = map(x -> x > 10E-5 ? "black" : "lightgrey", biomass_flux)
    node_colors = map(x -> x > 10E-5 ? "black" : "lightgrey", spp_biomass)

    edge_weights = map(x -> exp(x), biomass_flux)
    node_weights = map(x -> 15 + exp(x), spp_biomass)
    
    graphplot!(fp, g;

        layout = foodweb_layout, 
        node_size = node_weights,
        edge_width = edge_weights,
        arrow_size = 10*edge_weights,
        edge_color = edge_colors,
        # edge_attr = ,
        # node_attr = ,
        node_color = node_colors
    )

    return fp
end

""" 
Handles taking user input on edge weights and turning that into a Vector
defining the width of every edge.  
"""
function edge_weight_vector(fw, g, weight_attr)

    if weight_attr isa Dict

        spp = species(fw)
        weights = Vector{Float64}(undef, length(edges(g)))
    
        for (i, e) in enumerate(edges(g))
  
            sbj_sp = spp[src(e)]
            obj_sp = spp[dst(e)]

            weights[i] = weight_attr[(sbj_sp, obj_sp)].val
        end 
    
        return weights
    elseif weight_attr isa Number

        weights = [weight_attr for i ∈ 1:length(edges(g))]
        return weights
    else
        error("Improper format for edge weights")
    end

end

"""
Handles taking node weight input from the user and turning that into
a vector of sizes for each node.
"""
function node_weight_vector(fw, weight_attr)

    if weight_attr isa Dict

        spp = species(fw)
        weights = map(x -> weight_attr[x], spp) 
        return weights
    elseif weight_attr isa Number

        weights = [weight_attr for i ∈ 1:length(species(fw))]
        return weights
    else
        error("Improper format for node weights")
    end
end

function trophicpositions(fw)

    # Trophic level for all the species in the web.
    tls = distancetobase.(Ref(fw[]), species(fw[]), Ref(mean))

    # Set initial node y positions to their trophic level.
    spp_to_id = Dict(species(fw[]) .=> 1:richness(fw[]));
    init_pos = Dict{Int64, Point2f}();

    for spp ∈ species(fw[])

        init_pos[spp_to_id[spp]] = Point2f(1.0, tls[spp_to_id[spp]])
    end

    return init_pos
end