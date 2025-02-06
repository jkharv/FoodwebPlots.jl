@recipe(FoodwebPlot, foodweb) do scene

    Attributes(
        node_weights = 1.0, # Scalar or Dict(String => Float64)
        edge_weights = 1.0, # Scalar or Dict(Tuple(2, String) => Float64)
        trophic_levels = false,
        draw_loops = true
    )
end

function Makie.plot!(fp::FoodwebPlot)

    s = richness(fp.foodweb[]) 
    g = SimpleDiGraph(s)

    edge_types = Union{
        eltype(SpeciesInteractionNetworks.interactions(fp.foodweb[])),
        eltype(edges(g))
    }

    conversion_dict = Dict{edge_types, edge_types}()

    for i in SpeciesInteractionNetworks.interactions(fp.foodweb[])

        if isloop(i) & !fp.draw_loops[]

            continue
        else
            
            src_index = findfirst(x -> x == object(i), species(fp.foodweb[]))
            dst_index = findfirst(x -> x == subject(i), species(fp.foodweb[]))
       
            e = Edge(src_index, dst_index)  

            conversion_dict[e] = i
            conversion_dict[i] = e

            add_edge!(g, e)
        end
    end

    initialpos, pin = nothing, nothing
    if fp.trophic_levels[]

        initialpos, pin = trophic_pin(fp.foodweb)
    else

        pin = [false for i in 1:s]
        initialpos = [Point2(rand(2)) for i in 1:s]
    end

    edge_weight_vec = nothing
    if !(typeof(fp.edge_weights[]) <: Number)

        edge_weight_vec = process_edge_weights(fp.edge_weights[], g, conversion_dict)
    end

    node_weight_vec = nothing
    if !(typeof(fp.node_weights[]) <: Number)

        node_weight_vec = process_node_weights(fp.node_weights[], g, fp.foodweb)
    end

    layout = Stress(;
        pin, 
        initialpos
    )

    graphplot!(fp, g; 
        layout = layout,
        edge_width = edge_weight_vec,
        arrow_size = 3*edge_weight_vec,
        node_size = node_weight_vec
    )

end

function process_edge_weights(edge_weights, g, conversion_dict)

    vec = Vector{Float32}()

    for e in edges(g)

        val = edge_weights[conversion_dict[e]]

        append!(vec, val)
    end

    return rescale(vec, 1, 8.0)
end

function process_node_weights(node_weights, g, foodweb)

    vec = Vector{Float32}()
    spp = species(foodweb[])

    for n in vertices(g)
              
        sp = spp[n]
       
        push!(vec, node_weights[sp])
    end

    return rescale(vec, 10.0, 30.0)
end

function rescale(vec, a, b)

    new_vec = log.(abs.(deepcopy(vec)))

    min = minimum(new_vec)
    max = maximum(new_vec)

    for i in eachindex(new_vec)

        new_vec[i] = a + (new_vec[i] - min)*(b - a)/(max - min)
    end

    return new_vec
end

function trophic_pin(foodweb)

    s = richness(foodweb[])

    yinit = distancetobase.(Ref(foodweb[]), species(foodweb[]), Ref(mean))
    xinit = [i * 1/s for i in 1:s]

    initialpos = Dict(collect(1:s) .=> Point2.(xinit, yinit))
    pin = Dict([i => (false, true) for i in 1:s])

    return initialpos, pin
end
