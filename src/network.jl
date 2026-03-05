@recipe(FoodwebPlot, foodweb) do scene

    Attributes(
        node_weights = 1.0, # Float64 or Dict{Edge, Float64} 
        edge_weights = 1.0, # Float64 or Dict{Edge, Float64} 
        edge_colours = :black,
        node_colours = :black,
        trophic_levels = true,
        draw_loops = false 
    )
end

function Makie.plot!(
    fp::FoodwebPlot{Tuple{<:SpeciesInteractionNetwork{Unipartite{Symbol}, Quantitative{Float64}}}})
  
    intxs = SpeciesInteractionNetworks.interactions(fp.foodweb[])   
    weights = Dict(intxs .=> [last(x) for x in intxs])

    println(weights)

    x = render(Binary, fp.foodweb[])

    println(typeof(x))

    foodwebplot(x,
        node_weights = fp.node_weights[], 
        edge_weights = weights,
        edge_colours = fp.edge_colours[],
        node_colours = fp.node_colours[],
        trophic_levels = fp.trophic_levels[],
        draw_loops = fp.draw_loops[]
    )    
end
 
function Makie.plot!(fp::FoodwebPlot{Tuple{<:SpeciesInteractionNetwork{Unipartite{Symbol}, Binary{Bool}}}})
    
    s = richness(fp.foodweb[])

    g, conversion_dict = convert_network_type(fp.foodweb[], fp.draw_loops[])

    edge_weights = process_edge_weights(fp, g, conversion_dict)
    node_weights = process_node_weights(fp, g)

    edge_colours = process_edge_colours(fp, g, conversion_dict)
    node_colours = process_node_colours(fp, g)

    # Whether to lock the y position to the trophic level of a species.
    initialpos, pin = nothing, nothing
    if fp.trophic_levels[]

        initialpos, pin = trophic_pin(fp.foodweb)
    else

        pin = [false for i in 1:s]
        initialpos = [Point2(rand(2)) for i in 1:s]
    end

    layout = Stress(;
        pin, 
        initialpos
    )

    graphplot!(fp, g; 
        layout = layout,
        edge_width = edge_weights,
        arrow_size = 5*edge_weights,
        node_size = node_weights,
        edge_color = edge_colours,
        node_color = node_colours
    )

end

function process_edge_colours(fp, g, conversion_dict)

    if fp.edge_colours[] isa Symbol  

        return [fp.edge_colours[] for e in edges(g)]
    end

    colours = Vector{Symbol}()

    for e in edges(g)

        if haskey(fp.edge_colours[], [conversion_dict[e]])

            push!(weights, fp.edge_colours[][conversion_dict[e]])
        else
            push!(weights, :black)
        end
    end

    return colours 
end

function process_edge_weights(fp, g, conversion_dict)

    if fp.edge_weights[] isa Number 

        return [fp.edge_weights[] for e in edges(g)]
    end

    weights = Vector{Float32}()
    mean_weight = (mean ∘ values)(fp.edge_weights[])

    for e in edges(g)

        if haskey(fp.edge_weights[], [conversion_dict[e]])

            push!(weights, fp.edge_weights[][conversion_dict[e]])
        else
            push!(weights, mean_weight)
        end
    end

    return rescale(weights, 0.2, 5.0)
end

function process_node_weights(fp, g)

    if fp.node_weights[] isa Number 

        return [fp.node_weights[] for v in vertices(g)]
    end

    vec = Vector{Float32}()
    mean_weight = (mean ∘ values)(fp.node_weights[])
    spp = species(fp.foodweb[])

    for n in vertices(g)
              
        sp = spp[n]
        
        if haskey(fp.node_weights[], sp) 

            push!(vec, fp.node_weights[][sp])
        else

            push!(vec, mean_weight)
        end
    end

    return rescale(vec, 5.0, 30.0)
end

function process_node_colours(fp, g)

    if fp.node_colours[] isa Symbol

        return [fp.node_colours[] for v in vertices(g)]
    end

    vec = Vector{Symbol}()
    spp = species(fp.foodweb[])

    for n in vertices(g)
              
        sp = spp[n]
        
        if haskey(fp.node_colours[], sp) 

            push!(vec, fp.node_colours[][sp])
        else

            push!(vec, :black)
        end
    end

    return vec
end

function convert_network_type(web::SpeciesInteractionNetwork, loops::Bool)

    s = richness(web)
    g = SimpleDiGraph(s)
  
    edge_types = Union{
        eltype(SpeciesInteractionNetworks.interactions(web)),
        eltype(edges(g))
    }

    conversion_dict = Dict{edge_types, edge_types}()

    for i in SpeciesInteractionNetworks.interactions(web)

        if (i[1] == i[2]) & !loops

            continue
        else

            src_index = findfirst(x -> x == i[2], species(web))
            dst_index = findfirst(x -> x == i[1], species(web))
       
            e = Edge(src_index, dst_index)  

            conversion_dict[e] = i
            conversion_dict[i] = e

            add_edge!(g, e)
        end
    end   

    return (g, conversion_dict)
end

function rescale(vec, a, b)

    new_vec = log.(abs.(copy(vec)))

    min = minimum(new_vec)
    max = maximum(new_vec)

    for i in eachindex(new_vec)

        new_vec[i] = a + (new_vec[i] - min)*(b - a)/(max - min)
    end

    return new_vec
end

function trophic_pin(foodweb)

    s = richness(foodweb[])

    yinit = distancetobase.(Ref(foodweb[]), species(foodweb[]), Ref(mean)) .+ 1.0
    xinit = [i * 1/s for i in 1:s]

    for (i, sp) in (enumerate ∘ species)(foodweb[])

        if generality(foodweb[], sp) == 0        

            yinit[i] = 1
        end
    end

    initialpos = Dict(collect(1:s) .=> Point2.(xinit, yinit))
    pin = Dict([i => (false, true) for i in 1:s])

    return initialpos, pin
end