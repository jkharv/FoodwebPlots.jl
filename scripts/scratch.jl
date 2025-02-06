using FoodwebPlots
using HigherOrderFoodwebs
using DifferentialEquations
using SpeciesInteractionNetworks
using Statistics
using ModelingToolkit
using DataFrames
using GLMakie

web = nichemodel(30, 0.15)
fwm = FoodwebModel(web)

# ---------------------------------------------------------------- #
# Calculate some traits for each species to parametrize the model. #
# ---------------------------------------------------------------- #

traits = DataFrame(:species => species(fwm))
traits.vulnerability = vulnerability.(Ref(fwm.hg), traits.species);
traits.generality = generality.(Ref(fwm.hg), traits.species);
traits.trophic_level = distancetobase.(Ref(fwm.hg), traits.species, mean);

# Body Mass
mass_ratio = 1000;
f_mass(tl) = mass_ratio^tl;
traits.mass = f_mass.(traits.trophic_level);

# ------------------------------- #
# Create species-level parameters #
# ------------------------------- #

traits.metabolic_rate    = add_param!.(Ref(fwm), :x, traits.species,
    0.314 * traits.mass.^(-0.25)
);
traits.growth_rate       = add_param!.(Ref(fwm), :r, traits.species, 1.0);
traits.carrying_capacity = add_param!.(Ref(fwm), :k, traits.species, 1.0);
traits.max_consumption   = add_param!.(Ref(fwm), :y, traits.species, 4.0);

# --------------------------------------------------------- #
# Subset the interactions for different parts of the model. #
# --------------------------------------------------------- #

growth = filter(isloop, HigherOrderFoodwebs.interactions(fwm));
producer_growth = filter(x -> isproducer(fwm, subject(x)), growth);
consumer_growth = filter(x -> isconsumer(fwm, subject(x)), growth);
trophic = filter(!isloop, HigherOrderFoodwebs.interactions(fwm));

# --------------------------------------------------------- #
#  Set up some eij parameters for each trophic interaction  #
# --------------------------------------------------------- #

assimilation_efficiencies = Dict{Interaction, Num}();

F(B_focal, B_resources, a_resources, b0) = 
    (B_focal * a_resources[1]) / (b0 + sum(a_resources .* B_resources));

for intx ∈ trophic

    sym = Symbol("$(subject(intx))_$(object(intx))")
   
    s = subject(intx)
    o = object(intx)

    e_value = 0.45

    p = add_param!(fwm, :e, sym, e_value)
    
    assimilation_efficiencies[intx] = p
end

# ----------------------------------------- #
# Set up the dynamical aspects of the model #
# ----------------------------------------- #

for i ∈ producer_growth

    sbj = subject(i)
    s = fwm.conversion_dict[sbj]

    r = traits[traits.species .== sbj, :growth_rate][1]
    k = traits[traits.species .== sbj, :carrying_capacity][1]

    fwm.dynamic_rules[i] = DynamicRule(
        s*r*(1 - s/k)
    )
end

for i ∈ consumer_growth

    sbj = subject(i)
    s = fwm.conversion_dict[sbj]

    x = traits[traits.species .== sbj, :metabolic_rate][1]

    fwm.dynamic_rules[i] = DynamicRule(
        -x * s
    )
end

for i ∈ trophic

    s = fwm.conversion_dict[subject(i)]
    o = fwm.conversion_dict[object(i)]
    m = [fwm.conversion_dict[x] for x in with_role(:AF_modifier, i)]
    r = [o, m...]

    x = traits[traits.species .== subject(i), :metabolic_rate][1]
    y = traits[traits.species .== subject(i), :max_consumption][1]
    
    e = assimilation_efficiencies[i]

    ω = [1/length(r) for i in 1:length(r)]
    
    fwd = x * y * s * F(o, r, ω, 0.5) * e
    bwd = -fwd / e 

    fwm.dynamic_rules[i] = DynamicRule(fwd, bwd)
end

fwm = assemble_foodweb(fwm, Rosenbrock23());
prob = ODEProblem(fwm);

# Set up the callbacks
et = ExtinctionThresholdCallback(fwm, 1e-20);

# Simulate
sol = solve(prob, Rosenbrock23();
    reltol = 1e-3, 
    abstol = 1e-3,
    callback = et, 
    force_dtmin = true,
    maxiters = 1e6,
    tspan = 5000
);

foodwebtimeseries(fwm, sol)

net, fluxes, biomasses  = pairwise_network_sample(fwm, sol, tspan = 0:2000)
foodwebplot(net; 
    draw_loops = false, 
    trophic_levels = true,
    edge_weights = fluxes,
    node_weights = biomasses
)


function pairwise_network_sample(fwm::FoodwebModel, sol::ODESolution; 
    tspan = first(sol.t):last(sol.t))

    edge_weights  = Dict{AnnotatedHyperedge{Symbol}, Float64}()
    node_weights  = Dict{Symbol, Float64}()
    intxs = Vector{AnnotatedHyperedge{Symbol}}()
    spp = Set{Symbol}() # Set enforces uniqueness

    for intx in HigherOrderFoodwebs.interactions(fwm)
    
        s = subject(intx)
        o = object(intx)
        eq = fwm.dynamic_rules[intx].forwards_function

        # Average interaction strength over the sample period.
        edge_weight = mean(sol(tspan, idxs = eq))

        if abs(edge_weight) > 0.0
            
            push!(intxs, AnnotatedHyperedge([s,o], [:subject, :object]))
            edge_weights[intxs[end]] = edge_weight
        end
    end

    for sp in species(fwm)

        biomass = mean(sol(tspan, idxs = sp))

        if biomass > 0.0

            push!(spp, sp)
            node_weights[sp] = biomass
        end
    end

    return (
        SpeciesInteractionNetwork((Unipartite ∘ collect)(spp), intxs),
        edge_weights,
        node_weights
    )
end

