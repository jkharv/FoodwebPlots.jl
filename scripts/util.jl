function flux_t(hg::DynamicalHypergraph, sys::ODESystem, sol::ODESolution, t)

    cm = deepcopy(CommunityMatrix(hg))
    ints = EcologicalHypergraphs.interactions(hg)
   
    # Create a Param -> Val Dict
    pars = Dict()
    for p in sys.defaults

        if fix_other_peoples_dumb_ommissions(first(p), parameters(sys))

            pars[first(p)] = last(p)
        end
    end 

    # substitute in all the actual parameters values
    for i ∈ ints

        sbj = species(subject(i))[1]
        obj = species(object(i))[1] 

        cm[sbj, obj] = substitute(cm[sbj, obj], pars)
        cm[obj, sbj] = substitute(cm[obj, sbj], pars)
    end

    # Get all the state vars at t
    t_vals = Dict(collect(sys.states) .=> sol[sys.states, t])

    # substitute in all the actual parameters values
    for i ∈ ints

        sbj = species(subject(i))[1]
        obj = species(object(i))[1] 

        cm[sbj, obj] = substitute(cm[sbj, obj], t_vals)
        cm[obj, sbj] = substitute(cm[obj, sbj], t_vals)
    end

    ints = EcologicalHypergraphs.interactions(hg)
    flux_dict = Dict()
    for i ∈ ints

        sbj = species(subject(i))[1]
        obj = species(object(i))[1]

        if sbj != obj

            flux_dict[(sbj, obj)] = cm[sbj, obj]
        end
    end

    return flux_dict
end

function fix_other_peoples_dumb_ommissions(x, itr)

    for y in itr
    
        if isequal(x, y)
            return true
        end
    end
    return false
end
