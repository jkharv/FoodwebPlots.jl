@recipe(FoodwebTimeSeries, foodweb, solution) do scene

    Attributes()
end

function Makie.plot!(ts::FoodwebTimeSeries)

    for sp in species(ts.foodweb[])

        # TODO Make this find every positive population interval and draw a line
        # _segment_ for it.
        ext = findfirst(x -> x == 0.0, ts.solution[][sp])
        if isnothing(ext)

            ext = length(ts.solution[][sp])
        end

        lines!(ts, ts.solution[].t[1:ext], ts.solution[][sp][1:ext])
        
        # TODO Find all the upcrossings and downcrossings and use different
        # arrows to distinguish them.
        scatter!(ts, ts.solution[].t[ext], -1.0, 
            marker = '↘', 
            markersize = 20)

        # TODO Optional labelling of all the invasion and extinction events with
        # the species name

    end

    return ts
end

