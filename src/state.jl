struct State
	membrane::Matrix{Int}
	cytosol::Vector{Int}
end

function State(M::Model, totmembrane::Vector, cytosol::Vector; rng)
    Nspecies = nspecies(M)
    membrane = zeros(nv(M.g), Nspecies)
    for m in axes(membrane, 2)
        for _ in 1:totmembrane[m]
            membrane[rand(rng, axes(membrane,1)), m] += 1
        end
    end
    State(membrane, cytosol)
end