struct State
	membrane::Matrix{Int}
	cytosol::Vector{Int}
end

function State(M::Model, totmembrane::Vector, cytosol::Vector; rng)
    Nspecies, Nsites = nspecies(M), nsites(M)
    membrane = zeros(nv(M.g), Nspecies)
    for m in 1:Nspecies
        for _ in 1:totmembrane[m]
            membrane[rand(rng, 1:Nsites), m] += 1
        end
    end
    State(membrane, cytosol)
end