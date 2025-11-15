Base.@kwdef struct Model{G, S, T, L1, L2, L3, L4}
    species::S
    g::G
    rea::L1
    att::L2
    det::L3
    dif::L4
    rho_0::T
end

nspecies(M::Model) = length(M.species)

Base.length(M::Model) = nv(M.g)
