Base.@kwdef struct Model{GT, S, T, L0, L1, L2, L3, L4}
    species::S
    G::GT
    rea::L0 = ()
    cat::L1 = ()
    att::L2 = ()
    det::L3 = ()
    dif::L4 = ()
    rho_0::T
end

macro species(ex...)
    return quote
        $(Expr(:tuple, (esc(x) for x in ex)...)) = Iterators.countfrom(1);
        $(Expr(:tuple, (esc(x) for x in ex)...))
    end
end

macro reaction(ex...)
    @assert length(ex) == 2
    constants = ex[1] isa Expr ? ex[1].args : ex[1]
    @assert length(ex[2].args) == 3
    ex1, ex2 = ex[2].args[2:3]
    substrates = ex1 isa Symbol ? ex1 : ex1.args[2:end] 
    products = ex2 isa Symbol ?  ex2 : ex2.args[2:end]
    if length(constants) == 1
        :($(Expr(:tuple, (esc(x) for x in substrates)...)), 
            $(Expr(:tuple, (esc(x) for x in products)...)),
            $(ex[1]))
    else # single molecule catalysis
        @assert length(products) == 2 && products[1] == substrates[1]
        :($((esc(x) for x in substrates)...), 
            $((esc(x) for x in products[2:end])...),
            $constants...)
    end

end

nspecies(M::Model) = length(M.species)

nsites(M::Model) = nv(M.G)
