struct Model{G}
    g::G
    posx::Vector{Float64}
    posy::Vector{Float64}
    dA::Float64
    dB::Float64
    dEA::Float64
    dEB::Float64
    kAc::Float64
    kBc::Float64
    kAa::Float64
    kAd::Float64
    kBa::Float64
    kBd::Float64
    KMM::Float64
    rho_0::Float64
end

Base.length(M::Model) = nv(M.g)

struct State{T}
	nA::Vector{T}
	nB::Vector{T}
	nEA::Vector{T}
	nEB::Vector{T}
	cytoEA::Base.RefValue{T}
	cytoEB::Base.RefValue{T}
end

Base.length(s::State) = length(s.nEA)


function plot(M::Model, s::State)
    L = floor(Int, sqrt(length(M)))
    color = RGB.(s.nA ./ 30, s.nB ./ 30, 0)
    x,y = M.posx./maximum(abs, M.posx), M.posy/maximum(abs, M.posy)
    compose(context(units=UnitBox(-1.2, -1.4, 2.4, 2.8)), ngon(x, y, fill(2/L/sqrt(3), length(M.posx)), fill(6, length(M.posx))), fill(color))
end
