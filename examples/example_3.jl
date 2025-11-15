using MembraneRD, ExponentialQueues
using Compose, MembraneRD.Colors
import MembraneRD: hexagon
using MembraneRD.Graphs

Base.@kwdef struct Species
    d::Float64 # diffusion of molecule in membrane
    kd::Float64 # detachment from membrane
end


Base.@kwdef struct ModelN{G, L1, L2}
    g::G
    posx::Vector{Float64}
    posy::Vector{Float64}
    species::Vector{Species}
    reactions::L1
    ka::L2
    rho_0::Float64
end

nspecies(M::ModelN) = length(M.species)

Base.length(M::ModelN) = nv(M.g)

struct StateN
	membrane::Matrix{Int}
	cytosol::Vector{Int}
end

function StateN(g, totmembrane::Vector, cytosol::Vector; rng)
    Nspecies = length(totmembrane)
    membrane = zeros(nv(g), Nspecies)
    for m in axes(membrane, 2)
        for _ in 1:totmembrane[m]
            membrane[rand(rng, axes(membrane,1)), m] += 1
        end
    end
    StateN(membrane, cytosol)
end



function run_RDN!(state::StateN, M::ModelN, T; 
        stats = (_, _)->nothing, 
        rng = Random.default_rng())

    Nspecies = nspecies(M)
    Qn = ((ExponentialQueue(length(M)) for _ in 1:Nspecies)...,)
    Qcat = [ExponentialQueue(length(M)) for _ in M.reactions]
    Qatt = [ExponentialQueue(length(M))*0.0 for _ in M.ka]
    Qdet = ((Qn[m]*M.species[m].kd for m in 1:Nspecies)...,)

    function update(i::Int)
        for ((e,s,_,_,km),q) in zip(M.reactions, Qcat)
            q[i] = state.membrane[i,e] * state.membrane[i,s] / (state.membrane[i,s] + km)
        end
        for ((m,m1,_),q) in zip(M.ka, Qatt)
            q.q[i] = state.membrane[i,m1]
            q.f[] = state.cytosol[m]
        end
        for m in 1:Nspecies
            Qn[m][i] = state.membrane[i,m]
        end
    end


    foreach(update, 1:length(M))

    #arrival is chosen uniformly between its neighbours
    rand_neighbor(i) = rand(rng, neighbors(M.g, i))

    Q = NestedQueue(
            ((:dif,m) => Qn[m] * M.species[m].d for m in 1:Nspecies)...,
            ((:att,m) => q*ka for ((m,m1,ka),q) in zip(M.ka,Qatt))...,
            ((:det,m) => Qdet[m] for m in 1:Nspecies)...,
            ((:cat,(s,p)) => kc*q for ((_,s,p,kc,_),q) in zip(M.reactions,Qcat))...,
        )

    #return Q

    println("starting simulation, $(length(Q)) events in the queue")

    t::Float64 = 0.0
    while !isempty(Q)
        ((ev, m::Union{Int,Tuple{Int,Int}}), i), dt = peek(Q; rng)
        t += dt
        t > T && break # reached end time for simulation
        stats(t, state)
        @inbounds if ev === :dif # diffusion
            j = rand_neighbor(i)
            state.membrane[i,m] -= 1
            state.membrane[j,m] += 1
            update(i)
            update(j)
        elseif ev === :cat # reaction
            s, p = m
            state.membrane[i,s] -= 1
            state.membrane[i,p] += 1
            update(i)
        elseif ev === :att #attachment to membrane
            state.cytosol[m] -= 1
            state.membrane[i,m] += 1
            update(i)
        elseif ev === :det #detachment from membrane
            state.membrane[i,m] -= 1
            state.cytosol[m] += 1
            update(i)
        end
    end
	stats(T, state)
end
	
using ProgressMeter, JLD, Random



function build_model_state(L; rng = Random.default_rng())
    g,posx,posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    unit_length = 1.0
    theta = 0.2126944621086619#Stot/V
    Ac = unit_length^2#area associated to each cell, set unit lengthscale
    Stot = N*Ac #area of each cell is assumed to be 1, setting unit lengthscale
    V = Stot/theta
    d_timescale = 0.01#this sets unit timescale
    dA,dB,dEA,dEB = d_timescale, d_timescale, d_timescale, d_timescale
    #rates in the theory (do not correspond excatly to those used in the simulations), some slight dimensional changes are needed
    kAa_th = 1.0
    kAd_th = 1.0*kAa_th
    kAc_th = 1.0
    kBa_th = 1.0
    kBd_th = 1.0*kBa_th
    kBc_th = 1.3*kAc_th
    KMMA_th = 1.0
    KMMB_th = 1.0
    #rates to implement
    kAc = kAc_th
    kBc = kBc_th
    kAa = kAa_th/V
    kAd = kAd_th
    kBa = kBa_th/V
    kBd = kBd_th
    KMMA = KMMA_th*Ac
    KMMB = KMMB_th*Ac

    A = Species(d=dA,kd=0.0)
    B = Species(d=dB,kd=0.0)
    EA = Species(d=dEA,kd=kAd)
    EB = Species(d=dEB,kd=kBd)


    reactions = ((3,2,1,kAc,KMMA),(4,1,2,kBc,KMMB))
    ka = ((3,1,kAa), (4,2,kBa))

    M = ModelN(; g, posx, posy, reactions, ka, species=[A,B,EA,EB], rho_0 = 0.0)

    totmol = N * 10
    totA, totB = floor(Int, totmol / 2), floor(Int, totmol / 2)
    totEA, totEB = floor(Int, 0.1*N), floor(Int, 0.1*N)
    memEA = floor(Int, totEA*(theta/(kAd_th/kAa_th))*(totA/Stot)/(1+((theta/(kAd_th/kAa_th))*(totA/Stot))))
    memEB = floor(Int, totEB*(theta/(kBd_th/kBa_th))*(totB/Stot)/(1+((theta/(kBd_th/kBa_th))*(totB/Stot))))
    cytoEA, cytoEB = totEA - memEA, totEB - memEB
    s = StateN(g, [totA, totB, memEA, memEB], [0.0, 0.0, cytoEA, cytoEB]; rng)
    M,s
end



function plot(M::ModelN, s::StateN)
    posx, posy = mm .* M.posx, mm .* M.posy
    (x0,x1),(y0,y1) = extrema(posx), extrema(posy)
    set_default_graphic_size(x1-x0+3mm, y1-y0+3mm)
    compose(context(), hexagon(posx, posy, 1mm), fill(RGB.(0.0, s.membrane[:,1] ./ 30, s.membrane[:,2] ./ 30))    )
end

function Plotter(M::ModelN)
    stats(_, s) = display(plot(M, s))
end



T = 20000.0
Tmeas = 10.0
Nsave = 10
L = 40

#for reproducibility
seed = 22
ran_ng = Random.Xoshiro(seed)
M,s = build_model_state(L,rng=ran_ng)
p = ProgressShower(T)
#m = Measurer(M; name="test_example_1", Nsave)
pl = TimeFilter(Plotter(M); times=Tmeas:200:T)
stats = TimeFilter(p; times=Tmeas:Tmeas:T)

@time run_RDN!(s, M, T; stats, rng=ran_ng) 