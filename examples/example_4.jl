using MembraneRD, ExponentialQueues
using Compose, MembraneRD.Colors
import MembraneRD: hexagon
using MembraneRD.Graphs


Base.@kwdef struct ModelN{G, S, T, L1, L2, L3, L4}
    species::S
    g::G
    rea::L1
    att::L2
    det::L3
    dif::L4
    rho_0::T
end

nspecies(M::ModelN) = length(M.species)

Base.length(M::ModelN) = nv(M.g)

struct StateN
	membrane::Matrix{Int}
	cytosol::Vector{Int}
end

function StateN(M::ModelN, totmembrane::Vector, cytosol::Vector; rng)
    Nspecies = nspecies(M)
    membrane = zeros(nv(M.g), Nspecies)
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

    N = length(M)
    Nspecies = nspecies(M)
    Qn = [ExponentialQueue(N) for _ in M.species]
    Qcat = [ExponentialQueue(N) for _ in M.rea]
    Qatt = [ExponentialQueue(N)*0.0 for _ in M.att]
    
    function update(i::Int)
        for ((e,s,_,_,km),q) in zip(M.rea, Qcat)
            q[i] = state.membrane[i,e] * state.membrane[i,s] / (state.membrane[i,s] + km)
        end
        for ((m,m1,_),q) in zip(M.att, Qatt)
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
            ((:dif,m) => Qn[m] * d for (m,d) in M.dif)...,
            ((:att,m) => q*ka for ((m,_,ka),q) in zip(M.att, Qatt))...,
            ((:det,m) => Qn[m] * kd for (m,kd) in M.det)...,
            ((:cat,(s,p)) => q*kc for ((_,s,p,kc,_),q) in zip(M.rea, Qcat))...,
        )

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



function build_model(L)
    g,posx,posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    unit_length = 1.0
    theta = 0.2126944621086619 #Stot/V
    Ac = unit_length^2 #area associated to each cell, set unit lengthscale
    Stot = N*Ac #area of each cell is assumed to be 1, setting unit lengthscale
    V = Stot/theta
    d_timescale = 0.01 #this sets unit timescale
    dA,dB,dEA,dEB = Iterators.repeated(d_timescale)
    #rates in the theory (do not correspond excatly to those used in the simulations), some slight dimensional changes are needed
    kAa_th = 1.0; kAd_th = 1.0*kAa_th; kAc_th = 1.0
    kBa_th = 1.0; kBd_th = 1.0*kBa_th; kBc_th = 1.3*kAc_th
    KMMA_th = 1.0; KMMB_th = 1.0
    #rates to implement
    kAc = kAc_th; kBc = kBc_th; kAa = kAa_th/V
    kAd = kAd_th; kBa = kBa_th/V; kBd = kBd_th
    KMMA = KMMA_th*Ac; KMMB = KMMB_th*Ac

    species = ("A","B","EA","EB")
    A,B,EA,EB = Iterators.countfrom(1)
    rea = ((EA,B,A,kAc,KMMA), (EB,A,B,kBc,KMMB))
    att = ((EA,A,kAa), (EB,B,kBa))
    det = ((EA,kAd),(EB,kBd))
    dif = ((A,dA),(B,dB),(EA,dEA),(EB,dEB))

    M = ModelN(; species, g, rea, att, det, dif, rho_0 = 0.0)

    totmol = N * 10
    totA, totB = floor(Int, totmol/2), floor(Int, totmol/2)
    totEA, totEB = floor(Int, 0.1*N), floor(Int, 0.1*N)
    #memEA = floor(Int, totEA*(theta/(kAd_th/kAa_th))*(totA/Stot)/(1+theta/(kAd_th/kAa_th)*totA/Stot))
    #memEB = floor(Int, totEB*(theta/(kBd_th/kBa_th))*(totB/Stot)/(1+theta/(kBd_th/kBa_th)*totB/Stot))
    memEA, memEB = 0.0, 0.0
    cytoEA, cytoEB = totEA - memEA, totEB - memEB
    (; M, mem = [totA, totB, memEA, memEB], cyto = [0.0, 0.0, cytoEA, cytoEB], posx, posy)
end


function build_model3(L)
    g,posx,posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    species = ("A","B","EA","EB","C","EC")
    A,B,C,EA,EB,EC = Iterators.countfrom(1)
    rea = ((EA,B,A,1.0,1.0), (EB,A,B,1.0,1.0),(EC,A,C,1.0,1.0),(EA,C,A,1.0,1.0))
    att = ((EA,A,1.0), (EB,B,1.0), (EC,C,1.0))
    det = ((EA,1.0),(EB,1.0),(EC,1.0))
    dif = ((A,1.0),(B,1.0),(EA,1.0),(EB,1.0),(C,1.0),(EC,1.0))
    M = ModelN(; species, g, rea, att, det, dif, rho_0 = 0.0)
    (; M, mem = [10*N, 10*N, 10*N, 0, 0, 0], 
        cyto = [fill(0,3); fill(floor(Int, 0.5*N),3)], 
        posx, posy)
end





T = 2000.0
Tmeas = 10.0
Nsave = 10
L = 70

#for reproducibility
seed = 22
rng = Random.Xoshiro(seed)
(; M, mem, cyto, posx, posy) = build_model3(L)
s = StateN(M.g, mem, cyto; rng)

p = ProgressShower(T)
#m = Measurer(M; name="test_example_1", Nsave)
pl = TimeFilter(PlotterN(posx, posy); times=Tmeas:50:T)
stats = TimeFilter(p, pl; times=Tmeas:Tmeas:T)

@time run_RDN!(s, M, T; stats, rng) 