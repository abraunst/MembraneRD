using MembraneRD

#here define two examples of measures
	
using ProgressMeter, JLD, Random

struct Measure
	time::Float64
    phi_av::Float64
    binder_cumulant::Float64
    cytoEA::Float64
    cytoEB::Float64
    rho::Float64
end

#make measures here
function Measurer(M::Model; name, Nsave)
	measures = Measure[]
    Nc = length(M)
    next::Int = 1
    function take_measure(t, s)
        #calculate phi at each cell (normalized, goes from -1, to 1)
        ϕ(i) = (s.nB[i] - s.nA[i])/(s.nA[i] + s.nB[i] + 1e-10)
        totA, totB = sum(s.nA), sum(s.nB)
        #measure phi
        ϕav = (totB - totA) / (totA + totB)
        #in some denominators I added + 10^-10 in order to avoid NaNs
        m2 = sum((ϕ(i)-ϕav)^2 for i in 1:Nc) / Nc
        m4 = sum((ϕ(i)-ϕav)^4 for i in 1:Nc) / Nc
        #measure Binder cumulant
        bc = 1 - (m4 / (m2 ^ 2 + 1e-10)) / 3			
        #the following is to check for convergence of rho to 1
        rho = M.rho_0 * ((M.kAd / M.kAa) + totA) / ((M.kBd / M.kBa) + totB)
        push!(measures, Measure(t,ϕav, bc, s.cytoEA[], s.cytoEB[], rho))
        if next % Nsave ==0
            println("T = $t and <ϕ>/c = $(ϕav)")
            println("Binder cumulant: $bc")
            save("$name/config_T=$(round(t,digits=3)).jld",compress=true, "state", s)
            save("$name/measures.jld",compress=true, "measures", measures)
        end
        next += 1
    end
end


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
	kAs=0.0
	kBs=0.0

    M = Model(g, posx, posy, dA, dB, dEA, dEB, kAc, kBc, kAa, kAd, kBa, kBd,kAs,kBs, KMMA, KMMB, 0.0)

    totmol = N * 10
    totA, totB = floor(Int, totmol / 2), floor(Int, totmol / 2)
    totEA, totEB = floor(Int, 0.1*N), floor(Int, 0.1*N)
    memEA = floor(Int, totEA*(theta/(kAd_th/kAa_th))*(totA/Stot)/(1+((theta/(kAd_th/kAa_th))*(totA/Stot))))
    memEB = floor(Int, totEB*(theta/(kBd_th/kBa_th))*(totB/Stot)/(1+((theta/(kBd_th/kBa_th))*(totB/Stot))))
    cytoEA, cytoEB = totEA - memEA, totEB - memEB
    s = State(M, totA, totB, memEA, memEB, cytoEA, cytoEB; rng)
    M,s
end

T = 20000.0
Tmeas = 10.0
Nsave = 10
L = 100

#for reproducibility
seed = 22
ran_ng = Random.Xoshiro(seed)
M,s = build_model_state(L,rng=ran_ng)
p = ProgressShower(T)
m = Measurer(M; name="test_example_1", Nsave)
pl = TimeFilter(MembraneRD.Plotter(M); times=Tmeas:200:T)
stats = TimeFilter(m, p, pl; times=Tmeas:Tmeas:T)

@time run_RD!(s, M, T; stats, rng=ran_ng) 