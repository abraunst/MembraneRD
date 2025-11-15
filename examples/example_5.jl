using MembraneRD, Random, Colors, Compose
import Cairo, Fontconfig

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

    M = Model(; species, g, rea, att, det, dif, rho_0 = 0.0)

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
    g, posx, posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    species = ("A","B","EA","EB","C","EC")
    A,B,C,EA,EB,EC = Iterators.countfrom(1)
    rea = ((EA,B,A,1.0,1.0), (EB,A,B,1.0,1.0),(EC,A,C,1.0,1.0),(EA,C,A,1.0,1.0))
    att = ((EA,A,1.0), (EB,B,1.0), (EC,C,1.0))
    det = ((EA,1.0),(EB,1.0),(EC,1.0))
    dif = ((A,1.0),(B,1.0),(EA,1.0),(EB,1.0),(C,1.0),(EC,1.0))
    M = Model(; species, g, rea, att, det, dif, rho_0 = 0.0)
    (; M, mem = [10*N, 10*N, 10*N, 0, 0, 0], 
        cyto = [fill(0,3); fill(floor(Int, 0.5*N),3)], 
        posx, posy)
end

T = 5000.0
Tmeas = 10.0
Nsave = 10
L = 200

#for reproducibility
seed = 22
rng = Random.Xoshiro(seed)
(; M, mem, cyto, posx, posy) = build_model3(L)
s = State(M, mem, cyto; rng)

p = ProgressShower(T)
#m = Measurer(M; name="test_example_1", Nsave)

function PNGSaver(name, i::Int = 1)
    function f(s)
        draw(PNG(@printf "%s_%06d.png" name i), s)
        i += 1
        nothing
    end
end


#imagestacker = ImgStacker(backend=()->PNG(Cairo.CairoImageSurface(30, 30,1)))
pl = TicFilter(Plotter(PNGSaver("images/image"), posx, posy; colors=[RGB(m==1,m==2,m==3)/30 for m in 1:nspecies(M)]); seconds=5.0)
stats = TimeFilter(p, pl; times=Tmeas:Tmeas:T)

@time run_RD!(s, M, T; stats, rng)


