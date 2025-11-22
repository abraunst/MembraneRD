using MembraneRD, Random, Colors, VideoIO, Compose, Images, Printf
import Cairo, Fontconfig


function build_model3(L)
    G, posx, posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    species = @species A B C EBA ECA EAB EAC
    cat = (@reaction (1.0,1.0) EBA + B => EBA + A),
        (@reaction (1.0,1.0) EAB + A => EAB + B),
        (@reaction (1.0,1.0) EAC + A => EAC + C),
        (@reaction (1.0,1.0) ECA + C => ECA + A)
    att = (EBA,A,1.0), (EAB,B,1.0), (ECA,A,1.0), (EAC,C,1.0)
    det = (EBA,1.0), (EAB,1.0), (ECA,1.0), (EAC,1.0)
    dif = (A,0.1), (B,0.1), (C,0.1), 
        (EBA,0.1), (ECA,0.1), (EAB,0.1), (EAC,0.1)
    M = Model(; species, G, cat, att, det, dif, rho_0 = 0.0)
    (; M, mem = [2N,10N,10N,0,0,0,0], 
        cyto = floor.(Int,[0,0,0,0.5N,0.5N,0.5N,0.5N]), 
        posx, posy)
end

T = 2000.0
Nsave = 10
L = 50

#for reproducibility
seed = 22
rng = Random.Xoshiro(seed)
(; M, mem, cyto, posx, posy) = build_model3(L)
s = State(M, mem, cyto; rng)


#m = Measurer(M; name="test_example_1", Nsave)

times = 0:100:T
saver = Pusher(Tuple{Float64,State})
colors = [RGB(m==1,m==2,m==3)/30 for m in 1:nspecies(M)]
plotter = Plotter(posx, posy; colors)
displayer = StopWatchFilter(display ∘ plotter; seconds=5.0)


stats = TimeFilter(ProgressShower(T), saver, displayer; times)
@profview run_RD!(s, M, T; stats, rng)
savevideo("video2.mp4", saver.stack, plotter)


