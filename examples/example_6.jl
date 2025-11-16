using MembraneRD, Random, Colors, VideoIO, Compose, ImageMagick, Printf
import Cairo, Fontconfig

function build_model3(L)
    g, posx, posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    species =  1:9
    A,B,C,EAB,EAC,EBA,EBC,ECA,ECB = 1:9
    rea = (EBA,B,A,1.0,1.0), (EAB,A,B,1.0,1.0),(EAC,A,C,1.0,1.0),(ECA,C,A,1.0,1.0),(EBC,B,C,1.0,1.0),(ECB,C,B,1.0,1.0)
    att = (EBA,A,1.0), (EBC,C,1.0), (EAB,B,1.0), (EAC,C,1.0), (ECA,A,1.0), (ECB,B,1.0)
    det = (EBA,1.0), (EBC,1.0), (EAB,1.0), (EAC,1.0), (ECA,1.0), (ECB,1.0)
    dif = Tuple((i,1.0) for i in eachindex(species))
    M = Model(; species, g, rea, att, det, dif, rho_0 = 0.0)
    (; M, mem = [fill(10*N, 3); fill(0, 6)], 
        cyto = [fill(0,3); fill(floor(Int, 0.5*N),6)], 
        posx, posy)
end

T = 6000.0
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
        draw(PNG((@sprintf "%s_%06d.png" name i)), s)
        i += 1
        nothing
    end
end


pl = TicFilter(Plotter(PNGSaver("images/image"), posx, posy; colors=[RGB(m==1,m==2,m==3)/30 for m in 1:nspecies(M)]); seconds=5.0)
stats = TimeFilter(p, pl; times=Tmeas:Tmeas:T)

@time run_RD!(s, M, T; stats, rng)


