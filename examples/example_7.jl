using MembraneRD, Random, Colors, VideoIO, Compose, Images, Printf
import Cairo, Fontconfig


function build_model3(L)
    g, posx, posy = MembraneRD.gen_hex_lattice(L)
    N = length(posx)
    species = @species A B C EBA ECA EAB EAC
    cat = (@reaction (1.0,1.0) EBA + B => EBA + A),
        (@reaction (2.0,1.0) EAB + A => EAB + B),
        (@reaction (2.0,1.0) EAC + A => EAC + C),
        (@reaction (1.0,1.0) ECA + C => ECA + A)
    att = (EBA,A,1.0), (EAB,B,1.0), (EAC,C,1.0), (ECA,A,1.0)
    det = (EBA,1.0), (EAB,1.0), (EAC,1.0), (ECA,1.0)
    dif = (A,0.2), (B,0.1), (C,0.1), (EBA,0.1), (ECA,0.1), (EAB,0.1), (EAC,0.02)
    M = Model(; species, g, cat, att, det, dif, rho_0 = 0.0)
    (; M, mem = [2N,10N,10N,0,0,0,0], 
        cyto = floor.(Int,[0,0,0,0.5N,0.5N,0.5N,0.5N]), 
        posx, posy)
end

T = 10000.0
Tmeas = 10.0
Nsave = 10
L = 150

#for reproducibility
seed = 22
rng = Random.Xoshiro(seed)
(; M, mem, cyto, posx, posy) = build_model3(L)
s = State(M, mem, cyto; rng)

p = ProgressShower(T)
#m = Measurer(M; name="test_example_1", Nsave)

function PNGSaver()
    stack = Matrix{RGB{N0f8}}[]
    function f(s)
        io = IOBuffer()
        draw(PNG(io; emit_on_finish=false), s)
        A = Images.load(io)
        m, n = 2 .* (size(A) .÷ 2)
        push!(stack, [RGB{N0f8}(A[i,j]) for i ∈ 1:m, j ∈ 1:n])
        nothing
    end
end

function savevideo(saver, filename)
    encoder_options = (crf=23, preset="medium")
    VideoIO.save(filename, saver.stack; framerate=10, encoder_options)
end

saver = PNGSaver()
pl = TicFilter(Plotter(saver, posx, posy; colors=[RGB(m==1,m==2,m==3)/30 for m in 1:nspecies(M)]); seconds=5.0)
pl2 = TicFilter(Plotter(display, posx, posy; colors=[RGB(m==1,m==2,m==3)/30 for m in 1:nspecies(M)]); seconds=5.0)
stats = TimeFilter(p, pl, pl2; times=Tmeas:Tmeas:T)

@time run_RD!(s, M, T; stats, rng)
savevideo(saver, "prova.mp4")




