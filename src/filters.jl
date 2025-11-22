using ColorVectorSpace, Images, ProgressMeter, Colors, Compose
import VideoIO

"""
`ProgressShower(T; kwds...)`

A filter to show progress

`T`: final time
"""
function ProgressShower(T; kwds...)
    p = Progress(100; kwds...)
    stat(t, _) = update!(p, floor(Int, 100*t/T))
end

"""
`TimeFilter(callbacks...; times)`

A filter that calls each filter in `callbacks` once per each time in `times`. 
This is intended to run in the innest cycle of `run_RD` and it's very fast
"""
function TimeFilter(callbacks...; times)
    next = Ref(1)
    function stats(t, s)
        while next[] <= lastindex(times) && times[next[]] <= t
            for cb in callbacks
                cb(times[next[]], s)
            end
            next[] += 1
        end
    end
end

"""
`StopWatchFilter(callbacks...; seconds=1)`

A filter that calls `callbacks` with at least `seconds` wall-clock 
delay between calls (intended e.g. for iterative plotting)
"""
function StopWatchFilter(callbacks...; seconds=1)
    last = time()
    function stats(t, s)
        next = time()
        if next - last > seconds
            last = next
            for cb in callbacks
                cb(t, s)
            end
        end
    end
end

function hexagon(x, y, r)
    polygon([[(x[i]+r*cos(π/3*(k+3/2)),
               y[i]+r*sin(π/3*(k+3/2))) for k in 0:6] 
                    for i in eachindex(x,y)])
end

"""
`composed(s::State; posx, posy, colors, Δ)`

Generates a `Compose` image of an hexagonal grid from `s::State` 
"""
function composed(s::State; posx, posy, colors, Δ)
    X, Y = Δ .* posx, Δ .* posy
    (x0,x1),(y0,y1) = extrema(X), extrema(Y)
    set_default_graphic_size(x1-x0+3Δ, y1-y0+3Δ)
    compose(context(), 
        hexagon(X, Y, 1Δ),
        fill(s.membrane * colors)
    ) 
end

"""
`Plotter(posx, posy; colors, Δ=mm)`

A filter to generate a `Compose` image from each state `s`. Use e.g. the
`display ∘ Plotter(posx, posy; colors)` filter for iterative display, etc.
"""
function Plotter(posx, posy; colors, Δ=mm)
    f(_, s::State) = composed(s; posx, posy, colors, Δ)
end

"""
`raster(c::Compose.Context) -> Matrix{<:RGB}`

Converts a `Context` object into a raster image.
"""
function raster(c)
    io = IOBuffer()
    draw(PNG(io; emit_on_finish=false), c)
    A = Images.load(io)
    m, n = 2 .* (size(A) .÷ 2)
    #crop to even dimensions and N0f8 colorspace
    RGB{N0f8}.(@view A[1:m,1:n])
end

"""
Pushes something into a vector
"""
function Pusher(::Type{T} = Any) where T
    stack = T[]
    pusher(x...) = push!(stack, x)
end

"""
`savevideo(fn, stack, plotter; kwd)`

Saves a .MP4 file generated from a stack of either `State`s or `(::Float64,::State)` tuples 

`fn`: filename
`stack::AbstractVector`: stack of images
`plotter`: converts a `t,s` into an image (e.g. a `Plotter` works)
"""
function savevideo(fn, stack::AbstractVector{Tuple{Float64,State}}, plotter; framerate=10, encoder_options = (crf=23, preset="medium"))
    savevideo(fn, [plotter(t, s) for (t,s) in stack]; framerate, encoder_options)
end

function savevideo(fn, stack::AbstractVector{State}, plotter; framerate=10, encoder_options = (crf=23, preset="medium"))
    savevideo(fn, [plotter(0.0, s) for s in stack]; framerate, encoder_options)
end

"""
`savevideo(fn, stack, plotter; kwd)`

Saves a .MP4 file generated from a stack of `RGB` images or `Compose.Context` images

`fn`: filename
`stack::AbstractVector`: stack of images
"""
function savevideo(fn, stack::AbstractVector{<:Compose.Context}; framerate=10, encoder_options = (crf=23, preset="medium"))
    savevideo(fn, raster.(stack); framerate, encoder_options)
end

function savevideo(fn, stack::AbstractVector{<:AbstractMatrix}; framerate=10, encoder_options = (crf=23, preset="medium"))
    VideoIO.save(fn, stack; framerate, encoder_options)
end

"""
`OpenVideo(f, filename, s::State; plotter, kwd...)`

This function is intended for iterative video generation, without saving
intermediate raster images and thus saving memory/disk space

use e.g. as follows:

```julia
OpenVideo("video.mp4", s; plotter) do framesaver
    stats = TimeFilter(framesaver; times=Tmeas:Tmeas:T)
    run_RD!(s, M, T; stats, rng)
end
```

where `plotter` takes `t,s` and produces a bitmap matrix.
"""
function OpenVideo(f, fn, s::State; plotter, framerate=10, encoder_options = (crf=23, preset="medium"))
    VideoIO.open_video_out(w->f((t,s)->write(w,plotter(t,s))), fn, plotter(0.,s); framerate, encoder_options)
end