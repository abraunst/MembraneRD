using ColorVectorSpace

function ProgressShower(T)
    p = Progress(100)
    stat(t, _) = update!(p, floor(Int, 100*t/T))
end

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

function TicFilter(callbacks...; seconds=1)
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

function Plotter(f, posx, posy; colors, Δ=mm)
    function plot(t::Float64, s::State)
        X, Y = Δ .* posx, Δ .* posy
        (x0,x1),(y0,y1) = extrema(X), extrema(Y)
        set_default_graphic_size(x1-x0+3Δ, y1-y0+3Δ)
        compose(context(), 
            hexagon(X, Y, 1Δ),
                fill(s.membrane * colors)
        ) |> f
    end
end