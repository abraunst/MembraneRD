using ProgressMeter

struct Measure
    phi_av::Float64
    binder_cumulant::Float64
    cytoEA::Float64
    cytoEB::Float64
    rho::Float64
end

#make measures here
function Measurer(M::Model; times)
    isempty(times) && return
	measures = Measure[]
    Nc = length(M)
	next = 1
    function take_measure(t, s)
        #calculate phi at each cell (normalized, goes from -1, to 1)
        ϕ(i) = (s.nB[i] - s.nA[i])/(s.nA[i] + s.nB[i] + 1e-10)
        while next <= lastindex(times) && times[next] <= t 
            totA, totB = sum(s.nA), sum(s.nB)
            #measure phi
            ϕav = (totB - totA) / (totA + totB)
            println("T = $(times[next]) and <ϕ>/c = $(ϕav)")
            #in some denominators I added + 10^-10 in order to avoid NaNs
            m2 = sum((ϕ(i)-ϕav)^2 for i in 1:Nc) / Nc
            m4 = sum((ϕ(i)-ϕav)^4 for i in 1:Nc) / Nc
            #measure Binder cumulant
            bc = 1 - (m4 / (m2 ^ 2 + 1e-10)) / 3
            println("Binder cumulant: $bc")			
            #the following is to check for convergence of rho to 1
            rho = M.rho_0 * ((M.kAd / M.kAa) + totA) / ((M.kBd / M.kBa) + totB)
            push!(measures, Measure(ϕav, bc, s.cytoEA[], s.cytoEB[], rho))
            next += 1
        end
    end
end

function ProgressMeasurer(T)
    p = Progress(100)
    oldval = 0
    function stat(t, _)
        newval = floor(Int, 100*t/T)
        if newval > oldval
            update!(p, newval)
            oldval = newval
        end
    end
end


function PlotMeasurer(M; times)
    next::Int = 1
    function stats(t, s)
        if t > times[next]
            display(plot(M,s))
            display(times[next])
            next += 1
        end
    end
end