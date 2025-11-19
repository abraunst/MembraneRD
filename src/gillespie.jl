@enum Event evdif evatt evdet evcat evrea


function run_RD!(state::State, M::Model, T; 
        stats = (_, _)->nothing, 
        rng = Random.default_rng())

    N, Nspecies  = nsites(M), nspecies(M)
    Qn = [StaticExponentialQueue(N) for _ in M.species]
    Qcat = [StaticExponentialQueue(N) for _ in M.cat]
    Qrea = [StaticExponentialQueue(N) for _ in M.rea]
    Qatt = [StaticExponentialQueue(N)*0.0 for _ in M.att]
    
    function update(i::Int)
        for ((e,s,_,_,km),q) in zip(M.cat, Qcat)
            q[i] = state.membrane[i,e]  / (1 + km/state.membrane[i,s])
        end
        for ((s,_,_), q) in zip(M.rea, Qrea)
            q[i] = prod(state.membrane[i,m] for m in s)
        end
        for ((m,m1,_),q) in zip(M.att, Qatt)
            q.q[i] = state.membrane[i,m1]
            q.f[] = state.cytosol[m]
        end
        for m in 1:Nspecies
            Qn[m][i] = state.membrane[i,m]
        end
    end

    foreach(update, 1:N)

    Q = NestedQueue(
            ((evdif,m) => Qn[m] * d for (m,d) in M.dif)...,
            ((evatt,m) => q*ka for ((m,_,ka),q) in zip(M.att, Qatt))...,
            ((evdet,m) => Qn[m] * kd for (m,kd) in M.det)...,
            ((evcat,r) => q*kc for (r,(_,_,_,kc,_),q) in zip(Iterators.countfrom(1),M.cat, Qcat))...,
            ((evrea,r) => q*k for ((_,_,k),q) in zip(M.rea, Qrea))...,
        )

    println("starting simulation, $(length(Q)) events in the queue")

    t::Float64 = 0.0
    while !isempty(Q)
        ((ev::Event,m::Int),i::Int),dt::Float64 = peek(Q; rng)
        t += dt
        t > T && break # reached end time for simulation
        stats(t, state)
        @inbounds if ev == evdif # diffusion
            #arrival is chosen uniformly between its neighbours
            j = rand(rng, neighbors(M.G, i))
            state.membrane[i,m] -= 1
            state.membrane[j,m] += 1
            update(i)
            update(j)
        elseif ev == evcat # catalytic reaction
            _, s, p = M.cat[m]
            state.membrane[i,s] -= 1
            state.membrane[i,p] += 1
            update(i)
        elseif ev == evatt #attachment to membrane
            state.cytosol[m] -= 1
            state.membrane[i,m] += 1
            update(i)
        elseif ev == evdet # detachment from membrane
            state.membrane[i,m] -= 1
            state.cytosol[m] += 1
            update(i)
        else # ev == evrea # reaction
            s, p = M.rea[m]
            for m in s
                state.membrane[i,m] -= 1
            end
            for m in p
                state.membrane[i,m] += 1
            end
            update(i)
        end
    end
	stats(T, state)
end