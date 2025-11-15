function run_RD!(state::State, M::Model, T; 
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
    L,C = LinearIndices((Nspecies,Nspecies)), CartesianIndices((Nspecies,Nspecies))

    Q = NestedQueue(
            ((:dif,m) => Qn[m] * d for (m,d) in M.dif)...,
            ((:att,m) => q*ka for ((m,_,ka),q) in zip(M.att, Qatt))...,
            ((:det,m) => Qn[m] * kd for (m,kd) in M.det)...,
            ((:cat,L[s,p]) => q*kc for ((_,s,p,kc,_),q) in zip(M.rea, Qcat))...,
        )

    println("starting simulation, $(length(Q)) events in the queue")

    t::Float64 = 0.0
    while !isempty(Q)
        ((ev, m), i), dt = peek(Q; rng)
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
            s, p = Tuple(C[m])
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