function run_RD!(s::State, M::Model, T; 
        stats = (x...)->nothing, 
        rng = Random.default_rng())

    QA,QB,QEA,QEB,QcatA,QcatB = (ExponentialQueue(length(M)) for _ in 1:6)
    QattEA = QA * 0.0
    QattEB = QB * 0.0

    function update(i)
        QA[i] = s.nA[i]
        QB[i] = s.nB[i]
        QcatB[i] = s.nEB[i] * s.nA[i] / (s.nA[i] + M.KMM)
        QcatA[i] = s.nEA[i] * s.nB[i] / (s.nB[i] + M.KMM)
        QEA[i] = s.nEA[i]
        QEB[i] = s.nEB[i]
        QattEA.f[] = s.cytoEA[] * M.kAa
        QattEB.f[] = s.cytoEB[] * M.kBa
    end

    foreach(update, 1:length(M))

    #arrival is chosen uniformly between its neighbours
    rand_neighbor(i) = rand(rng, neighbors(M.g, i))

    Q = NestedQueue(
            :difA => QA * M.dA,
            :difB => QB * M.dB,
            :catA => QcatA * M.kAc,
            :catB => QcatB * M.kBc,
            :attEA => QattEA,
            :attEB => QattEB,
            :difEA => QEA * M.dEA,
            :difEB => QEB * M.dEB,
            :detEA => QEA * M.kAd,
            :detEB => QEB * M.kBd
        )

    println("starting simulation, $(length(Q)) events in the queue")

    t::Float64 = 0.0
    while !isempty(Q)
        (ev, i), dt = peek(Q; rng)
        t += dt
        t > T && break
        stats(t, s)
        @inbounds if ev === :difA #diffusion of specie A
            j = rand_neighbor(i)
            s.nA[i] -= 1
            s.nA[j] += 1
            update(i)
            update(j)
        elseif ev === :difB #diffusion of specie B
            j = rand_neighbor(i)
            s.nB[i] -= 1
            s.nB[j] += 1
            update(i)
            update(j)
        elseif ev === :catA #B+EA->A+EA
            s.nB[i] -= 1
            s.nA[i] += 1
            update(i)
        elseif ev === :catB #A+EB->B+EB
            s.nA[i] -= 1
            s.nB[i] += 1
            update(i)
        elseif ev === :attEA #attachment of EA from cytosol
            s.cytoEA[] -= 1
            s.nEA[i] += 1
            update(i)
        elseif ev === :attEB #attachment of EB from cytosol
            s.cytoEB[] -= 1
            s.nEB[i] += 1
            update(i)
        elseif ev === :detEA #detachment of EA
            s.nEA[i] -= 1
            s.cytoEA[] += 1
            update(i)
        elseif ev === :detEB #detachment of EB
            s.nEB[i] -= 1
            s.cytoEB[] += 1
            update(i)
        elseif ev === :difEA #diffusion of EA
            j = rand_neighbor(i)
            s.nEA[i] -= 1
            s.nEA[j] += 1
            update(i)
            update(j)
        elseif ev === :difEB #diffusion of EB
            j = rand_neighbor(i)
            s.nEB[i] -= 1
            s.nEB[j] += 1
            update(i)
            update(j)
        end
    end
end
