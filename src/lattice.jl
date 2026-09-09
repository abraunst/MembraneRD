using LinearAlgebra, SparseArrays, Graphs

Open(L) = spdiagm(1=>trues(L-1))
Closed(L) = spdiagm(1=>trues(L-1), -L+1 => trues(1))
sym(x) = x .| x'
rect(L; boundary=Closed) = sym(boundary(L))
rect(L1, L...; boundary=Closed) = kron(rect(L1; boundary),I(prod(L))) .| kron(I(L1), rect(L...; boundary))
dashed(L; boundary=Closed) = spdiagm(1 => isodd.(1:L-1), -1 => iseven.(1:L-1), L-1 => [boundary == Closed]) |> dropzeros!
hexa(L1, L2; boundary=Closed) = rect(L1,L2; boundary) .| sym(kron(dashed(L1; boundary),boundary(L2)))
gen_square_lattice(L1, L2 = L1) = (; g = SimpleGraph(rect(L1, L2)), layout=(cellshape=:square, posx = mod1.(1:L1*L2, L1), posy = fld1.(1:L1*L2, L1)))
function gen_hex_lattice(L1, L2 = L1; boundary=Closed)
    g = SimpleGraph(hexa(L1, L2; boundary))
    x, y = mod1.(1:L1*L2, L1), fld1.(1:L1*L2, L1)
    (; g = g, layout = (cellshape=:hexagon, posx = (x .+ 0.5 .* isodd.(y)) .* 2cos(π/6), posy = y .*(1 + cos(π/3))))
end

