using LinearAlgebra, SparseArrays, Graphs, Compose, Colors

Open(L) = spdiagm(1=>trues(L-1))
Closed(L) = spdiagm(1=>trues(L-1), -L+1 => trues(1))
sym(x) = x .| x'
rect(L; boundary=Closed) = sym(boundary(L))
rect(L1, L...; boundary=Closed) = kron(rect(L1; boundary),I(prod(L))) .| kron(I(L1), rect(L...; boundary))
hexa(L1, L2; boundary=Closed) = rect(L1,L2; boundary) .| sym(kron(boundary(L1),boundary(L2)))
gen_square_lattice(L) = SimpleGraph(rect(L, L)), mod1.(1:L^2, L), fld1.(1:L^2, L)
function gen_hex_lattice(L)
    g = SimpleGraph(hexa(L, L))
    x, y = mod1.(1:L^2, L), fld1.(1:L^2, L)
    g, mod1.((y .- 1)/2 + x .- 1, L), (y .- 1) .* sqrt(3) ./ 2 
end


