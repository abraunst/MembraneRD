using LinearAlgebra, SparseArrays, Graphs

Open(L) = spdiagm(1=>trues(L-1))
Closed(L) = spdiagm(1=>trues(L-1), -L+1 => trues(1))
sym(x) = x .| x'
rect(L; boundary=Closed) = sym(boundary(L))
rect(L1, L...; boundary=Closed) = kron(rect(L1; boundary),I(prod(L))) .| kron(I(L1), rect(L...; boundary))
hexa(L1, L2; boundary=Closed) = rect(L1,L2; boundary) .| sym(kron(boundary(L1),boundary(L2)))
gen_square_lattice(L) = SimpleGraph(rect(L, L)), mod1.(1:L^2, L), fld1.(1:L^2, L)
gen_hex_lattice(L) = SimpleGraph(hexa(L, L)), mod1.(1:L^2, L) .- 0.5 .* fld1.(1:L^2, L), 1.0.*fld1.(1:L^2, L)
