using LinearAlgebra, SparseArrays, Graphs

Open(L) = spdiagm(1=>trues(L-1))
Closed(L) = spdiagm(1=>trues(L-1), -L+1 => trues(1))
sym(x) = x .| x'
rect(L; boundary=Closed) = sym(boundary(L))
rect(L1, L...; boundary=Closed) = kron(rect(L1; boundary),I(prod(L))) .| kron(I(L1), rect(L...; boundary))
dashed(L; boundary=Closed) = (A = spdiagm(1 => isodd.(1:L-1), -1 => iseven.(1:L-1)) |> dropzeros!; boundary === Closed && (A[1,end] = true); A) 
hexa(L1, L2; boundary=Closed) = rect(L1,L2; boundary) .| sym(kron(dashed(L1; boundary),boundary(L2)))
gen_square_lattice(L) = SimpleGraph(rect(L, L)), mod1.(1:L^2, L), fld1.(1:L^2, L)
function gen_hex_lattice(L; boundary=Closed)
    g = SimpleGraph(hexa(L, L; boundary))
    x, y = mod1.(1:L^2, L), fld1.(1:L^2, L)
    g, x .+ 0.5 .* isodd.(y), y
end
