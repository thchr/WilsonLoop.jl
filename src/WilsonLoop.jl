module WilsonLoop
using StaticArrays
using LinearAlgebra

export wilsonloop, berryphase, wannierhamiltonian, wannierstate, wilsonspectrum, spectrum

include("Hamiltonians.jl") # defines Hamiltonians module with example Hamiltonians

""" 
    wilsonloop(statefunc, kpath, coords=nothing)

Calculates the Wilson loop W around a path `kpath`, for a function `statefunc` that returns 
a selection (i.e. a multiplet) of (normalized!) states when evaluated at a (vectorial)
**k**-point (see e.g. [`pickeigvecs!`](@ref)); the states are assumed to associate with a
Hermitian Hamiltonian (though it does not need to be explicitly parametrizable).

## Arguments
- `statefunc`: a single-argument function that returns an set of `Ns` states (as a complex 
  `Array` of size dim(H) × `Ns`) evaluable at a (vectorial) **k**-point.
- `kpath`: a vector of k-points, specifying the path of the Wilson loop; it is expected to 
  span a range of 2π, i.e. is specified in normalized units.
- `coords`: a vector of coordinates associated with the elements of the Hamiltonian or 
  eigenstates; generally, a vector of vectors (dim H × dim R). If `nothing`, all coordinates
  are set to the origin.
"""
function wilsonloop(statefunc::Function, kpath::AbstractVector, coords::Union{AbstractVector,Nothing}=nothing)
    Nk = length(kpath)
    
    uk1 = statefunc(kpath[1]) # pickeigvecs!(H(kpath[1]), states); 
    ukjm1 = uk1 # don't have to worry about array aliasing (i.e., no need to copy(...)), since we always do full overwrites, not elementwise  
    Ns = size(uk1,2)
    M = Matrix{eltype(eltype(uk1))}(I, Ns, Ns)
    Mijwork = similar(M)

    # All overlaps <u(k_{j})|u(k_{j-1})> for j = 2, 3, .., N-1
    for j = 2:Nk-1
        ukj = statefunc(kpath[j]) # pickeigvecs!(H(kpath[j]), states)
        M = overlap!(Mijwork, ukj, ukjm1)*M
        ukjm1 = ukj
    end
    
    # Last overlap <u(k_N)|u(k_{N-1})>; u(k_N) is obtained from u(k_1), since k_1 = k_N mod G, cf. periodicity
    if !isnothing(coords) # if coords is nothing, we assume it means all coordinates equal the origin, i.e. phase term is unity
        G = kpath[end] .- kpath[1]
        phase = exp.(-1im*map(r -> sum(r.*G), coords)) # the map(...) is just a fancy way of writing dot(G,r) when r is a vector of coordinates
        for uk1_n in eachcol(uk1) # create u(k_N) from u(k_1) by multiplying appropriate phase factor; overwrite into uk1
            uk1_n .*= phase
        end
    end
    M = overlap!(Mijwork, uk1, ukjm1)*M
    
    return M
end


""" 
    wilsonloop((H::Function, states), kpath, coords=nothing)

Return the Wilson loop around a **k**-path `kpath`, for `states` of a Hamiltonian function
`H` (provided as a tuple (H, states))

## Arguments:
- `H`: a Hamiltonian, specified as a function evaluable at a (vectorial) **k**-point.
- `states`: a selection of states from the `H` (an AbstractVector).
- `kpath`: a vector of k-points, specifying the path of the Wilson loop; it is expected to
  span a range of 2π, i.e. is specified in normalized units.
- `coords`: a vector of coordinates associated with the elements of the Hamiltonian or
  eigenstates; generally, a vector of vectors (dim H × dim R). If `nothing`, all coordinates
  are set to the origin.
"""
function wilsonloop(hamiltonianstatetuple::Tuple{T, AbstractVector} where T<:Function, 
                    kpath::AbstractVector, coords::Union{AbstractVector,Nothing}=nothing)
    H, states = hamiltonianstatetuple
    wilsonloop(kvec -> pickeigvecs!(H(kvec), states), 
               kpath, coords)
end

function pickeigvecs!(A, states)
    if A isa SArray
        F = eigen(Hermitian(A)) 
    else
        F = eigen!(Hermitian(A)) # this is destructive, i.e. overwrites A
    end
    idx = sortperm(F.values, by=real) # make sure eigenvectors are sorted (automatic in julia-1.2.x, but not guaranteed in previous versions)
    u = F.vectors[:, @view idx[states]]
    return normalizevecs!(u)
end

function normalizevecs!(u)
    for un in eachcol(u)
        un ./= norm(un,2)
    end
    return u
end

@docs raw""" 
    overlap(uki, ukj) 

Computes the overlap matrix
```math
M_{nm}^{ij} = \langle u_n(k_i) | u_m(k_j) \rangle
```
given a set of (orthonormal) Bloch states ``{u_n}`` evaluated at (usually adjacent)
**k**-points ``k_i`` and ``k_j``. A *set* of Bloch states is expected.
"""
overlap(uki, ukj) = ( S = size(uki, 2); return overlap!(Matrix{eltype(uki)}(undef, S,S), uki, ukj) )

function overlap!(Mij, uki, ukj)  
    mul!(Mij, uki', ukj)
    unitarize!(Mij)
    return Mij
end

function unitarize!(Mij)
    F = svd!(Mij) # allow overwriting Mij, since we re-overwrite below again (saves some allocations)
    return mul!(Mij, F.U, F.Vt)
end

""" 
    wannierhamiltonian

Identify a topologically equivalent Hamiltonian H_edge from the Wilson loop W via the
identification ``W = \\exp(iH_{\\text{edge}})``
"""
function wannierhamiltonian(statefunc, kpath::AbstractVector, coords::Union{AbstractVector,Nothing}=nothing)
    M = wilsonloop(statefunc, kpath, coords)
    u = statefunc(kpath[1])
    return -1im.*u*log(M)*u' # why is this basis change needed?
end
function wannierhamiltonian(hamiltonianstatetuple::Tuple{T, AbstractVector} where T<:Function, 
                            kpath::AbstractVector, coords::Union{AbstractVector,Nothing}=nothing)
    H, states = hamiltonianstatetuple
    return wannierhamiltonian(kvec -> pickeigvecs!(H(kvec), states), kpath, coords)
end

function wannierstate(statefunc, kpath::AbstractVector, sectorstates, coords::Union{AbstractVector,Nothing}=nothing)
    _, v = wilsonspectrum(wilsonloop(statefunc, kpath, coords)) # eigenstates of the wilson loop; wilson states
    u    = statefunc(kpath[1]) # energy eigenstate at the initial k-point of the wilson loop
    w = u*v[:,sectorstates] # wilson states in the basis of energy states (equivalent to the more tedius sum above)
    return w
end
function wannierstate(hamiltonianstatetuple::Tuple{T, AbstractVector} where T<:Function,
                      kpath::AbstractVector, sectorstates, coords::Union{AbstractVector,Nothing}=nothing)
    H, states = hamiltonianstatetuple
    return wannierstate(kvec -> pickeigvecs!(H(kvec), states), kpath, sectorstates, coords)
end

function spectrum(H::Function, kvec)
    E = eigvals(Hermitian(H(kvec)))
end


domain2π(φ) = φ < 0 ? φ + 2π : φ
angle2π(x) = domain2π(angle(x))
"""
    wilsonspectrum(M)

Return the (possibly non-Abelian) Berry phase spectrum associated with a Wilson loop
matrix `M` via its Berry phases φ (∈[0,2π]) and the associated Wilson loop eigenfunctions
sorted according to φ.
"""
function wilsonspectrum(M)
    F = eigen(M)
    φ = angle.(F.values)
    idx = sortperm(φ)
    return φ[idx], normalizevecs!(F.vectors[:,idx])
end

""" 
    berryphase(M)

Return the (possibly non-Abelian) Berry phases φ (∈[0,2π]) of a Wilson loop matrix `M`, 
in sorted order.
"""
function berryphase(M)
    φ = sort(angle.(eigvals(M)))
end

end # module