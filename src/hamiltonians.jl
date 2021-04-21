using StaticArrays

const σx = @SArray [0 1; 1 0] # use SArrays to make functions of σi inferrable inside functions (e.g. kron)
const σy = @SArray [0 1im; -1im 0]
const σz = @SArray [1 0; 0 -1]
const σ0 = @SArray [1 0; 0 1]

""" 
    HChern(kx, ky, m)

2D Hamiltonian with nonzero Chern number, from
Benalcazar's PRB, Eq. (4.28). `m` is analogous, but
not quite a mass term; seems we have a Dirac model 
when `|m| = 2` and when `m = 0`
"""
function HChern(kx, ky, m)
    sin(kx)*σx  .+  sin(ky)*σy  .+  (m + cos(kx) + cos(ky))*σz
end
HChern(kvec, m) = HChern(kvec[1], kvec[2], m)

""" 
    HQSH(kx, ky, m)

"""
function HQSH(kx, ky, m)
    (sin(kx)*(kron(σz,σx) .+ kron(σx,σx)) .+
    sin(ky)*(kron(σy,σx) .+ kron(σ0,σy)) .+
    (2 - m - cos(kx) - cos(ky))*kron(σ0,σz))
end
HQSH(kvec, m) = HQSH(kvec[1], kvec[2], m)

const Γ0 =  kron(σz, σ0) # not used at the moment...
const Γ1 = -kron(σy, σx)
const Γ2 = -kron(σy, σy)
const Γ3 = -kron(σy, σz)
const Γ4 =  kron(σx, σ0)
""" 
    Hq(kx, ky)

Eq. (6.29) from Benalcazar's paper; has a quadrupole moment
for default parameter choices
"""
function Hq(kx, ky; γx=0.5, γy=0.5, λx=1.0, λy=1.0)
    H = (γx + λx*cos(kx))*Γ4 .+ (λx*sin(kx))*Γ3 .+ (γy + λy*cos(ky))*Γ2 .+ (λy*sin(ky))*Γ1
    return H
end