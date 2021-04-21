

# --- load packages (just include for now) ----
using Plots
pyplot()

using WilsonLoop
using WilsonLoop.Hamiltonians

# --- scripting ... ---
Nk = 50; 
kxy = range(-π,π,length=Nk)
chooseH = "QSH"
if chooseH == "Chern"
    states = 1:1
    H(kvec) = HChern(kvec[1],kvec[2], 1)
    dimH = 2    
elseif chooseH == "QSH"
    states = 1:2
    H(kvec) = HQSH(kvec[1],kvec[2], 3)
    dimH = 4
elseif chooseH == "Hq" || chooseH == "HOTI"
    states = 1:2
    H(kvec) = Hq(kvec[1],kvec[2])
    dimH = 4
    if chooseH == "HOTI"
        states_loop1 = 1:2
        states_loop2 = 1:1

        # "edge" Hamiltonian approach
        HWx(kx) = wannierhamiltonian((H, states_loop1), [[kx, ky] for ky in kxy])
        HWy(ky) = wannierhamiltonian((H, states_loop1), [[kx, ky] for kx in kxy])
        φxtheny1 = berryphase(wilsonloop((HWx, states_loop2), collect(kxy)))
        φythenx1 = berryphase(wilsonloop((HWy, states_loop2), collect(kxy)))
        @show φxtheny1
        @show φythenx1

        # direct wilson state approach, in u basis
        wx(kx) = wannierstate((H, states_loop1), [[kx, ky] for ky in kxy], states_loop2)
        wy(ky) = wannierstate((H, states_loop1), [[kx, ky] for kx in kxy], states_loop2)
        φxtheny2 = berryphase(wilsonloop(wx, collect(kxy)))
        φythenx2 = berryphase(wilsonloop(wy, collect(kxy)))
        @show φxtheny2
        @show φythenx2
    end
end


φ = Matrix{Float64}(undef, Nk, length(states))
for kk = 1:Nk
    kpath= [[kxy[kk], ky] for ky in kxy] # path along ky for fixed kx (=>HWx(kx))
    φ[kk, :] = berryphase(wilsonloop((H, states), kpath))
end

pW=plot(kxy, φ);
xlims!(extrema(kxy))
ylims!((-π,π))

# spectrum 
Nkiir = 100;
kirrx = [collect(range(0,stop=pi,length=Nkiir))[1:end-1]
         collect(range(float(pi), stop=pi, length=Nkiir))[1:end-1]
         collect(range(float(pi), stop=0, length=Nkiir))]
kirry = [collect(range(0,stop=0,length=Nkiir))[1:end-1]
         collect(range(0,stop=pi,length=Nkiir))[1:end-1]
         collect(range(float(pi), stop=0, length=Nkiir))]
kirrl = vec([0 cumsum(sqrt.((kirrx[2:end]-kirrx[1:end-1]).^2 + (kirry[2:end]-kirry[1:end-1]).^2))...])
E = Matrix{Float64}(undef, length(kirrl), dimH)
for kk=1:length(kirrx);
    E[kk, :] = spectrum(H, [kirrx[kk], kirry[kk]])
end
pE=plot(kirrl, E);
plot(pW,pE,layout=2)
#plot(kirrx,kirry)
