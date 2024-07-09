using Plots
using Roots

cmap = cgrad(:oslo, rev=true)

function dispersion(k, N, ϕ0, c0, χ, σ, λ, D)
    a = 2*(-λ^2*N+λ^2*N*ϕ0)
    b = -(k^2*λ^2 + 2*c0*D*N + D*k^2*λ^2*N - k^2*λ^2*ϕ0 - 2*c0*D*N*ϕ0 +
    k^2*λ^2*N*ϕ0 - 2*χ*k^2*λ^2*N*ϕ0 - D*k^2*λ^2*N*ϕ0 +
    k^4*λ^2*N*ϕ0 + 2*χ*k^2*λ^2*N*ϕ0^2 - k^4*λ^2*N*ϕ0^2 +
    D*N*ϕ0*σ - D*N*ϕ0^2*σ + N*ϕ0*σ^2 - N*ϕ0^2*σ^2)
    c = -2*c0*D*k^2 - D*k^4*λ^2 + 2*c0*D*k^2*ϕ0 +
        D*k^4*λ^2*ϕ0 - 2*c0*D*k^2*N*ϕ0 +
        4*χ*c0*D*k^2*N*ϕ0 - 2*c0*D*k^4*N*ϕ0 -
        D*k^4*λ^2*N*ϕ0 + 2*χ*D*k^4*λ^2*N*ϕ0 - D*k^6*λ^2*N*ϕ0 -
        4*χ*c0*D*k^2*N*ϕ0^2 + 2*c0*D*k^4*N*ϕ0^2 -
        2*χ*D*k^4*λ^2*N*ϕ0^2 + D*k^6*λ^2*N*ϕ0^2 - D*k^2*ϕ0*σ +
        D*k^2*ϕ0^2*σ - D*k^2*N*ϕ0^2*σ + 2*χ*D*k^2*N*ϕ0^2*σ -
        D*k^4*N*ϕ0^2*σ - 2*χ*D*k^2*N*ϕ0^3*σ + D*k^4*N*ϕ0^3*σ -
        D*k^2*N*ϕ0*σ^2 + D*k^2*N*ϕ0^2*σ^2
    
    return (-b - safesqrt(b^2-4*a*c))/(2*a)
end

function safesqrt(x)
    if x>0
        return sqrt(x)
    else
        return 0
    end
end

k = 0:0.01:100
σstep = 0.001
Nlist = 1:1:100
σValues = 0:σstep:0.3

kmat = zeros(Float64, length(Nlist), length(σValues))

for N in Nlist
    for σ in 1:length(σValues)
        sFunc(k) = dispersion(k, N, 0.1, 0.001, 0.8, (σ-1)*σstep, 0.6, sqrt(N))
        intervals = find_zeros(sFunc, k)
        len = length(intervals)
        if len==1
            kmat[N, σ] = 5
        elseif len==2 && sFunc(intervals[2]/2) > 0
            kmat[N, σ] = 0.0005
        else
            kmat[N, σ] = intervals[2]
        end
    end
    println(N)
end

kmatcopydefault = copy(kmat)
kmatΦHigh = copy(kmat)
kmatΦLow = copy(kmat)
kmatcHigh = copy(kmat)
kmatcLow = copy(kmat)
kmatχHigh = copy(kmat)
kmatχLow = copy(kmat)

λmat = 5 ./ kmat
λmat = log10.(λmat)

display(surface(σValues, Nlist, λmat, 
zlim=(-0.0001, 4.0001), camera=(0, 90), c=cmap,
xticks=0:0.1:0.5, yticks=(0:20:100,(0:20:100...,:right)), zticks=false))
png("χLow.png")


