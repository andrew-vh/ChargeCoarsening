using Plots
using Roots

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

k = 0:0.01:5
N = 50
sList = dispersion.(k, N, 0.1, 0.001, 1, 0.1, 0.6, sqrt(N))
sFunction(k) = dispersion(k, N, 0.1, 0.001, 1, 0.1, 0.6, sqrt(N))
gr()
display(plot(k, sList, ylims=(-0.003, 0.003)))

intervals = find_zeros(sFunction, k)
println(intervals)
println(length(intervals))

σstep = 0.001
Nlist = 1:1:100
σValues = 0:σstep:0.3

wavelength = Matrix{Float64}(undef, length(Nlist), length(σValues))

for N in Nlist
    for σ in 1:length(σValues)
        sFunc(k) = dispersion(k, N, 0.1, 0.001, 1, (σ-1)*σstep, 0.6, sqrt(N))
        intervals = find_zeros(sFunc, k)
        len = length(intervals)
        if len==1
            wavelength[N, σ] = 0.7
        elseif len==2
            wavelength[N, σ] = -0.05
        else
            wavelength[N, σ] = intervals[2]
        end
    end
    println(N)
end

cmap = cgrad(:oslo, rev=true)

display(surface(σValues, Nlist, wavelength, 
zlim=(-0.0501,0.7001), camera=(0, 90), c=cmap, xlabel="σ", ylabel="N", title="Instability Phase Diagram"))
png("phasediagram.png")