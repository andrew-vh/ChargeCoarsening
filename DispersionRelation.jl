using Plots
using Roots

cmap = cgrad(:oslo, rev=true)

# Returns s as a function of k and parameters
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

# Returns the critical sigma between the blue and white regions
function sigmawhite(Nlist, ϕ0, c0, χ, λ)
    return (λ^2 .- λ^2*ϕ0 .- 2*c0*ϕ0 .* Nlist .+ λ^2*ϕ0 .* Nlist .- 2*χ*λ^2*ϕ0 .* Nlist .+
            2*c0*ϕ0^2 .* Nlist .+ 2*χ*λ^2*ϕ0^2 .* Nlist) ./ (((2*λ-ϕ0)*(-1+ϕ0)*ϕ0) .* Nlist)
end

# Returns the critical sigma between the blue and black regions
function sigmablack(Nlist, ϕ0, c0, χ, λ)
    return 1/2 .* (1 .- 1 ./ Nlist .- 1/(1-ϕ0) .+ 2*χ*ϕ0 .+ 
            safesqrt.(4*c0 .* (4*χ .- 2/(1-ϕ0) .- 2/ϕ0 ./ Nlist).+(-1 .+ 1 ./ Nlist .+ 1/(1-ϕ0) .- 2*χ*ϕ0).^2))
end

# Returns the critical sigma between the blue and white regions for large N
function limsigmawhite(ϕ0, c0, χ, λ)
    return (-2*c0 + λ^2 - 2*χ*λ^2 + 2*ϕ0*c0 + 2*χ*λ^2*ϕ0) / ((2*λ - ϕ0) * (ϕ0 - 1))
end

# Returns the critical sigma between the blue and black regions for large N
function limsigmablack(ϕ0, c0, χ, λ)
    return 1/2 * (1 - 1/(1-ϕ0) + 2*χ*ϕ0 + safesqrt(4*c0*(4*χ - 2/(1-ϕ0)) + (-1+1/(1-ϕ0)-2*χ*ϕ0)^2))
end

# Generates a phase diagram of N and σ given parameters
function phasediagram(ϕ0, c0, χ, λ)
    k = 0:0.01:100
    σstep = 0.001
    Nlist = 1:1:300
    σValues = 0:σstep:0.4

    kmat = zeros(Float64, length(Nlist), length(σValues))

    for N in Nlist
        for σ in 1:length(σValues)
            sFunc(k) = dispersion(k, N, ϕ0, c0, χ, (σ-1)*σstep, λ, sqrt(N))
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

    λmat = 5 ./ kmat
    λmat = log10.(λmat)

    surface(Nlist, σValues, λmat', 
    zlim=(-0.0001, 4.0001), xlim=(0, 300), ylim=(0, 0.4), camera=(0, 90), c=cmap,
    xticks=0:100:300, yticks=0:0.1:0.5, zticks=false)
    plot!(Nlist, sigmawhite(Nlist, ϕ0, c0, χ, λ), 4*ones(length(Nlist)), color=:black, linewidth=3)
    plot!(Nlist, sigmablack(Nlist, ϕ0, c0, χ, λ), 4*ones(length(Nlist)), color=:white, linewidth=3)
end

phasediagram(0.1, 0.01, 1, 0.6)
png("Nσdiagram.png")

# Generates a phase diagram of ϕ0 and σ in the large N limit
function phasediagramϕ(c0, χ, λ)
    k = 0:0.01:100
    ϕvalues = 0:0.002:0.6
    σValues = 0:0.001:0.4

    kmat = zeros(Float64, length(ϕvalues), length(σValues))

    for ϕ0 in 1:length(ϕvalues)
        for σ in 1:length(σValues)
            sFunc(k) = dispersion(k, 10000, (ϕ0-1)*0.002, c0, χ, (σ-1)*0.001, λ, 100)
            intervals = find_zeros(sFunc, k)
            len = length(intervals)
            if len==1
                kmat[ϕ0, σ] = 5
            elseif len==2 && sFunc(intervals[2]/2) > 0
                kmat[ϕ0, σ] = 0.0005
            else
                kmat[ϕ0, σ] = intervals[2]
            end
        end
        println(ϕ0)
    end

    λmat = 5 ./ kmat
    λmat = log10.(λmat)

    surface(ϕvalues, σValues, λmat', 
    zlim=(-0.0001, 4.0001), xlim=(0, 0.6), ylim=(0, 0.4), camera=(0, 90), c=cmap, clim=(0, 4),
    xticks=0:0.2:0.6, yticks=0:0.1:0.4, zticks=false)
    plot!(ϕvalues, limsigmawhite.(ϕvalues, c0, χ, λ), 4*ones(length(ϕvalues)), color=:black, linewidth=3)
    plot!(ϕvalues, limsigmablack.(ϕvalues, c0, χ, λ), 4*ones(length(ϕvalues)), color=:white, linewidth=3)
end

phasediagramϕ(0.001, 1, 0.6)
png("ϕpσdiagram.png")

# Generates a phase diagram of c0 and σ in the large N limit
function phasediagramc(ϕ0, χ, λ)
    k = 0:0.01:100
    cvalues = 0:0.0005:0.05
    σValues = 0:0.001:0.4

    kmat = zeros(Float64, length(cvalues), length(σValues))

    for c0 in 1:length(cvalues)
        for σ in 1:length(σValues)
            sFunc(k) = dispersion(k, 10000, ϕ0, (c0-1)*0.0005, χ, (σ-1)*0.001, λ, 100)
            intervals = find_zeros(sFunc, k)
            len = length(intervals)
            if len==1
                kmat[c0, σ] = 5
            elseif len==2 && sFunc(intervals[2]*0.99) >= 0
                kmat[c0, σ] = 0.0005
            else
                kmat[c0, σ] = intervals[2]
            end
        end
        println(c0)
    end

    print(maximum(kmat))
    print(minimum(kmat))
    λmat = 5 ./ kmat
    λmat = log10.(λmat)

    surface(cvalues, σValues, λmat', 
    zlim=(-0.0001, 4.0001), xlim=(0, 0.05), ylim=(0, 0.4), camera=(0, 90), c=cmap, clim=(0, 4),
    xticks=0:0.01:0.05, yticks=0:0.1:0.4, zticks=false)
    plot!(cvalues, limsigmawhite.(ϕ0, cvalues, χ, λ), 4*ones(length(cvalues)), color=:black, linewidth=3)
    plot!(cvalues, limsigmablack.(ϕ0, cvalues, χ, λ), 4*ones(length(cvalues)), color=:white, linewidth=3)
end

phasediagramc(0.1, 1, 0.6)
png("ϕ+σdiagram.png")

# Generates a phase diagram of χ and σ in the large N limit
function phasediagramχ(ϕ0, c0, λ)
    k = 0:0.01:100
    χvalues = 0.5:0.004:1.5
    σValues = 0:0.001:0.5

    kmat = zeros(Float64, length(χvalues), length(σValues))

    for χ0 in 1:length(χvalues)
        for σ in 1:length(σValues)
            sFunc(k) = dispersion(k, 10000, ϕ0, c0, 0.5+(χ0-1)*0.004, (σ-1)*0.001, λ, 100)
            intervals = find_zeros(sFunc, k)
            len = length(intervals)
            if len==1
                kmat[χ0, σ] = 5
            elseif len==2 && sFunc(intervals[2]*0.99) >= 0
                kmat[χ0, σ] = 0.0005
            else
                kmat[χ0, σ] = intervals[2]
            end
        end
        println(χ0)
    end

    print(maximum(kmat))
    print(minimum(kmat))
    λmat = 5 ./ kmat
    λmat = log10.(λmat)

    surface(χvalues, σValues, λmat', 
    zlim=(-0.0001, 4.0001), xlim=(0.5, 1.5), ylim=(0, 0.4), camera=(0, 90), c=cmap, clim=(0, 4),
    xticks=0.5:0.5:1.5, yticks=0:0.1:0.4, zticks=false)
    plot!(χvalues, limsigmawhite.(ϕ0, c0, χvalues, λ), 4*ones(length(χvalues)), color=:black, linewidth=3)
    plot!(χvalues, limsigmablack.(ϕ0, c0, χvalues, λ), 4*ones(length(χvalues)), color=:white, linewidth=3)
end

phasediagramχ(0.1, 0.001, 0.6)
png("χσdiagram.png")

# kList = zeros(Float64, length(Nlist))
# for N in Nlist
#     sFunc(k) = dispersion(k, N, 0.1, 0.001, 0.8, 0.1, 0.6, sqrt(N))
#     intervals = find_zeros(sFunc, k)
#     len = length(intervals)
#         if len==1
#             kList[N] = 5
#         elseif len==2 && sFunc(intervals[2]/2) > 0
#             kList[N] = 0.0005
#         else
#             kList[N] = intervals[2]
#         end
#     println(N)
# end

# println(kList)

# λlist = 5 ./ kList
# λlist = log10.(λlist)

# display(plot(Nlist, λlist, seriestype=:scatter))

k = 0:0.01:0.5
s = dispersion.(k, 10000, 0.1, 0.1, 1, 0.469, 0.6, 100)
plot(k, s)