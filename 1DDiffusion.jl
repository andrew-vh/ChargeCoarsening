using Pkg
#Pkg.add("DifferentialEquations")
using DifferentialEquations
using LinearAlgebra
using Plots
gr()

# Determines the centers of the grid sites.
function initialize_grid(n)
    dx = 1/n
    centers = dx/2:dx:1-dx/2
    return dx, centers
end

# Gives the concentration at a position given a Gaussian distribution.
function gaussian(x, μ, σ)
    return exp(-(x-μ)^2/(2*σ^2))
end

# Sets the initial concentration of each grid site with a Gaussian distribution.
function fill_gaussian(centers, n, μ, σ)
    ϕ = Array{Float64}(undef, n)
    for i in 1:n
        ϕ[i] = gaussian(centers[i], μ, σ)
    end
    return ϕ
end

# Defines the differential equation that describes diffusion.
function ode!(dϕ, ϕ, p, t)
    n = Int64(p[1])
    dϕ[1] = (ϕ[2] + ϕ[n] - 2*ϕ[1]) / (p[2])^2
    for i in 2:n-1
        dϕ[i] = (ϕ[i+1] + ϕ[i-1] - 2*ϕ[i]) / (p[2])^2
    end
    dϕ[n] = (ϕ[n-1] + ϕ[1] - 2*ϕ[n]) / (p[2])^2
end

function main(n, μ, σ, t, name)
    dx, centers = initialize_grid(n)
    ϕ = fill_gaussian(centers, n, μ, σ)
    tspan = (0, t)
    f = ODEFunction(ode!)
    p = [n, dx]
    prob = ODEProblem(f, ϕ, tspan, p)
    sol = solve(prob, Rosenbrock23())

    anim =  @animate for solution in sol.u
        plot(centers, solution, ylims=(0,1))
    end
    gif(anim, "./$name.gif", fps=4)
end

main(100, 0.5, 0.1, 0.1, "periodic")

