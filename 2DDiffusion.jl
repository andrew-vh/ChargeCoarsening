using DifferentialEquations
using LinearAlgebra
using Plots
gr()

# Determines the centers of the grid sites.
function initialize_grid(nx, ny, h)
    dx = 1/nx
    dy = h/ny
    x_centers = dx/2:dx:1-dx/2
    y_centers = dy/2:dy:h-dy/2
    return dx, dy, x_centers, y_centers
end

dx, dy, x_centers, y_centers = initialize_grid(50, 50, 1)

# Gives the concentration at a position given a Gaussian distribution.
function gaussian(x, y, μx, μy, σ)
    return exp(-((x-μx)^2+(y-μy)^2)/(2*σ^2))
end

# Sets the initial concentration of each grid site with a Gaussian distribution.
function fill_gaussian(x_centers, y_centers, nx, ny, μx, μy, σ)
    ϕ = Matrix{Float64}(undef, nx, ny)
    for i in 1:nx, j in 1:ny
        ϕ[i, j] = gaussian(x_centers[i], y_centers[j], μx, μy, σ)
    end
    return ϕ
end

# Defines the differential equation that describes diffusion.
function ode!(dϕ, ϕ, p, t)
    nx = Int64(p[1])
    ny = Int64(p[2])
    for i in 1:nx, j in 1:ny
        sum_x = 0
        (i==1) ? sum_x += ϕ[nx, j] : sum_x += ϕ[i-1, j]
        (i==nx) ? sum_x += ϕ[1, j] : sum_x += ϕ[i+1, j]
        sum_x -= (2*ϕ[i, j])

        sum_y = 0
        (j==1) ? sum_y += ϕ[i, ny] : sum_y += ϕ[i, j-1]
        (j==ny) ? sum_y += ϕ[i, 1] : sum_y += ϕ[i, j+1]
        sum_y -= (2*ϕ[i, j])

        dϕ[i, j] = sum_x / (p[3])^2 + sum_y / (p[4])^2
    end
end

function main(nx, ny, h, μx, μy, σ, t)
    dx, dy, x_centers, y_centers = initialize_grid(nx, ny, h)
    ϕ = fill_gaussian(x_centers, y_centers, nx, ny, μx, μy, σ)
    tspan = (0, t)
    f = ODEFunction(ode!)
    p = [nx, ny, dx, dy]
    prob = ODEProblem(f, ϕ, tspan, p)
    sol = solve(prob, Rosenbrock23())
    return sol
end

sol = main(50, 50, 1, 0.3, 0.3, 0.3, 1)

for solution in sol.u
    display(surface(x_centers, y_centers, solution, camera=(0, 90), c=:viridis, cbar=(0,1), clims=(0,1)))
end

function animate(sol, name)
    anim =  @animate for solution in sol.u
        surface(x_centers, y_centers, solution, zlim=(0,1), camera=(0, 90), c=:viridis, cbar=(0,1), clims=(0,1))
    end
    gif(anim, "./$name.gif", fps=5)
end

animate(sol, "2d")