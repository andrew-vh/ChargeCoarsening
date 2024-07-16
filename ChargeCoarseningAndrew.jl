using DifferentialEquations
using Random
using Statistics
using Plots
using LinearAlgebra
using Interpolations
using ProgressLogging
using ForwardDiff
using BenchmarkTools
using QuadGK
using FFTW
using SparseDiffTools
gr()

# Computes dx/dy and cell locations
function initialize_grid(nx, ny, L)
    dx = L / nx
    dy = L / ny
    x_centers = range(dx/2, stop=L-dx/2, length=nx) |> collect
    y_centers = range(dy/2, stop=L-dy/2, length=ny) |> collect
    return dx, dy, x_centers, y_centers
end

# Initializes the 3 phis, u, and psi and combines them into one vector
# u = electrostatic potential
# psi = -Laplacian of phi
function initialize_solution(x_centers, y_centers, p)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal= p
    nx = length(x_centers)
    ny = length(y_centers)
    phi0_vals = phi0 * (ones(Float64, nx*ny) .+ 0.01*randn(nx*ny))
    phicat0_vals = phicat0 * (ones(Float64, nx*ny) .+ 0.01*randn(nx*ny))
    phian0_vals = phian0 * (ones(Float64, nx*ny) .+ 0.01*randn(nx*ny)) # arbitrarily initializes phian0
    phian0_vals = phian0_vals .+ (sum(sigma*phi0_vals .+ phicat0_vals .- phian0_vals))/(nx*ny) # corrects phian0 for electroneutrality
    u0_vals = zeros(Float64, nx*ny) # unimportant guess
    psi0_vals = zeros(Float64, nx*ny) # unimportant guess
    z0 = vcat(phi0_vals, phicat0_vals, phian0_vals, u0_vals, psi0_vals)
    return z0
end

# Differential algebraic equation describing the system's dynamics
function dae(dz, z, p, t)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal= p
    nx=convert(Int, nx)
    ny=convert(Int, ny)

    # converts z vector into matrices for 5 variables
    phi=reshape(z[1:nx*ny], (nx,ny))
    phicat=reshape(z[nx*ny+1:2*nx*ny], (nx,ny))
    phian=reshape(z[2*nx*ny+1:3*nx*ny], (nx,ny))
    u=reshape(z[3*nx*ny+1:4*nx*ny], (nx,ny))
    psi=reshape(z[4*nx*ny+1:5*nx*ny], (nx,ny))
    dz=reshape(dz, (nx,ny, 5))

    Jx=similar(phi)
    Jy=similar(phi)
    Jcatx=similar(phicat)
    Jcaty=similar(phicat)
    Janx=similar(phian)
    Jany=similar(phian)
    Jux=similar(u)
    Juy=similar(u)
    Jpsix=similar(psi)
    Jpsiy=similar(psi)
    
    mu=similar(phi)
    mu=log.(abs.(phi).^(1/N)./abs.(1 .-phi).+1e-6).-2*chi*phi.+ psi .+ sigma*u
    
    # calculates all fluxes for each site
    for i in 1:nx, j in 1:ny
        ip1 = (i==nx) ? 1 : i+1
        jp1 = (j==ny) ? 1 : j+1

        phiavgx = (phi[i,j] + phi[ip1,j])/2
        phiavgy = (phi[i,j] + phi[i,jp1])/2
        Jx[i,j] = -phiavgx * (mu[ip1, j] - mu[i, j])/dx
        Jy[i,j] = -phiavgy * (mu[i, jp1] - mu[i, j])/dy
        # Jx[i,j] = - phiavgx * ((1/(N*phi[i,j]) + 1/(1-phi[i,j]) - 2*chi)*(phi[ip1, j] - phi[i, j]) -
        #             (psi[ip1, j] - psi[i, j]) + sigma*(u[ip1, j] - u[i, j]))/dx
        # Jy[i,j] = - phiavgy * ((1/(N*phi[i,j]) + 1/(1-phi[i,j]) - 2*chi)*(phi[i, jp1] - phi[i, j]) -
        #             (psi[i, jp1] - psi[i, j]) + sigma*(u[i, jp1] - u[i, j]))/dy

        Jcatx[i,j] = - D*(phicat[ip1, j] - phicat[i, j])/dx -D* (phicat[i,j] + phicat[ip1,j])/2 *
                        (u[ip1, j] - u[i, j])/dx
        Jcaty[i,j] = - D*(phicat[i, jp1] - phicat[i, j])/dy -D* (phicat[i,j] + phicat[i,jp1])/2 *
                        (u[i, jp1] - u[i, j])/dy

        Janx[i,j] = - D*(phian[ip1, j] - phian[i, j])/dx - D*(phian[i,j] + phian[ip1,j])/2 *
                        (u[ip1, j] - u[i, j])/dx
        Jany[i,j] = - D*(phian[i, jp1] - phian[i, j])/dy - D*(phian[i,j] + phian[i,jp1])/2 *
                        (u[i, jp1] - u[i, j])/dy

        # Jcatx[i,j] = -D*(phicat[ip1, j] - phicat[i, j])/dx-D*(phicat[ip1,j]+phicat[i,j])/2 * (u[ip1, j] - u[i, j])/dx
        # Jcaty[i,j] = -D*(phicat[i, jp1] - phicat[i, j])/dy-D*(phicat[i,jp1]+phicat[i,j])/2 * (u[i, jp1] - u[i, j])/dy

        # Janx[i,j] = -D*(phian[ip1, j] - phian[i, j])/dx+D*(phian[ip1,j]+phian[i,j])/2 * (u[ip1, j] - u[i, j])/dx
        # Jany[i,j] = -D*(phian[i, jp1] - phian[i, j])/dy+D*(phian[i,jp1]+phian[i,j])/2 * (u[i, jp1] - u[i, j])/dy

        Jux[i,j] = -lambda^2*(u[ip1, j] - u[i, j])/dx
        Juy[i,j] = -lambda^2*(u[i, jp1] - u[i, j])/dy

        Jpsix[i,j] = -(phi[ip1, j] - phi[i, j])/dx
        Jpsiy[i,j] = -(phi[i, jp1] - phi[i, j])/dy
    end

    # determines dz for all variables
    for i in 1:nx, j in 1:ny
        im1 = (i==1) ? nx : i-1
        jm1 = (j==1) ? ny : j-1

        dz[i,j,1] = -(Jx[i,j]-Jx[im1,j])/dx - (Jy[i,j]-Jy[i,jm1])/dy
        dz[i,j,2] = -(Jcatx[i,j]-Jcatx[im1,j])/dx - (Jcaty[i,j]-Jcaty[i,jm1])/dy
        dz[i,j,3] = -(Janx[i,j]-Janx[im1,j])/dx - (Jany[i,j]-Jany[i,jm1])/dy
        dz[i,j,4] = -(Jux[i,j]-Jux[im1,j])/dx - (Juy[i,j]-Juy[i,jm1])/dy + sign((phicat[i,j]-phian[i,j]+sigma*phi[i,j]))*maximum((abs(phicat[i,j]-phian[i,j]+sigma*phi[i,j]), 1e-8))
        dz[i,j,5] = -(Jpsix[i,j]-Jpsix[im1,j])/dx - (Jpsiy[i,j]-Jpsiy[i,jm1])/dy + psi[i,j]
    end

    dz = reshape(z, (5*nx*ny, 1))
    println(t)
    return dz
end

# Runs the simulation given input parameters
function run_simulation(p)
    @time begin
        L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal=p 
        nx=convert(Int,nx)
        ny=convert(Int, ny)
        dx, dy, x_centers, y_centers = initialize_grid(nx, ny,L)
        z0 = initialize_solution(x_centers, y_centers, p)
        tspan = (0.0, tfinal)

        # Solve DAE
        M = I(5*nx*ny)
        M[3*nx*ny+1:end,3*nx*ny+1:end]=zeros((2*nx*ny, 2*nx*ny))
        sparse_jac = (J, dz, z, p, t) -> SparseDiffTools.jacobian!(J, dae, dz, z, p, t)
        f = ODEFunction(dae, mass_matrix=M, jac_prototype=sparse_jac)
        prob = ODEProblem(f, z0, tspan, p)
        sol = solve(prob,Rosenbrock23(),reltol=1e-6,abstol=1e-6, progress = true)
    end
    return x_centers, y_centers, sol
end

function plot_solution(x_centers, y_centers, z, title_str, p)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal=p 
    nx=convert(Int,nx)
    ny=convert(Int, ny)
    surface(x_centers, y_centers, z, xlabel="x", ylabel="y", zlabel="u", title=title_str,  camera=(0, 90), c=:viridis, zlims=(0, 1), clim=(0, 1))
end

# Runs and displays the simulation as a gif
function main(L, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, tfinal, nt, name)
    println("START")
    p = [L, L/nx, L/ny, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, sigma*phi0 + phicat0, tfinal]
    x_centers, y_centers, sol = run_simulation(p)
    t_values = sol.t
    interpolation_times = 0:tfinal/nt:tfinal
    anim = @animate for j=1:nt+1
        println(interpolation_times[j])
        interpolated_solution = sol(interpolation_times[j])
        phi_vals = reshape(interpolated_solution[1:nx*ny],(nx,ny))
        phicat_vals = reshape(interpolated_solution[nx*ny+1:2*nx*ny],(nx,ny))
        phian_vals = reshape(interpolated_solution[2*nx*ny+1:3*nx*ny], (nx,ny))
        println(sum(phi_vals))
        println(sum(phicat_vals))
        println(sum(phian_vals))
        plot_solution(x_centers, y_centers, phi_vals,string(interpolation_times[j]),p)
    end
    gif(anim, "./$name.gif", fps=4)
    println("FINISH")
end

main(10, 10, 10, 1.44, 0.1, 0.6, sqrt(10), 10, 0.1, 0.002, 1000, 100, "test4")