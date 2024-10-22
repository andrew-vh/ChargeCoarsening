using DifferentialEquations, Random, Statistics
using Plots
using LinearAlgebra
using Interpolations
using ProgressLogging
using ForwardDiff
using BenchmarkTools
using QuadGK
using FFTW
using SparseDiffTools
using SparseArrays
using Printf
using Symbolics
using CSV
using DataFrames
using DelimitedFiles
using JLD2

ENV["GKSwstype"]="nul"
 
gr()
default(display_type=:inline)


function initialize_grid(nx, ny, L)
    dx = L / nx
    dy = L / ny
    x_centers = range(dx/2, stop=L-dx/2, length=nx) |> collect
    y_centers = range(dy/2, stop=L-dy/2, length=ny) |> collect
    return dx, dy, x_centers, y_centers
end

function fill_mat(nx, ny, phi)
    phi_mat=Matrix{Any}(undef, nx+2, ny+2)

    phi_mat[2:nx+1,2:ny+1]=phi
    phi_mat[1,2:ny+1]=phi[end,:]
    phi_mat[end,2:ny+1]=phi[1,:]
    phi_mat[2:nx+1,1]=phi[:,end]
    phi_mat[2:nx+1, end]=phi[:,1]
    phi_mat[1,1]=(phi_mat[1,2]+phi_mat[2,1])/2
    phi_mat[end,end]=(phi_mat[end,end-1]+phi_mat[end-1,end])/2
    phi_mat[1,end]=(phi_mat[2,end]+phi_mat[1,end-1])/2
    phi_mat[end,1]=(phi_mat[end,2]+phi_mat[end-1,1])/2
    return phi_mat

end

function initialize_solution(x_centers, y_centers, p)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal= p
    nx=length(x_centers)
    ny=length(y_centers)
    phi0_vals=reshape(phi0*(ones(Int,nx, ny) .+0.01*randn((nx, ny))), (nx*ny,1))
    phicat0_vals=reshape(phicat0*(ones(Int,nx, ny) .+0.01*randn((nx, ny))), (nx*ny,1))
    phian0_vals=reshape(phian0*(ones(Int,nx, ny) .+0.01*randn((nx, ny))), (nx*ny,1))
    phian0_vals=phian0_vals .+(sum(sigma*phi0_vals .+ phicat0_vals .-phian0_vals))/(nx*ny)
    println((sum(sigma*phi0_vals .+ phicat0_vals .-phian0_vals))/(nx*ny))
    u0_vals=-(sigma*phi0_vals .+ phicat0_vals .-phian0_vals)*0
    psi0_vals=-phi0_vals.+phi0
    z0 =vcat(phi0_vals, phicat0_vals,phian0_vals, u0_vals, psi0_vals)
    return z0
end


# function dae(dz, z, p, t)  
#     L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal= p
#     nx=convert(Int,nx)
#     ny=convert(Int, ny)
    
#     # Reshape z into matrices
#     @views phi = reshape(z[1:nx*ny], nx, ny)
#     @views phicat = reshape(z[nx*ny+1:2*nx*ny], nx, ny)
#     @views phian = reshape(z[2*nx*ny+1:3*nx*ny], nx, ny)
#     @views u = reshape(z[3*nx*ny+1:4*nx*ny], nx, ny)
#     @views psi = reshape(z[4*nx*ny+1:5*nx*ny], nx, ny)

#     # Preallocate matrices
#     @views phi_mat = fill_mat(nx, ny, phi)
#     @views phicat_mat = fill_mat(nx, ny, phicat)
#     @views phian_mat = fill_mat(nx, ny, phian)
#     @views u_mat = fill_mat(nx, ny, u)
#     @views psi_mat = fill_mat(nx, ny, psi)

#     # Reshape dz into 3D array
#     @views dz = reshape(dz, nx, ny, 5)

#     @views mu_mat=similar(phi_mat)

#     @views mu=log.(abs.(phi).^(1/N)./abs.(1 .-phi) .+1e-6).-2*chi*phi.+ psi
#     @views mu_mat=fill_mat(nx,ny, mu)


#     @views phi_flux_x=-(phi_mat[2:end,2:end-1].+phi_mat[1:end-1,2:end-1])/2 .*(mu_mat[2:end,2:end-1] .-mu_mat[1:end-1,2:end-1])/dx
#     @views phi_flux_y=-(phi_mat[2:end-1,2:end].+phi_mat[2:end-1,1:end-1])/2 .*(mu_mat[2:end-1,2:end] .-mu_mat[2:end-1,1:end-1])/dy


#     @views phicat_flux_x=-D*(phicat_mat[2:end,2:end-1] .-phicat_mat[1:end-1,2:end-1])/dx-D*(phicat_mat[2:end,2:end-1].+phicat_mat[1:end-1,2:end-1])/2 .*(u_mat[2:end,2:end-1] .-u_mat[1:end-1,2:end-1])/dx
#     @views phicat_flux_y=-D*(phicat_mat[2:end-1,2:end] .-phicat_mat[2:end-1,1:end-1])/dy-D*(phicat_mat[2:end-1,2:end].+phicat_mat[2:end-1,1:end-1])/2 .*(u_mat[2:end-1,2:end] .-u_mat[2:end-1,1:end-1])/dy


#     @views phian_flux_x=-D*(phian_mat[2:end,2:end-1] .-phian_mat[1:end-1,2:end-1])/dx+D*(phian_mat[2:end,2:end-1].+phian_mat[1:end-1,2:end-1])/2 .*(u_mat[2:end,2:end-1] .-u_mat[1:end-1,2:end-1])/dx
#     @views phian_flux_y=-D*(phian_mat[2:end-1,2:end] .-phian_mat[2:end-1,1:end-1])/dy+D*(phian_mat[2:end-1,2:end].+phian_mat[2:end-1,1:end-1])/2 .*(u_mat[2:end-1,2:end] .-u_mat[2:end-1,1:end-1])/dy

#     @views u_flux_x=-lambda^2*(u_mat[2:end,2:end-1] .-u_mat[1:end-1,2:end-1])/dx
#     @views u_flux_y=-lambda^2*(u_mat[2:end-1,2:end] .-u_mat[2:end-1,1:end-1])/dy


#     @views psi_flux_x=-(phi_mat[2:end,2:end-1] .-phi_mat[1:end-1,2:end-1])/dx
#     @views psi_flux_y=-(phi_mat[2:end-1,2:end] .-phi_mat[2:end-1,1:end-1])/dy
    
#     @views dz[:,:,1]=-(phi_flux_x[2:end,:]-phi_flux_x[1:end-1,:])/dx - (phi_flux_y[:,2:end]-phi_flux_y[:,1:end-1])/dy
#     @views dz[:,:,2]=-(phicat_flux_x[2:end,:]-phicat_flux_x[1:end-1,:])/dx - (phicat_flux_y[:,2:end]-phicat_flux_y[:,1:end-1])/dy
#     @views dz[:,:,3]=-(phian_flux_x[2:end,:]-phian_flux_x[1:end-1,:])/dx - (phian_flux_y[:,2:end]-phian_flux_y[:,1:end-1])/dy
#     @views dz[:,:,4]=-(u_flux_x[2:end,:]-u_flux_x[1:end-1,:])/dx - (u_flux_y[:,2:end]-u_flux_y[:,1:end-1])/dy.+sign.(phicat.-phian.+sigma*phi).*max.(abs.(phicat.-phian.+sigma*phi),1e-8 .*ones((nx,ny)))
#     @views dz[:,:,5]=-(psi_flux_x[2:end,:]-psi_flux_x[1:end-1,:])/dx - (psi_flux_y[:,2:end]-psi_flux_y[:,1:end-1])/dy.+psi
#     dz =reshape(dz, (5*nx*ny,1))

# end

function dae(dz, z, p, t)
    println(t)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal= p
    nx=convert(Int,nx)
    ny=convert(Int, ny)
    
    # converts z vector into matrices for 5 variables
    phi=reshape(z[1:nx*ny], (nx,ny))
    phicat=reshape(z[nx*ny+1:2*nx*ny], (nx,ny))
    phian=reshape(z[2*nx*ny+1:3*nx*ny], (nx,ny))
    u=reshape(z[3*nx*ny+1:4*nx*ny], (nx,ny))
    psi=reshape(z[4*nx*ny+1:5*nx*ny], (nx,ny))

    dz=reshape(dz, (nx,ny, 5))
    
    phi_flux_x=similar(phi)
    phi_flux_y=similar(phi)

    phicat_flux_x=similar(phicat)
    phicat_flux_y=similar(phicat)

    phian_flux_x=similar(phian)
    phian_flux_y=similar(phian)

    u_flux_x=similar(u)
    u_flux_y=similar(u)

    psi_flux_x=similar(psi)
    psi_flux_y=similar(psi)

    mu=similar(phi)
    mu=log.(abs.(phi).^(1/N)./abs.(1 .-phi) .+1e-6).-2*chi*phi.+ psi.+sigma*u

    for i in 1:nx, j in 1:ny
        ip1 = (i == nx) ? 1 : i + 1  # Periodic boundary condition
        im1 = (i == 1) ? nx : i - 1  # Periodic boundary condition
        jp1 = (j == ny) ? 1 : j + 1  # Periodic boundary condition
        jm1 = (j == 1) ? ny : j - 1  # Periodic boundary condition

        phi_flux_x[i,j] = -(phi[ip1,j]+phi[i,j])/2 * (mu[ip1, j] - mu[i, j])/dx
        phi_flux_y[i,j] = -(phi[i,jp1]+phi[i,j])/2 * (mu[i, jp1] - mu[i, j])/dy 

        phicat_flux_x[i,j] = -D*(phicat[ip1, j] - phicat[i, j])/dx-D*(phicat[ip1,j]+phicat[i,j])/2 * (u[ip1, j] - u[i, j])/dx
        phicat_flux_y[i,j] = -D*(phicat[i, jp1] - phicat[i, j])/dy-D*(phicat[i,jp1]+phicat[i,j])/2 * (u[i, jp1] - u[i, j])/dy

        phian_flux_x[i,j] = -D*(phian[ip1, j] - phian[i, j])/dx+D*(phian[ip1,j]+phian[i,j])/2 * (u[ip1, j] - u[i, j])/dx
        phian_flux_y[i,j] = -D*(phian[i, jp1] - phian[i, j])/dy+D*(phian[i,jp1]+phian[i,j])/2 * (u[i, jp1] - u[i, j])/dy

        u_flux_x[i,j]=-lambda^2*(u[ip1, j] - u[i, j])/dx
        u_flux_y[i,j]=-lambda^2*(u[i, jp1] - u[i, j])/dy


        psi_flux_x[i,j]=-(phi[ip1, j] - phi[i, j])/dx
        psi_flux_y[i,j]=-(phi[i, jp1] - phi[i, j])/dy
    end
    for i in 1:nx, j in 1:ny
        ip1 = (i == nx) ? 1 : i + 1  # Periodic boundary condition
        im1 = (i == 1) ? nx : i - 1  # Periodic boundary condition
        jp1 = (j == ny) ? 1 : j + 1  # Periodic boundary condition
        jm1 = (j == 1) ? ny : j - 1  # Periodic boundary condition
        dz[i, j, 1] =  -(phi_flux_x[i,j]-phi_flux_x[im1,j])/dx - (phi_flux_y[i,j]-phi_flux_y[i,jm1])/dy
        dz[i, j, 2] =  -(phicat_flux_x[i,j]-phicat_flux_x[im1,j])/dx - (phicat_flux_y[i,j]-phicat_flux_y[i,jm1])/dy
        dz[i, j, 3] =  -(phian_flux_x[i,j]-phian_flux_x[im1,j])/dx - (phian_flux_y[i,j]-phian_flux_y[i,jm1])/dy
        # dz[i, j, 4]= -(u_flux_x[i,j]-u_flux_x[im1,j])/dx-(u_flux_y[i,j]-u_flux_y[i,jm1])/dy+sign((phicat[i,j]-phian[i,j]+sigma*phi[i,j]))*maximum((abs(phicat[i,j]-phian[i,j]+sigma*phi[i,j]), 1e-8))
        dz[i, j, 4]= -(u_flux_x[i,j]-u_flux_x[im1,j])/dx-(u_flux_y[i,j]-u_flux_y[i,jm1])/dy+sign((phicat[i,j]-phian[i,j]+sigma*phi[i,j]))*maximum((abs(phicat[i,j]-phian[i,j]+sigma*phi[i,j])))
        dz[i, j, 5]= -(psi_flux_x[i,j]-psi_flux_x[im1,j])/dx-(psi_flux_y[i,j]-psi_flux_y[i,jm1])/dy+psi[i,j]

    end
    dz =reshape(dz, (5*nx*ny,1))
end

function run_simulation(p)
    @time begin
        L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal=p 
        nx=convert(Int,nx)
        ny=convert(Int, ny)
        dx, dy, x_centers, y_centers = initialize_grid(nx, ny,L)
        z0 = initialize_solution(x_centers, y_centers, p)
        
        tspan = (0.0, tfinal)

        # Solve the DAE using ode15s
        M=I(5*nx*ny)
        M[3*nx*ny+1:end,3*nx*ny+1:end]=zeros((2*nx*ny, 2*nx*ny))
        
        sparse_bool=true
        
        if sparse_bool
            dz0=z0*0.0
            jac_sparsity = Symbolics.jacobian_sparsity((dz, z) -> dae(dz, z, p, 0.0), dz0, z0)
            f = ODEFunction(dae, mass_matrix=M, jac_prototype=float.(jac_sparsity))
        else
            f = ODEFunction(dae, mass_matrix=M)
        end
        prob = ODEProblem(f, z0, tspan,p)
        cb = ContinuousCallback(condition, affect!)
        sol = solve(prob,Rosenbrock23(autodiff=false), progress = true, callback=cb)
    end
    return x_centers, y_centers, sol
end

function condition(z, t, integrator)
    p=integrator.p
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal=p 
    phi=z[1:Int(nx*ny)]
    return minimum(phi)
end

function affect!(integrator)
    p=integrator.p
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal=p 
    integrator.u[1:Int(nx*ny)] .= max.(integrator.u[1:Int(nx*ny)], 0.0)
end

function plot_solution(x_centers, y_centers, z, title_str, zlims, p)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal=p 
    nx=convert(Int,nx)
    ny=convert(Int, ny)
    # surface(x_centers, y_centers, z, xlabel="x", ylabel="y", zlabel="u", title=title_str,  camera=(0, 90), c=:viridis, zlims=(0, 1), clim=(0, 1))
    surface(x_centers, y_centers, z, xlabel="x", ylabel="y", zlabel="u", title=title_str,  camera=(0, 90), c=:viridis, clim=(0, 1), zlims=(0,1))
end

function compute_structure_factor(field::Matrix{Float64})
    # Perform 2D Fourier Transform of the field
    F = fft(field.-mean(field))

    # Compute the magnitude squared (power spectrum)
    S = abs2.(F)

    return S
end

# Function to compute radial average
function radial_average(S::Matrix{Float64})
    nx, ny = size(S)
    center_x, center_y = div(nx, 2), div(ny, 2)
    radial_bins = collect(0:0.1:maximum(hypot.(1:nx .- center_x, 1:ny .- center_y)))
    radial_profile = zeros(length(radial_bins))

    bin_counts = zeros(length(radial_bins))

    for i in 1:nx
        for j in 1:ny
            r = hypot(i - center_x, j - center_y)
            bin = searchsortedfirst(radial_bins, r)
            if bin > 0 && bin <= length(radial_bins)
                radial_profile[bin] += S[i, j]
                bin_counts[bin] += 1
            end
        end
    end

    # Avoid division by zero by using only non-zero bin counts
    radial_profile = radial_profile ./ max.(bin_counts, 1)

    return radial_bins, radial_profile
end

function trapezoidal_integration(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    integral = 0.0
    for i in 1:n-1
        integral += 0.5 * (x[i+1] - x[i]) * (y[i+1] + y[i])
    end
    return integral
end

function find_R(phi)
    S = compute_structure_factor(phi)
    radial_bins, radial_profile = radial_average(S)
    integral1 = trapezoidal_integration(radial_bins, radial_profile)
    integral2=trapezoidal_integration(radial_bins, radial_bins.*radial_profile)
    R=integral2/integral1
    return R
end


# Example usage
nx::Int = 20  # Number of spatial grid points in x-direction
ny::Int = 20 # Number of spatial grid points in y-direction
L=20


N=100

dx = L/ nx
dy = L / ny

D=sqrt(N)
lambda=0.6
sigma=0.1
chi=1.2
phi0=0.1
phicat0=0.001
phian0=phicat0+sigma*phi0
tfinal=100
dt=tfinal/100

p=[L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal]

x_centers, y_centers, sol = run_simulation(p)

mkdir("test3")
cd("test3")

# Time points where you want to interpolate
interpolation_times = 0:dt:tfinal
Rt=similar(interpolation_times)
anim = @animate for j=1:length(interpolation_times)
    interpolated_solution=sol(interpolation_times[j])
    phi_values = reshape(interpolated_solution[1:nx*ny],(nx,ny))
    zlims = (0, 1)
    plot_solution(x_centers, y_centers, phi_values,string(interpolation_times[j]), zlims, p)
    Rt[j]=find_R(phi_values)
end 
gif(anim, "animation.gif", fps=4)

parameters = "N=$N, σ=$sigma, ϕ0=$phi0, ϕ+=$phicat0, χ=$chi"
open("parameters.txt", "w") do file
    write(file, parameters)
end

@save "solution.jld2" sol

# Create the plot
plot(interpolation_times, Rt, xlabel="t", ylabel="R(t)", title="Plot of Rt vs t", legend=false)
savefig("Rplot.pdf")

writedlm("Rdata.csv", Rt, ',')

anim2 = @animate for j=1:length(interpolation_times)
    interpolated_solution=sol(interpolation_times[j])
    u_values = reshape(interpolated_solution[3*nx*ny+1:4*nx*ny],(nx,ny))
    zlims = (-1, 1)
    plot_solution(x_centers, y_centers, u_values,string(interpolation_times[j]),p, (-1, 1))
end 
gif(anim2, "potential.gif", fps=4)

anim3 = @animate for j=1:length(interpolation_times)
    interpolated_solution=sol(interpolation_times[j])
    salt_values = reshape(interpolated_solution[nx*ny+1:2*nx*ny],(nx,ny)).+reshape(interpolated_solution[2*nx*ny+1:3*nx*ny],(nx,ny))
    zlims = (0, 0.1)
    plot_solution(x_centers, y_centers, salt_values,string(interpolation_times[j]),p, (0, 0.1))
end 
gif(anim3, "salt.gif", fps=4)