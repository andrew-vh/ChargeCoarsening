using DifferentialEquations, Random, Statistics
using Plots
using LinearAlgebra
using Interpolations
using ProgressLogging
using ForwardDiff
gr()

function initialize_grid(nx, ny, L)
    dx = L/ nx
    dy = L / ny
    x_centers = dx/2:dx:L-dx/2
    y_centers = dy/2:dy:L-dy/2
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

function initialize_solution(x_centers, y_centers,L, p)
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

function dae(dz, z, p,t)
    # println(t)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal= p
    nx=convert(Int,nx)
    ny=convert(Int, ny)
    
    phi=reshape(z[1:nx*ny], (nx,ny))
    phicat=reshape(z[nx*ny+1:2*nx*ny], (nx,ny))
    phian=reshape(z[2*nx*ny+1:3*nx*ny], (nx,ny))
    u=reshape(z[3*nx*ny+1:4*nx*ny], (nx,ny))
    psi=reshape(z[4*nx*ny+1:5*nx*ny], (nx,ny))

    phi_mat=fill_mat(nx,ny, phi)
    phicat_mat=fill_mat(nx,ny, phicat)
    phian_mat=fill_mat(nx,ny, phian)
    u_mat=fill_mat(nx,ny, u)
    psi_mat=fill_mat(nx,ny, psi)

    dz=reshape(dz, (nx,ny, 5))

    mu=similar(phi)

    mu_mat=similar(phi_mat)

    mu=log.(abs.(phi).^(1/N)./abs.(1 .-phi) .+1e-6).-2*chi*phi.+ psi

    mu_mat=log.(abs.(phi_mat).^(1/N)./abs.(1 .-phi_mat) .+1e-6).-2*chi*phi_mat.+ psi_mat

    phi_flux_x=zeros((nx+1, ny))
    phi_flux_y=zeros((nx, ny+1))

    phi_flux_x=-(phi_mat[2:end,2:end-1].+phi_mat[1:end-1,2:end-1])/2 .*(mu_mat[2:end,2:end-1] .-mu_mat[1:end-1,2:end-1])/dx
    phi_flux_y=-(phi_mat[2:end-1,2:end].+phi_mat[2:end-1,1:end-1])/2 .*(mu_mat[2:end-1,2:end] .-mu_mat[2:end-1,1:end-1])/dy
    dz[:,:,1]=-(phi_flux_x[2:end,:]-phi_flux_x[1:end-1,:])/dx - (phi_flux_y[:,2:end]-phi_flux_y[:,1:end-1])/dy

    phicat_flux_x=zeros((nx+1, ny))
    phicat_flux_y=zeros((nx, ny+1))

    phicat_flux_x=-D*(phicat_mat[2:end,2:end-1] .-phicat_mat[1:end-1,2:end-1])/dx-D*(phicat_mat[2:end,2:end-1].+phicat_mat[1:end-1,2:end-1])/2 .*(u_mat[2:end,2:end-1] .-u_mat[1:end-1,2:end-1])/dx
    phicat_flux_y=-D*(phicat_mat[2:end-1,2:end] .-phicat_mat[2:end-1,1:end-1])/dy-D*(phicat_mat[2:end-1,2:end].+phicat_mat[2:end-1,1:end-1])/2 .*(u_mat[2:end-1,2:end] .-u_mat[2:end-1,1:end-1])/dy
    dz[:,:,2]=-(phicat_flux_x[2:end,:]-phicat_flux_x[1:end-1,:])/dx - (phicat_flux_y[:,2:end]-phicat_flux_y[:,1:end-1])/dy

    phian_flux_x=zeros((nx+1, ny))
    phian_flux_y=zeros((nx, ny+1))

    phian_flux_x=-D*(phian_mat[2:end,2:end-1] .-phian_mat[1:end-1,2:end-1])/dx+D*(phian_mat[2:end,2:end-1].+phian_mat[1:end-1,2:end-1])/2 .*(u_mat[2:end,2:end-1] .-u_mat[1:end-1,2:end-1])/dx
    phian_flux_y=-D*(phian_mat[2:end-1,2:end] .-phian_mat[2:end-1,1:end-1])/dy+D*(phian_mat[2:end-1,2:end].+phian_mat[2:end-1,1:end-1])/2 .*(u_mat[2:end-1,2:end] .-u_mat[2:end-1,1:end-1])/dy
    dz[:,:,3]=-(phian_flux_x[2:end,:]-phian_flux_x[1:end-1,:])/dx - (phian_flux_y[:,2:end]-phian_flux_y[:,1:end-1])/dy

    u_flux_x=zeros((nx+1, ny))
    u_flux_y=zeros((nx, ny+1))

    u_flux_x=-lambda^2*(u_mat[2:end,2:end-1] .-u_mat[1:end-1,2:end-1])/dx
    u_flux_y=-lambda^2*(u_mat[2:end-1,2:end] .-u_mat[2:end-1,1:end-1])/dy
    dz[:,:,4]=-(u_flux_x[2:end,:]-u_flux_x[1:end-1,:])/dx - (u_flux_y[:,2:end]-u_flux_y[:,1:end-1])/dy.+sign.(phicat.-phian.+sigma*phi).*max.(abs.(phicat.-phian.+sigma*phi),1e-8 .*ones((nx,ny)))

    psi_flux_x=zeros((nx+1, ny))
    psi_flux_y=zeros((nx, ny+1))

    psi_flux_x=-(phi_mat[2:end,2:end-1] .-phi_mat[1:end-1,2:end-1])/dx
    psi_flux_y=-(phi_mat[2:end-1,2:end] .-phi_mat[2:end-1,1:end-1])/dy
    dz[:,:,5]=-(psi_flux_x[2:end,:]-psi_flux_x[1:end-1,:])/dx - (psi_flux_y[:,2:end]-psi_flux_y[:,1:end-1])/dy.+psi

    dz =reshape(dz, (5*nx*ny,1))
    println(t)

end

function run_simulation(p)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal=p 
    nx=convert(Int,nx)
    ny=convert(Int, ny)
    dx, dy, x_centers, y_centers = initialize_grid(nx, ny,L)
    z0 = initialize_solution(x_centers, y_centers,L, p)
    
    tspan = (0.0, tfinal)

    # Solve the DAE using ode15s
    M=I(5*nx*ny)
    M[3*nx*ny+1:end,3*nx*ny+1:end]=zeros((2*nx*ny, 2*nx*ny))
    f = ODEFunction(dae,mass_matrix=M)
    prob = ODEProblem(f, z0, tspan,p)
    sol = solve(prob,Rosenbrock23(),reltol=1e-6,abstol=1e-6, progress = true)
    return x_centers, y_centers, sol
end

function plot_solution(x_centers, y_centers, z, title_str,p)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal=p 
    nx=convert(Int,nx)
    ny=convert(Int, ny)
    # surface(x_centers, y_centers, z, xlabel="x", ylabel="y", zlabel="u", title=title_str,  camera=(0, 90), c=:viridis, zlims=(0, 1), clim=(0, 1))
    surface(x_centers, y_centers, z, xlabel="x", ylabel="y", zlabel="u", title=title_str,  camera=(0, 90), c=:viridis)
end

# Example usage
nx::Int = 50  # Number of spatial grid points in x-direction
ny::Int = 50 # Number of spatial grid points in y-direction
L=50
N=10

dx = L/ nx
dy = L / ny

D=sqrt(N)
lambda=0.6
sigma=0.0
chi=(1+1/sqrt(N))^2/1.2
phi0=1/(1+sqrt(N))
phicat0=0.002
phian0=0.002+sigma*phi0
tfinal=100.0
dt=10.0

p=[L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal]

x_centers, y_centers, sol = run_simulation(p)
t_values = sol.t  # Time values of simulation

# Time points where you want to interpolate
interpolation_times = 0:dt:tfinal
anim = @animate for j=1:length(interpolation_times)
    println(interpolation_times[j])
    interpolated_solution=sol(interpolation_times[j])
    z_vals=reshape(interpolated_solution[1:nx*ny],(nx,ny))
    # z_vals=reshape(sol[3*nx*ny+1:4*nx*ny,1,j],(nx,ny))
    plot_solution(x_centers, y_centers, z_vals,string(interpolation_times[j]),p)
end 
gif(anim, "./animation.gif", fps=4)
