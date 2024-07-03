using Pkg
Pkg.activate("/home/jsouza/.julia/project1")
Pkg.instantiate()

using DifferentialEquations, Random, Statistics
using Plots
using LinearAlgebra
using Interpolations
using ProgressLogging
using Base.Threads

ENV["GKSwstype"]="nul"
 
gr()
default(display_type=:inline)


function initialize_grid(nx, ny, L)
    dx = L/ nx
    dy = L / ny
    x_centers = dx/2:dx:L-dx/2
    y_centers = dy/2:dy:L-dy/2
    return dx, dy, x_centers, y_centers
end

function initialize_solution(x_centers, y_centers,L, p)
    L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal= p
    nx=length(x_centers)
    ny=length(y_centers)
    phi0_vals=reshape(phi0*(ones(Int,nx, ny) .+0.001*randn((nx, ny))), (nx*ny,1))
    # phi0_vals= reshape(2*ones(Int,length(x_centers), length(y_centers))+sin.(2π * x_centers/L) * transpose(cos.(2π * y_centers/L))+0*randn(nx, ny), (nx*ny,1))
    # phicat0_vals=phi0_vals
    # phian0_vals=phi0_vals
    # u0_vals=phi0_vals
    # psi0_vals=phi0_vals
    phicat0_vals=reshape(phicat0*(ones(Int,nx, ny) .+0.01*randn((nx, ny))), (nx*ny,1))
    phian0_vals=reshape(phian0*(ones(Int,nx, ny) .+0.01*randn((nx, ny))), (nx*ny,1))
    phian0_vals=phian0_vals .+(sum(sigma*phi0_vals .+ phicat0_vals .-phian0_vals))/(nx*ny)
    # phian0_vals=sigma*phi0_vals .+ phicat0_vals
    println((sum(sigma*phi0_vals .+ phicat0_vals .-phian0_vals))/(nx*ny))
    u0_vals=-(sigma*phi0_vals .+ phicat0_vals .-phian0_vals)*0
    psi0_vals=-phi0_vals.-phi0
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

 
   

    mu=-2*chi*phi .+psi .+sigma*u

    # mu=log.(abs.(phi).^(1/N) .+1e-6) .-log.(abs.(1 .-phi).+1e-6) .-2*chi*phi.+ psi

    mu=log.(abs.(phi).^(1/N)./abs.(1 .-phi) .+1e-6).-2*chi*phi.+ psi



    for i in 1:nx, j in 1:ny
        ip1 = (i == nx) ? 1 : i + 1  # Periodic boundary condition
        im1 = (i == 1) ? nx : i - 1  # Periodic boundary condition
        jp1 = (j == ny) ? 1 : j + 1  # Periodic boundary condition
        jm1 = (j == 1) ? ny : j - 1  # Periodic boundary condition

        # phi_flux_x[i,j] = -(1-(phi[ip1,j]+phi[i,j])/2)/N*(phi[ip1, j] - phi[i, j])/dx-((phi[ip1,j]+phi[i,j])/2)*(phi[ip1, j] - phi[i, j])/dx-(1-(phi[ip1,j]+phi[i,j])/2)*(phi[ip1,j]+phi[i,j])/2 * (mu[ip1, j] - mu[i, j])/dx
        # phi_flux_y[i,j] = -(1-(phi[i,jp1]+phi[i,j])/2)/N*(phi[i, jp1] - phi[i, j])/dy-((phi[i,jp1]+phi[i,j])/2)*(phi[i, jp1] - phi[i, j])/dy-(1-(phi[i,jp1]+phi[i,j])/2)*(phi[i,jp1]+phi[i,j])/2 * (mu[i, jp1] - mu[i, j])/dy 
        
        # phi_flux_x[i,j] = -(1-(phi[ip1,j]+phi[i,j])/2)*(phi[ip1,j]+phi[i,j])/2 * (mu[ip1, j] - mu[i, j])/dx
        # phi_flux_y[i,j] = -(1-(phi[i,jp1]+phi[i,j])/2)*(phi[i,jp1]+phi[i,j])/2 * (mu[i, jp1] - mu[i, j])/dy 

        phi_flux_x[i,j] = -(phi[ip1,j]+phi[i,j])/2 * (mu[ip1, j] - mu[i, j])/dx
        phi_flux_y[i,j] = -(phi[i,jp1]+phi[i,j])/2 * (mu[i, jp1] - mu[i, j])/dy 

        # phi_flux_x[i,j] = -(1-phi[i,j])/N*(phi[ip1, j] - phi[i, j])/dx-phi[i,j]*(phi[ip1, j] - phi[i, j])/dx-(1-phi[i,j])*phi[i,j]* (mu[ip1, j] - mu[i, j])/dx
        # phi_flux_y[i,j] = -(1-phi[i,j])/N*(phi[i, jp1] - phi[i, j])/dy-phi[i,j]*(phi[i, jp1] - phi[i, j])/dy-(1-phi[i,j])*phi[i,j]* (mu[i, jp1] - mu[i, j])/dy 

        # phi_flux_x[i,j] = -(phi[ip1, j] - phi[i, j])/dx
        # phi_flux_y[i,j] = -(phi[i, jp1] - phi[i, j])/dy

        # phicat_flux_x[i,j] = -(phicat[ip1, j] - phicat[i, j])/dx
        # phicat_flux_y[i,j] = -(phicat[i, jp1] - phicat[i, j])/dy

        # phian_flux_x[i,j] = -(phian[ip1, j] - phian[i, j])/dx
        # phian_flux_y[i,j] = -(phian[i, jp1] - phian[i, j])/dy

        phicat_flux_x[i,j] = -D*(phicat[ip1, j] - phicat[i, j])/dx-D*(phicat[ip1,j]+phicat[i,j])/2 * (u[ip1, j] - u[i, j])/dx
        phicat_flux_y[i,j] = -D*(phicat[i, jp1] - phicat[i, j])/dy-D*(phicat[i,jp1]+phicat[i,j])/2 * (u[i, jp1] - u[i, j])/dy

        phian_flux_x[i,j] = -D*(phian[ip1, j] - phian[i, j])/dx+D*(phian[ip1,j]+phian[i,j])/2 * (u[ip1, j] - u[i, j])/dx
        phian_flux_y[i,j] = -D*(phian[i, jp1] - phian[i, j])/dy+D*(phian[i,jp1]+phian[i,j])/2 * (u[i, jp1] - u[i, j])/dy

        # phicat_flux_x[i,j] = -D*(phicat[ip1, j] - phicat[i, j])/dx-D*phicat[i,j] * (u[ip1, j] - u[i, j])/dx
        # phicat_flux_y[i,j] = -D*(phicat[i, jp1] - phicat[i, j])/dy-D*phicat[i,j] * (u[i, jp1] - u[i, j])/dy

        # phian_flux_x[i,j] = -D*(phian[ip1, j] - phian[i, j])/dx+D*phian[i,j] * (u[ip1, j] - u[i, j])/dx
        # phian_flux_y[i,j] = -D*(phian[i, jp1] - phian[i, j])/dy+D*phian[i,j] * (u[i, jp1] - u[i, j])/dy
        
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
        dz[i, j, 4]= -(u_flux_x[i,j]-u_flux_x[im1,j])/dx-(u_flux_y[i,j]-u_flux_y[i,jm1])/dy+sign((phicat[i,j]-phian[i,j]+sigma*phi[i,j]))*maximum((abs(phicat[i,j]-phian[i,j]+sigma*phi[i,j]), 1e-8))
        dz[i, j, 5]= -(psi_flux_x[i,j]-psi_flux_x[im1,j])/dx-(psi_flux_y[i,j]-psi_flux_y[i,jm1])/dy+psi[i,j]
        # dz[i, j, 4]= phicat[i,j]-phian[i,j]+sigma*phi[i,j]

        # dz[i,j,4]=u[i,j] 
        # dz[i, j, 4]= -(u_flux_x[i,j]-u_flux_x[im1,j])/dx-(u_flux_y[i,j]-u_flux_y[i,jm1])/dy

        # dz[i,j,5]=psi[i,j]
        # dz[i, j,4] =  -(u_flux_x[i,j]-u_flux_x[im1,j])/dx - (u_flux_y[i,j]-u_flux_y[i,jm1])/dy
        # dz[i, j, 5] =  -(psi_flux_x[i,j]-psi_flux_x[im1,j])/dx - (psi_flux_y[i,j]-psi_flux_y[i,jm1])/dy
        if i==1 & j==1
            dz[i, j, 4]= u[i,j]
        end
    end
    dz =reshape(dz, (5*nx*ny,1))
    println(t)
    # if mod(round(t),5)==0
    #     if mod(round(t*1000),5)==0
    #         println(t)
    #     end
    # end
    
    
   
    

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
    sol = solve(prob,Rosenbrock23(),reltol=1e-6,abstol=1e-6)
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
ny::Int = 50  # Number of spatial grid points in y-direction
L=50
N=10

dx = L/ nx
dy = L / ny

D=sqrt(N)
lambda=0.6
sigma=0.5
chi=(1+1/sqrt(N))^2/1.2
phi0=1/(1+sqrt(N))
phicat0=0.002
phian0=0.002+sigma*phi0
tfinal=3e2
dt=tfinal/100




p=[L, dx, dy, nx, ny, chi, sigma, lambda, D, N, phi0, phicat0, phian0, tfinal]



x_centers, y_centers, sol = run_simulation(p)
t_values = sol.t  # Time values for animation
println(t_values)


# Time points where you want to interpolate
interpolation_times = 0:dt:tfinal

interpolated_solution=sol(interpolation_times[end])
phi=reshape(interpolated_solution[1:nx*ny],(nx,ny))
phicat=reshape(interpolated_solution[nx*ny+1:2*nx*ny],(nx,ny))
phian=reshape(interpolated_solution[2*nx*ny+1:3*nx*ny],(nx,ny))
u=reshape(interpolated_solution[3*nx*ny+1:4*nx*ny],(nx,ny))
psi=reshape(interpolated_solution[4*nx*ny+1:5*nx*ny],(nx,ny))
z_vals=phi
plot_solution(x_centers, y_centers, z_vals,string(interpolation_times[end]),p)

anim = @animate for j=1:length(interpolation_times)
    println(interpolation_times[j])

    interpolated_solution=sol(interpolation_times[j])
    phi=reshape(interpolated_solution[1:nx*ny],(nx,ny))
    phicat=reshape(interpolated_solution[nx*ny+1:2*nx*ny],(nx,ny))
    phian=reshape(interpolated_solution[2*nx*ny+1:3*nx*ny],(nx,ny))
    u=reshape(interpolated_solution[3*nx*ny+1:4*nx*ny],(nx,ny))
    psi=reshape(interpolated_solution[4*nx*ny+1:5*nx*ny],(nx,ny))
    z_vals=phi
    # z_vals=reshape(interpolated_solution[4*nx*ny+1:5*nx*ny],(nx,ny))
    plot_solution(x_centers, y_centers, z_vals,string(interpolation_times[j]),p)
end 
gif(anim, "animation4.gif", fps=20)
