

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
    
    mu = log.(abs.(phi).^(1/N) ./ abs.(1 .-phi) .+ 1e-6) .- 2*chi*phi .+ psi + sigma*u
    mucat = log.(abs.(phicat) .+ 1e-6) + u
    muan = log.(abs.(phian) .+ 1e-6) - u

    # calculates all fluxes for each site
    for i in 1:nx, j in 1:ny
        ip1 = (i==nx) ? 1 : i+1
        jp1 = (j==ny) ? 1 : j+1

        Jx = - (phi[i,j] + phi[ip1,j])/2 * (mu[ip1,j] - mu[i,j])/dx
        Jy = - (phi[i,j] + phi[i,jp1])/2 * (mu[i,jp1] - mu[i,j])/dx

        Jcatx = - (phicat[i,j] + phicat[ip1,j])/2 * (mucat[ip1,j] - mucat[i,j])/dx
        Jcaty = - (phicat[i,j] + phicat[i,jp1])/2 * (mucat[ip1,j] - mucat[i,j])/dx

        Janx = - (phian[i,j] + phian[ip1,j])/2 * (muan[ip1,j] - muan[i,j])/dx
        Jany = - (phian[i,j] + phian[i,jp1])/2 * (muan[ip1,j] - muan[i,j])/dx

        # figure out fluxes for u and psi
    end

    # determines dz for all variables
    for i in 1:nx, j in 1:ny
        im1 = (i==1) ? nx : i-1
        jm1 = (j==1) ? ny : j-1

        dz[i,j,1] = -(Jx[i,j]-Jx[im1,j])/dx - (Jy[i,j]-Jy[i,jm1])/dy
        dz[i,j,2] = -(Jcatx[i,j]-Jcat[im1,j])/dx - (Jcaty[i,j]-Jcaty[i,jm1])/dy
        dz[i,j,3] = -(Janx[i,j]-Jan[im1,j])/dx - (Jany[i,j]-Jany[i,jm1])/dy
    end
end