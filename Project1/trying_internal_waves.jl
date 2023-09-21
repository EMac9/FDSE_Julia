using Oceananigans, Oceananigans.Grids, Plots, Printf
Nx = 128 # resolution
Lx = 2π  # domain extent

# Non-dimensional internal wave parameters
m = 16      # vertical wavenumber
k = 8       # horizontal wavenumber
N = 1       # buoyancy frequency
f = 0.2     # inertial frequency

ω² = (N^2 * k^2 + f^2 * m^2) / (k^2 + m^2)

# and thus
ω = sqrt(ω²)

U = k * ω   / (ω^2 - f^2)
V = k * f   / (ω^2 - f^2)
W = m * ω   / (ω^2 - N^2)
B = m * N^2 / (ω^2 - N^2)

# Some Gaussian parameters
A, x₀, z₀, δ = 1e-9, Lx/2, -Lx/2, Lx/15

# A Gaussian envelope
a(x, z) = A * exp( -( (x - x₀)^2 + (z - z₀)^2 ) / 2δ^2 )

u₀(x, y, z) = a(x, z) * U * cos(k*x + m*z)
v₀(x, y, z) = a(x, z) * V * sin(k*x + m*z)
w₀(x, y, z) = a(x, z) * W * cos(k*x + m*z)
b₀(x, y, z) = a(x, z) * B * sin(k*x + m*z) + N^2 * z

model = IncompressibleModel(    grid = RegularCartesianGrid(size=(Nx, 1, Nx), extent=(Lx, Lx, Lx)),
                             closure = ConstantIsotropicDiffusivity(ν=1e-6, κ=1e-6),
                            coriolis = FPlane(f=f),
                             tracers = :b,
                            buoyancy = BuoyancyTracer())

set!(model, u=u₀, v=v₀, w=w₀, b=b₀)


simulation = Simulation(model, Δt = 0.001 * 2π/ω, stop_iteration = 0,
                        progress_frequency = 20)

anim = @animate for i=0:100
    x, z = xnodes(Cell, model.grid)[:], znodes(Face, model.grid)[:]
    w = model.velocities.w

    contourf(x, z, w.data[1:Nx, 1, 1:Nx+1]',
             title=@sprintf("ωt = %.2f", ω*model.clock.time),
             levels=range(-1e-8, stop=1e-8, length=10),
             clims=(-1e-8, 1e-8),
             xlabel="x", ylabel="z",
             xlims=(0, Lx), ylims=(-Lx, 0),
             linewidth=0,
             c=:balance,
             legend=false,
             aspectratio=:equal)

    simulation.stop_iteration += 20
    run!(simulation)
end

mp4(anim, "internal_wave.mp4", fps = 15) # hide