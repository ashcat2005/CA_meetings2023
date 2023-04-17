#=-------------------------------------------------------------------
Find out more info on https://github.com/JackNarvaez/SPH.git
By Narváez J.
-------------------------------------------------------------------=#

include("ToyStar.jl")
include("Integrator.jl")
include("AuxiliaryFunctions.jl")

# Simulation parameters
seed = parse(UInt64, ARGS[1])  # Random seed
t0 = parse(Float64, ARGS[2])   # Initial time
tEnd = parse(Float64, ARGS[3]) # Final time
dt = parse(Float64, ARGS[4])   # Timestep
d  = parse(UInt8, ARGS[5])     # dimensions
M  = parse(Float64, ARGS[6])   # Star mass
R  = parse(Float64, ARGS[7])   # Star radius
k  = parse(Float64, ARGS[8])   # Pressure constant
n  = parse(Float64, ARGS[9])   # Polytropic index
nu = parse(Float64, ARGS[10])  # Viscosity coefficient
N  = parse(UInt64, ARGS[11])   # Number of Particles

m  = M/N    # Particle mass

h = 0.04/sqrt(N/1000) # Smoothing length
T = Int((tEnd-t0)/dt) # Time Steps
lmbda = Coeff_static_grav_potential(k, n, M, R)

rho = zeros(N) # Density
P = zeros(N) # Pressure
vel = zeros(N, 2) # Velocity
pos = zeros(N, 2) # Position
a = zeros(N, 2) # Acceleration

pos = Init_Dis(N, R, 0, 0, seed)

#Files
ioPos = open("./Files/StarPos.txt", "w");
ioRho = open("./Files/StarRho.txt", "w");
writedlm(ioPos, [pos[:, 1]])
writedlm(ioPos, [pos[:, 2]])
writedlm(ioRho, [rho])

Leap_Frog(T, dt, Acceleration, pos, vel, a, N, k, n, lmbda, nu, m, h, Gaussian_Kernel, Gradient_Gaussian_Kernel, ioPos, ioRho)

close(ioPos)
close(ioRho)
