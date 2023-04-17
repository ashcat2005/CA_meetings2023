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
k  = parse(Float64, ARGS[6])   # Pressure constant
n  = parse(Float64, ARGS[7])   # Polytropic index
nu = parse(Float64, ARGS[8])   # Viscosity coefficient

# Star 1
M1 = parse(Float64, ARGS[9])   # Star mass
R1 = parse(Float64, ARGS[10])  # Star radius
N1 = parse(UInt64, ARGS[11])   # Number of Particles
m1 = M1/N1  # Particle mass

# Star 2
M2 = parse(Float64, ARGS[12])  # Star mass
R2 = parse(Float64, ARGS[13])  # Star radius
N2 = parse(UInt64, ARGS[14])   # Number of Particles
m2 = M2/N2  # Particle mass

N = (N1+N2)
M = 0.5*(M1+M2)
m = 0.5*(m1+m2)
R = 0.5*(R1+R2)

h = 0.04/sqrt(0.5*N/1000) # Smoothing length
T = Int((tEnd-t0)/dt) # Time Steps
lmbda = Coeff_static_grav_potential(k, n, 2*M, R)

rho = zeros(N) # Density
P = zeros(N) # Pressure
vel = zeros(N, 2) # Velocity
pos = zeros(N, 2) # Position
a = zeros(N, 2) # Acceleration

pos[1:N1, :] = Init_Dis(N1, R1, -1, -1, seed)
pos[N1+1:end, :] = Init_Dis(N2, R2, 1, 1, seed);

#Files
ioPos = open("./Files/StarPos.txt", "w");
ioRho = open("./Files/StarRho.txt", "w");
writedlm(ioPos, [pos[:, 1]])
writedlm(ioPos, [pos[:, 2]])
writedlm(ioRho, [rho])

Leap_Frog(T, dt, Acceleration, pos, vel, a, N, k, n, lmbda, nu, m, h, Gaussian_Kernel, Gradient_Gaussian_Kernel, ioPos, ioRho)

close(ioPos)
close(ioRho)

record = 2
pmin = 0
pmax = 15

ioPos = open("./Files/StarPos.txt", "r");
ioRho = open("./Files/StarRho.txt", "r");
animation(T, record, dt, ioPos, ioRho, "Density", pmin, pmax, 2.5*R, "2Star")
close(ioPos)
close(ioRho)
ioPos = open("./Files/StarPos.txt", "r");
ioRho = open("./Files/StarRho.txt", "r");
animation2(T, record, dt, ioPos, pmin, pmax, 2.5*R, "2StarColor", N1)
close(ioPos)
close(ioRho)
