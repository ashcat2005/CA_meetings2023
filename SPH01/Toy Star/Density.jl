#=-------------------------------------------------------------------
Find out more info on https://github.com/JackNarvaez/SPH.git
By Narváez J.
-------------------------------------------------------------------=#
include("ToyStar.jl")
include("Integrator.jl")
include("AuxiliaryFunctions.jl")

t0 = parse(Float64, ARGS[1])   # Initial time
tEnd = parse(Float64, ARGS[2]) # Final time
dt = parse(Float64, ARGS[3])   # Timestep
N  = parse(UInt64, ARGS[4])   # Number of Particles
M  = parse(Float64, ARGS[5])   # Star mass
R  = parse(Float64, ARGS[6])   # Star radius
k  = parse(Float64, ARGS[7])   # Pressure constant
n  = parse(Float64, ARGS[8])   # Polytropic index
record = 2

h = 0.04/sqrt(N/1000) # Smoothing length
T = Int((tEnd-t0)/dt) # Time Steps
lmbda = Coeff_static_grav_potential(k, n, M, R)

ioPos = open("./Files/StarPos.txt", "r");
ioRho = open("./Files/StarRho.txt", "r");

pos = readdlm(ioPos, Float64, skipstart=2*T-2)
rho = readdlm(ioRho, Float64, skipstart=T-1)

close(ioPos)
close(ioRho)

pmin = minimum(rho) 
pmax = maximum(rho)

r_final = sqrt.(pos[1, :].^2 + pos[2, :].^2)
dens_Theo = DensityTheo(r_final, R, lmbda, k)

PlotDensity(r_final[:], rho[:], dens_Theo[:], R,  1.3*maximum(dens_Theo[:]))
