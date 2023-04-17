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
R  = parse(Float64, ARGS[4])   # Star radius
N  = parse(UInt64, ARGS[5])   # Number of Particles

pmin = 0.0 
pmax = 3.0

h = 0.04/sqrt(N/1000) # Smoothing length
T = Int((tEnd-t0)/dt) # Time Steps
record = 2

ioPos = open("./Files/StarPos.txt", "r");
ioRho = open("./Files/StarRho.txt", "r");

animation(T, record, dt, ioPos, ioRho, "Density", pmin, pmax, R, "1Star")

close(ioPos)
close(ioRho)
