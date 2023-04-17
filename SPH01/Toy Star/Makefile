# -------------------------------------------------------------------
# Find out more info on https://github.com/JackNarvaez/SPH.git
# By Narváez J.
# -------------------------------------------------------------------

#---------- General Simulation Parameters ----------#
seed = 1234 # Random seed
t0 = 0      # Initial time
tEnd = 20   # Final time
dt = 0.04   # Timestep
d  = 2      # dimensions

#-------- Simulation parameters for 1 Star ---------#
M  = 2      # Star mass
R  = 0.75   # Star radius
k  = 0.1    # Pressure constant
n  = 1      # Polytropic index
nu = 1      # Viscosity coefficient
N  = 100   # Number of Particles

#----- Simulation parameters for 2 Stars -----#
#Star 1
M1 = 4      # Star mass
R1 = 0.75   # Star radius
N1 = 400    # Number of Particles

# Star 2
M2 = 1      # star mass
R2 = 0.5      # star radius
N2 = 100   # Number of Particles

all: Main

Main: main.jl ToyStar.jl Integrator.jl AuxiliaryFunctions.jl
	julia $< ${seed} ${t0} ${tEnd} ${dt} ${d} ${M} ${R} ${k} ${n} ${nu} ${N}

1Star: main1Star.jl Density.jl ToyStar.jl Integrator.jl AuxiliaryFunctions.jl
	julia $< ${t0} ${tEnd} ${dt} ${R} ${N};\
	julia Density.jl ${t0} ${tEnd} ${dt} ${N} ${M} ${R} ${k} ${n}

2Star: main2Star.jl ToyStar.jl Integrator.jl AuxiliaryFunctions.jl
	julia $< ${seed} ${t0} ${tEnd} ${dt} ${d} ${k} ${n} ${nu} ${M1} ${R1} ${N1} ${M2} ${R2} ${N2}

clean:
	rm -f *.x *.txt *.png *.out *.gif *.mp4
