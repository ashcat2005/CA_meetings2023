#=-------------------------------------------------------------------
Find out more info on https://github.com/JackNarvaez/SPH.git
By Narváez J.
-------------------------------------------------------------------=#

using DelimitedFiles

function Euler_Cromer(T, dt, pos, vel, a, N, k, n, lmbda, nu, m, h, Acceleration, Kernel, Gradient_Kernel, ioPos, ioRho)
    #=---------------------------------------------------------------
    Euler Cromer method to integrate the system over time.
    -----------------------------------------------------------------
    Arguments:
    T:      Total time steps.
    dt:     Time step
    pos:    Position of particles.
    vel:    Velocity of particles.
    a:      Acceleration of particles.
    N:      Number of particles.
    k:      Pressure constant
    n:      Polytropic index
    lmbda:  Coefficient of static gravity potential
    nu:     Viscosity coefficient
    m:      Mass particle
    h:      Smoothing kernel length
    Acceleration: Acceleration function
    Kernel: Smoothing function
    Gradient_Kernel: Gradient of the smoothing function
    ---------------------------------------------------------------=#
    for i in 1:T-1
        pos.+= dt.*vel
        rho, P, a = Acceleration(pos, vel, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
        vel .+= dt.*a
        writedlm(ioPos, [pos[:, 1]])
        writedlm(ioPos, [pos[:, 2]])
        writedlm(ioRho, [rho])
        
    end
end;

function Verlet_Pos(T, dt, Acceleration, pos, vel, a, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel, ioPos, ioRho)
    #=---------------------------------------------------------------
    Position Verlet method to integrate the system over time.
    -----------------------------------------------------------------
    Arguments:
    T:      Total time steps.
    dt:     Time step
    pos:    Position of particles.
    vel:    Velocity of particles.
    a:      Acceleration of particles.
    N:      Number of particles.
    k:      Pressure constant
    n:      Polytropic index
    lmbda:  Coefficient of static gravity potential
    nu:     Viscosity coefficient
    m:      Mass particle
    h:      Smoothing kernel length
    Acceleration: Acceleration function
    Kernel: Smoothing function
    Gradient_Kernel: Gradient of the smoothing function
    ---------------------------------------------------------------=#
    for i in 1:T-1
        Temp_Pos = pos .+ 0.5*dt.*vel 
        rho, P, a = Acceleration(Temp_Pos, vel, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
        vel .+= dt.*a
        pos = Temp_Pos .+ 0.5*dt.*vel
        writedlm(ioPos, [pos[:, 1]])
        writedlm(ioPos, [pos[:, 2]])
        writedlm(ioRho, [rho])
    end
end;

function Leap_Frog(T, dt, Acceleration, pos, vel, a, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel, ioPos, ioRho)
    #=---------------------------------------------------------------
    Leap Frog method to integrate the system over time.
    -----------------------------------------------------------------
    Arguments:
    T:      Total time steps.
    dt:     Time step
    pos:    Position of particles.
    vel:    Velocity of particles.
    a:      Acceleration of particles.
    N:      Number of particles.
    k:      Pressure constant
    n:      Polytropic index
    lmbda:  Coefficient of static gravity potential
    nu:     Viscosity coefficient
    m:      Mass particle
    h:      Smoothing kernel length
    Acceleration: Acceleration function
    Kernel: Smoothing function
    Gradient_Kernel: Gradient of the smoothing function
    io:     File to save data
    ---------------------------------------------------------------=#
    rho, P, a = Acceleration(pos, vel, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
    Temp_Vel1 = vel
    for i in 1:T-1
        Temp_Vel2 = Temp_Vel1 .+ dt.*a
        pos .+= dt.*Temp_Vel2
        vel = 0.5.*(Temp_Vel1 .+ Temp_Vel2)
        Temp_Vel1 = Temp_Vel2
        rho, P, a = Acceleration(pos, vel, N, k, n, lmbda, nu, m, h, Kernel, Gradient_Kernel)
        writedlm(ioPos, [pos[:, 1]])
        writedlm(ioPos, [pos[:, 2]])
        writedlm(ioRho, [rho])
    end
end;
