#=-------------------------------------------------------------------
Find out more info on https://github.com/JackNarvaez/SPH.git
By Narváez J.
-------------------------------------------------------------------=#
using Plots
using LaTeXStrings
using KissSmoothing # Smooth data from particles' values
using DelimitedFiles

function animation(T, record, dt, ioPos, ioProp, Property, pmin, pmax, R, Name)
    #=---------------------------------------------------------------
    Generates a .mp4 video of the evolution of a system. It takes the
    values from the files ioPos and ioProp. The scale of the cmap is
    related to the magnitude of Property.
    ---------------------------------------------------------------=#
    anim = @animate for n in 1:record:T
        x = ConvertStringtoFloat(split(readline(ioPos), "\t"), N)
        y = ConvertStringtoFloat(split(readline(ioPos), "\t"), N)
        prop = ConvertStringtoFloat(split(readline(ioProp), "\t"), N)
        plot(x, y, seriestype=:scatter, aspect_ratio=:equal, cmap=:thermal, markerstrokewidth=0,
            marker_z = prop, ms=1.5, ma=0.5, xlim=(-4/3*R, 4/3*R), ylim=(-4/3*R, 4/3*R), legend=false, grid=false, 
            colorbar=true, colorbar_title="\n"*Property, clims=(pmin,pmax), dpi=300)
        xlabel!(L"x")
        ylabel!(L"y")
        annotate!(0, 4/3*R, text("t = $(round(n*dt;digits=1))",8))
    end
    gif(anim, Name*".mp4", fps = 10)
end;

function animation2(T, record, dt, ioPos, pmin, pmax, R, Name, N1)
    #=---------------------------------------------------------------
    Generates a .mp4 video of the evolution of a 2 body system. It 
    shows in red and blue the particles of first and second body, 
    respectively.
    ---------------------------------------------------------------=#
    anim = @animate for n in 1:record:T
        x = ConvertStringtoFloat(split(readline(ioPos), "\t"), N)
        y = ConvertStringtoFloat(split(readline(ioPos), "\t"), N)
        plot(x[1:N1], y[1:N1], seriestype=:scatter, aspect_ratio=:equal, mc=:red, markerstrokewidth=0,
            ms=1.5, ma=0.5, xlim=(-4/3*R, 4/3*R), ylim=(-4/3*R, 4/3*R), legend=false, grid=false, 
            dpi=300)
        plot!(x[N1+1:end], y[N1+1:end], seriestype=:scatter, aspect_ratio=:equal, mc=:blue, markerstrokewidth=0,
            ms=1.5, ma=0.5, xlim=(-4/3*R, 4/3*R), ylim=(-4/3*R, 4/3*R), legend=false, grid=false, 
            dpi=300)
        xlabel!(L"x")
        ylabel!(L"y")
        annotate!(0, 4/3*R, text("t = $(round(n*dt;digits=1))",8))
    end
    gif(anim, Name*".mp4", fps = 10)
end;

function PlotDensity(radius, dens, dens_Theo, R, PropMax, dots=true)
    #=---------------------------------------------------------------
    Plots the mass density vs the radius of a Toy Star. 
    ---------------------------------------------------------------=#
    P_r = sortperm(radius)
    plot(radius[P_r], dens[P_r],  seriestype=:scatter, ms=1.2, markerstrokewidth=0, color=:crimson, xlim=(0, R), ylim=(0, PropMax), labels=L"\rho_{SPH}(r)", grid=false, dpi=300)
    plot!(radius[P_r], dens_Theo[P_r], color=:blue, labels=L"\rho_{theo}(r)", grid=false, legend=true, dpi=300)
    xlabel!(L"r")
    ylabel!(L"\rho")
    savefig("Density.png")
end;

function ConvertStringtoFloat(String, N)
    #=---------------------------------------------------------------
    Converts an array with N entries from String to Float64 type.
    ---------------------------------------------------------------=#
    Array = zeros(N)
    for i in 1:N
        Array[i] = parse(Float64, String[i])
    end
    return Array
end;
