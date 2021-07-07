using AutoGrid

# global count = 0
func = function (x)
    # global count
    # count += 1
    return (abs(x)<pi ? sin(x) : x > 0 ? -(pi-x)^2 : (-pi-x)^2, abs(x-0.15)+2, mod(0.3*sin(x)+x, 1.5), abs(0.1/((im*exp10(x))^2 + 0.001*(0.1)*(im*exp10(x)) +(0.1)^2) ), )
    # return (abs(x-0.15)+2, mod(x, 1.5))
    #return (mod(x, 1.5),)
    # return (sin(x),)
end
lims = (-3, 3)
xscale = identity
xscale_inv = identity
funcConfigs = (FuncConfig(identity, find_peaks=false), FuncConfig(identity, find_peaks=false), WrappingFuncConfig((0.0, 1.5), find_peaks=false), FuncConfig(log10, find_peaks=true))
# funcConfigs = (FuncConfig(identity, find_peaks=false), WrappingFuncConfig((0.0, 1.5), find_peaks=false))
#funcConfigs = (WrappingFuncConfig(identity, (0.0, 1.5), find_peaks=false),)
agc = AutoGridConfig(func, lims, xscale, xscale_inv, funcConfigs, maxgridpoints=2000, abstol_dd=1e-2, reltol_dd=1.0, min_eps=1e-8)

values, grid = auto_grid(agc);
@show length(grid)
using Plots
plot(grid, values[1], m=:o)
plot!(grid, values[2], m=:o)
plot!(grid, values[3], m=:o)
plot!(grid, log10.(values[4]), m=:o)
plot!(grid, getindex.(func.(grid), 3), l=(:green, :dash))

values, grid = auto_grid(agc)
using BenchmarkTools
@btime values, grid = auto_grid($agc);

function bench1(grid)
    func.(grid)
end
@btime bench1(grid);

function bench2(grid)
    vec = [func(grid[1]),]
    for i in 2:length(grid)
        push!(vec, func(grid[i]))
    end
    return vec
end
@btime bench2(grid)
# No peaks
# 180.647 μs (2412 allocations: 91.73 KiB) Length: 85
# Both peaks
# 35.296 μs (61 allocations: 20.11 KiB) Length: 196
# 68.038 μs (59 allocations: 29.14 KiB) Length: 376

# Bench 1:
# 1.836 μs (4 allocations: 2.88 KiB)
# Bench 2:
# 8.815 μs (350 allocations: 16.31 KiB)