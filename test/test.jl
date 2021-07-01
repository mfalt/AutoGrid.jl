
using AutoGrid

# global count = 0
func = function (x)
    # global count
    # count += 1
    return (abs(x)<pi ? sin(x) : x > 0 ? -(pi-x)^2 : (-pi-x)^2, abs(x-0.15)+2)
end
lims = (-3, 3)
xscale = identity
xscale_inv = identity
funcConfigs = (FuncConfig(identity, find_peaks=true), FuncConfig(identity, find_peaks=true))
agc = AutoGridConfig(func, lims, xscale, xscale_inv, funcConfigs, min_eps=1e-6)

values, grid = auto_grid(agc);
@show length(grid)
using Plots
plot(grid, values[1], m=:o)
plot!(grid, values[2], m=:o)


values, grid = auto_grid(agc)
using BenchmarkTools
@btime values, grid = auto_grid($agc);

function bench1(grid)
    func.(grid)
end
@btime bench1(grid)

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
# 521.662 μs (6831 allocations: 256.01 KiB) Length: 171

# Bench 1:
# 1.836 μs (4 allocations: 2.88 KiB)
# Bench 2:
# 8.815 μs (350 allocations: 16.31 KiB)