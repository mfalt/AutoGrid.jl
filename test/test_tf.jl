using AutoGrid

func = function (w)
    s = im*w
    tf1 = 0.1*exp(-s)/(s^2 + 1e-4*s +(0.1)^2)
    tf2 = (10^2/(s^2 + s +(10)^2))*(10^2/(s^2 + s +(100)^2))
    return (abs.(tf1), atan.(imag.(tf1),real.(tf1)), abs.(tf2), atan.(imag.(tf2),real.(tf2)))
end
lims = (1e-2, 1e3)
xscale = log10
xscale_inv = exp10
funcConfigs = (FuncConfig(log10, find_peaks=true), WrappingFuncConfig((0.0, 2pi), find_peaks=false), FuncConfig(log10, find_peaks=true), WrappingFuncConfig((0.0, 2pi), find_peaks=false))
agc = AutoGridConfig(func, lims, xscale, xscale_inv, funcConfigs, abstol_dd=1e-2, reltol_dd=1.0, min_eps=1e-4)

values, grid = auto_grid(agc);
@show length(grid)
using Plots
plot(grid, values[1], m=:o, layout=(2,1), subplot=1, scale=:log10, yscale=:log10)
plot!(grid, values[3], m=:o, layout=(2,1), subplot=1, scale=:log10, yscale=:log10)
plot!(grid, values[2]*180/pi, m=:o, xscale=:log10, subplot=2)
plot!(grid, values[4]*180/pi, m=:o, xscale=:log10, subplot=2)

grid2 = 10 .^LinRange(-2,3,length(grid))
plot!(grid2, 0.01.*getindex.(func.(grid2), 1), l=(:green, :dash), subplot=1, m=:o)

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