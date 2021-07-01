module AutoGrid
export AutoGridConfig, AbstractFuncConfig, FuncConfig, SavingFuncConfig, WrappingFuncConfig, auto_grid

import LinearAlgebra: dot

abstract type AbstractAutoGridConfig end
abstract type AbstractFuncConfig{saving} end

nothingfun(x) = nothing

struct AutoGridConfig{N, FloatT<:Number, FuncT<:Function, XscaleT<:Function, XscaleInvT<:Function, FuncConfigs<:Tuple{Vararg{<:AbstractFuncConfig,N}}, LT<:Number} <: AbstractAutoGridConfig
    func::FuncT
    lims::NTuple{2,LT}
    xscale::XscaleT
    xscale_inv::XscaleInvT
    funcConfigs::FuncConfigs
    start_gridpoints::Int
    maxgridpoints::Int
    min_eps::FloatT
    reltol_dd::FloatT
    abstol_dd::FloatT
end

""" Standard config """
struct FuncConfig{saving, YscaleT<:Function, Diff2T, Diff3T} <: AbstractFuncConfig{saving}
    yscale::YscaleT
    find_peaks::Bool
    diff2::Diff2T
    diff3::Diff3T
    FuncConfig{saving}(yscale; find_peaks=true, diff2=normdiff2, diff3=normdiffdiff2) where saving =
        new{saving, typeof(yscale), typeof(diff2), typeof(diff3)}(yscale, find_peaks, diff2, diff3)
end
FuncConfig(args...; kwargs...) = FuncConfig{true}(args...; kwargs...)


""" Only saves the corresponding data, not used to computing grid """
struct SavingFuncConfig{saving} <: AbstractFuncConfig{saving}
    function SavingFuncConfig{saving}() where saving
        saving === true || throw(ArgumentError("SavingFuncConfig must have \"saving===true\" "))
        return new{true}()
    end
    SavingFuncConfig() = SavingFuncConfig{true}()
end
function Base.getproperty(sfc::SavingFuncConfig, sym)
    if sym === :yscale
        return nothingfun
    elseif sym === :find_peaks
        return false
    else
        return getfield(sfc, sym)
    end
end

struct WrappingFuncConfig{saving, YscaleT<:Function, Diff2T, Diff3T, T<:Number} <: AbstractFuncConfig{saving}
    yscale::YscaleT
    find_peaks::Bool
    wraplims::Tuple{T,T}
    diff2::Diff2T
    diff3::Diff3T
    WrappingFuncConfig{saving}(yscale, wraplims; find_peaks=true, diff2=normdiff2, diff3=normdiffdiff2) where saving=
        new{saving, typeof(yscale), typeof(wraplims), typeof(diff2), typeof(diff3)}(yscale, find_peaks, wraplims, diff2, diff3)
end
WrappingFuncConfig(args...; kwargs...) = WrappingFuncConfig{true}(args...; kwargs...)

AutoGridConfig(func::FuncT, lims::NTuple{2,LT}, xscale::XscaleT, xscale_inv::XscaleInvT, funcConfigs::FuncConfigs;
                start_gridpoints=30, maxgridpoints=2000,
                min_eps::FloatT=(xscale(lims[2])-xscale_inv(lims[1]))/10000,
                reltol_dd::FloatT=0.05, abstol_dd::FloatT=1e-2) where {N, FloatT<:Number, FuncT<:Function, XscaleT<:Function, XscaleInvT<:Function, FuncConfigs<:Tuple{Vararg{<:AbstractFuncConfig,N}}, LT<:Number} =
    AutoGridConfig{N, FloatT, FuncT, XscaleT, XscaleInvT, FuncConfigs, LT}(func, lims, xscale, xscale_inv, funcConfigs, start_gridpoints, maxgridpoints, min_eps, reltol_dd, abstol_dd)
"""
values, grid = auto_grid(func::Function, N, lims::Tuple, xscale=(identity, identity), FINCCONFIGSDEF; start_gridpoints, maxgridpoints, min_eps, reltol_dd, abstol_dd)
values, grid = auto_grid(agc::AutoGridConfig)
    
    Arguments:
    Compute a `grid` that tries to capture all features of the function `func` in the range `lim=(xmin,xmax)`.
    `func`: `Function` so that that returns a `Tuple` of values, e.g. `func(x) = (v1,v2,...,vN)` where each `vi` is a scalar or AbstractArray.
        note that the behavior for non scalar vi might be unintuitive.
    The algorithm will try to produce a grid so that features are captured in all of the values, but widely varying scaled might cause problems.
    `xscale`: `Tuple` containing a monotone function and its inverse corresponding to the default axis.
        Example: If you want to produce a plot with logarithmic x-axis, then you should set `xscale=(log10, exp10)`.
    `yscales`: `Tuple` containing one monotone function for each of the values `(v1,v2,...,vN)` corresponding to the default axis.
    
    Kwargs:
    `start_gridpoints`: The number of initial grid points.
    `maxgridpoints`: A warning will be thrown if this number of gridpoints are reached.
    `min_eps`: Lower limit on how small steps will be taken (in `xscale`).
        Example: with `xscale=(log10, exp10)` and `eps=1e-2`then `log10(grid[k+1])-log10(grid[k]) > (1e-2)/2`, i.e `grid[k+1] > 1.012grid[k]`.
    `reltol_dd = 0.05`, `abstol_dd = 0.01`: The criteria for further refinement is
        `norm(2*y2-y1-y3) < max(reltol_dd*max(norm(y2-y1), norm(y3-y2)), abstol_dd)`,
        where `y1,y2,y3` (in yscale) are 3 equally spaced points (in xscale).

    Output: 
    `values`: `Tuple` with one vector per output value of `func`.
    `grid`: `Vector` with the grid corresponding to the values, i.e `func(grid[i])[j] = values[j][i]`

    Example: The following code computes a grid and plots the functions: abs(sinc(x)) and cos(x)+1
        in lin-log and log-log plots respectively.
    
    func = x -> (abs(sinc(x)),cos(x)+1)
    lims = (0.1,7.0)
    xscale = (log10, exp10)
    yscales = (identity, log10)
    y, x = auto_grid(func, lims, xscale, yscales, min_eps=1e-3)
    
    plot(x, y[1], m=:o; xscale=:log10, layout=(2,1), yscale=:identity, subplot=1, lab="sinc(x)", size=(800,600))
    plot!(x, y[2], m=:o, xscale=:log10, subplot=2, yscale=:log10, lab="cos(x)+1", ylims=(1e-5,10))

"""
# function auto_grid(func, lims::Tuple, xscale::Tuple, yscales::Tuple;
#     start_gridpoints=30, maxgridpoints=2000,
#     min_eps=(xscale[1](lims[2])-xscale[1](lims[1]))/10000,
#     reltol_dd=0.05, abstol_dd=1e-2)
function auto_grid(agc::AbstractAutoGridConfig)
    
    # Linearly spaced grid in xscale
    init_grid = LinRange(agc.xscale(agc.lims[1]), agc.xscale(agc.lims[2]), agc.start_gridpoints)
    
    # Current count of number mindpoints
    num_gridpoints = agc.start_gridpoints

    func = agc.func
    funcConfigs = agc.funcConfigs
    # Current left point
    x1 = init_grid[1]
    x1_id = agc.xscale_inv(x1)
    
    y1 = apply_config_scales(funcConfigs, func(x1_id))
    prev_val_ref = Ref{Ref{typeof(y1)}}()

    y1_id = func(x1_id)

    # The full set of gridpoints
    grid = [x1_id,]
    # Tuple with list of all values (Faster than list of tuples + reshaping)
    #values = ntuple(i -> (saving ? [y1_id[i],] : nothing), length(yscales)) # Type unstable, but low performance cost, but affects output
    values = tuple_of_vectors(y1_id, funcConfigs)

    for (i,x2) in enumerate(init_grid)
        i == 1 && continue # The first interval starts at 1 -> 2
        # Needed when looking for local optima
        next_val_ref = Ref(Ref(apply_config_scales(funcConfigs,func(agc.xscale_inv(init_grid[i])+agc.min_eps/10))))
        # else
        #     next_val_ref = Ref{Ref{typeof(y1)}}()
        # end
        # Scale back to identity
        x2_id = agc.xscale_inv(x2)
        y2_id = func(x2_id)
        y2 = apply_config_scales(funcConfigs, y2_id)
        # Refine (if nessesary) section (x1,x2)
        num_gridpoints = refine_grid!(values, grid, x1, x2, y1, y2, y1_id, y2_id, prev_val_ref, next_val_ref, true, num_gridpoints, agc)
        # We are now done with [x1,x2]
        # Continue to next segment
        x1 = x2
        y1 = y2
        y1_id = y2_id
    end
    return values, grid
end

apply_scale(fc::AbstractFuncConfig, y_id) = fc.yscale(y_id)
apply_scale(fc::SavingFuncConfig, y_id) = nothing

function apply_config_scales(funcConfigs, y_id)
    map(apply_scale, funcConfigs, y_id)
end

# Given range (x1, x2) potentially add more gridpoints. x1 is assumed to be already added, x2 will be added
function refine_grid!(values, grid, x1, x2, y1, y2, y1_id, y2_id, prev_val_ref, next_val_ref, close_to_next_val_ref, num_gridpoints, agc::AbstractAutoGridConfig)
    # In scaled grid
    xm = (x1+x2)/2

    # Scaled back
    xm_id = agc.xscale_inv(xm)
    ym_id = agc.func(xm_id)
    ym = apply_config_scales(agc.funcConfigs, ym_id)
    # println("$x1, $xm, $x2")
    # isassigned(prev_val_ref) && isassigned(next_val_ref) && println("$(prev_val_ref[][]), $y1, $ym, $y2, $(next_val_ref[][])")
    (num_gridpoints >= agc.maxgridpoints) && @warn "Maximum number of gridpoints reached in refine_grid! at $xm_id, no further refinement will be made. Increase maxgridpoints to get better accuracy." maxlog=1
    #min_eps in scaled version, abs to avoid assuming monotonly increasing scale
    if (abs(x2 - x1) >= agc.min_eps) && (num_gridpoints < agc.maxgridpoints)
        refine_left, refine_right = should_refine(y1, ym, y2, prev_val_ref, next_val_ref, agc)
    else
        refine_left = refine_right = false
    end

    if refine_left
        # println("Ref left")
        refnext = Ref(Ref(y2))
        num_gridpoints = refine_grid!(values, grid, x1, xm, y1, ym, y1_id, ym_id, prev_val_ref, refnext, false, num_gridpoints, agc)
    else # We are done up until xm
        push!(grid, xm_id)
        map(save, values, ym_id, agc.funcConfigs)
        num_gridpoints += 1
    end

    if refine_right
        # println("Ref right")
        if !refine_left && !close_to_next_val_ref # Only occurs when looking for local optima
            # Could restrict further to only when right point if far away
            # Some extra evaluations here to avoid problems with right point being far away
            refnext = Ref(Ref(apply_config_scales(agc.funcConfigs,agc.func(agc.xscale_inv(x2)+agc.min_eps/10))))
            close_to_next_val_ref = true
        else
            refnext = next_val_ref
        end
        num_gridpoints = refine_grid!(values, grid, xm, x2, ym, y2, ym_id, y2_id, prev_val_ref, refnext, close_to_next_val_ref, num_gridpoints, agc)
    else # We are done up until x2
        push!(grid, agc.xscale_inv(x2))
        map(save, values, y2_id, agc.funcConfigs)
        num_gridpoints += 1
        # The last point before x2 was xm
        prev_val_ref[] = Ref(ym)
    end

    return num_gridpoints
end



# TODO We can scale when saving instead of recomputing scaling
# Svectors should also be faster in general
function should_refine(v1, v2, v3, prev_val_ref, next_val_ref, agc::AbstractAutoGridConfig)
    
    not_linear = !all(is_almost_linear.(v1, v2, v3, agc.funcConfigs, Ref(agc))) # Ref so struct is scalar
    not_linear && return true, true

    prev_val = isassigned(prev_val_ref) ? prev_val_ref[][] : nothing
    next_val = isassigned(next_val_ref) ? next_val_ref[][] : nothing
    pot_peaks = potential_peak.(v1, v2, v3, prev_val, next_val, agc.funcConfigs, Ref(agc))
    left_peaks = any(getindex.(pot_peaks, 1))
    right_peaks = any(getindex.(pot_peaks, 2))
    return left_peaks, right_peaks
end

function potential_peak(y1, ym, y2, prev_val, next_val, fc::AbstractFuncConfig, agc::AbstractAutoGridConfig)
    if fc.find_peaks
        #c1 = dot(ym-y1, y2-ym) < 0
        c1 = dotdiff(ym,y1,y2,ym) < 0
        c1 && return true, true
        # c2 = (prev_val !== nothing && dot(y1-prev_val, ym-y1) < 0)
        # c3 = (next_val !== nothing && dot(y2-ym, next_val-y2) < 0)
        c2 = (prev_val !== nothing && dotdiff(y1,prev_val,ym,y1) < 0)
        c3 = (next_val !== nothing && dotdiff(y2,ym,next_val,y2) < 0)
        return c2, c3
    else
        return false, false
    end
end

@inline dotdiff(v1::T,v2::T,v3::T,v4::T) where T<:Real = (v1-v2)*(v3-v4)
# TODO UNSAFE
@inline function dotdiff(v1::T,v2::T,v3::T,v4::T) where T
    s = zero(dot(first(v1)-first(v2), first(v3)-first(v4)))
    for i in eachindex(v1)
        @inbounds s += dot(v1[i]-v2[i], v3[i]-v4[i])
    end
    return s
end

@inline function is_almost_linear(y1, ym, y2, fc::AbstractFuncConfig, agc::AbstractAutoGridConfig)
    # We assume that x2-x1 \approx  x3-x2, so we need to check that ym-y1 approx y2-ym
    # Essentially low second derivative compared to derivative, so that linear approximation is good.
    # Second argument to avoid too small steps when derivatives are small
    lhs = normdiffdiff2(y1,ym,y2)
    rhs = max((agc.reltol_dd)^2*max(normdiff2(ym,y1), normdiff2(y2,ym)), (agc.abstol_dd)^2)
    # println("$y1, $ym, $y2")
    # println("$(lhs < rhs) lhs/rhs: $lhs, $rhs")
    return lhs < rhs
end
@inline is_almost_linear(y1, ym, y2, fc::SavingFuncConfig, agc::AbstractAutoGridConfig) = true

# This will probably never (almost) be called
@inline normdiff2(x::Number,y::Number) = abs2(x-y)

function normdiff2(x::A,y::A) where {T<:Number, N, A <: Union{NTuple{N,T},AbstractArray{T}}}
    length(x) == length(y) || throw(DimensionMismatch())
    sum = abs2(x[1]-y[1])
    @inbounds for i in 2:length(x)
        sum += abs2(x[i]-y[i])
    end
    return sum
end
function normdiff2(x,y)
    length(x) == length(y) || throw(DimensionMismatch())
    sum = normdiff2(x[1],y[1])
    @inbounds for i in 2:length(x)
        sum += normdiff2(x[i],y[i])
    end
    return sum
end
# @inline normdiff(x,y) = sqrt(normdiff2(x,y))


# This will probably never(almost) be called
@inline normdiffdiff2(x::Number,y::Number,z::Number) = abs2(y+y-x-z)

function normdiffdiff2(x::A,y::A,z::A) where {T<:Number, N, A <: Union{NTuple{N,T},AbstractArray{T}}}
    length(x) == length(y) || throw(DimensionMismatch())
    sum = abs2(y[1]+y[1]-x[1]-z[1])
    @inbounds for i in 2:length(x)
        sum += abs2(y[i]+y[i]-x[i]-z[i])
    end
    return sum
end
function normdiffdiff2(x,y,z)
    length(x) == length(y) || throw(DimensionMismatch())
    sum = normdiffdiff2(x[1],y[1],z[1])
    @inbounds for i in 2:length(x)
        sum += normdiffdiff2(x[i],y[i],z[i])
    end
    return sum
end
# @inline normdiffdiff(x,y,z) = sqrt(normdiffdiff2(x,y,z))

# @generated function apply_tuple3(fs::Tuple{Vararg{<:Any,N}}, x::Tuple{Vararg{<:Any,N}}) where N
#     vec = []
#     for i in 1:N
#         push!(vec, :(fs[$i].(x[$i])))
#     end
#     :(Core.tuple($(vec...)))
# end

@generated function tuple_of_vectors(y1_id::Tuple{Vararg{<:Any,N}}, funcConfigs::Tuple{Vararg{<:AbstractFuncConfig,N}}) where {N}
    vec = []
    for i in 1:N
        # println("Generating")
        # println(funcConfigs)
        if funcConfigs.types[i] <: AbstractFuncConfig{true}
            push!(vec, :([y1_id[$i],]))
        else
            push!(vec, :(nothing))
        end
    end
    :(Core.tuple($(vec...),))
end

save(values, ym_id, funcConfig::AbstractFuncConfig{false}) = nothing
function save(values, ym_id, funcConfig::AbstractFuncConfig{true})
    push!(values, ym_id)
    return nothing
end

end



# function auto_grid_alloc(agc::AbstractAutoGridConfig)
#     func = agc.func
#     yscales = agc.yscales

#     # Linearly spaced grid in xscale
#     init_grid = LinRange(agc.xscale(agc.lims[1]), agc.xscale(agc.lims[2]), agc.start_gridpoints)

#     # Current count of number mindpoints
#     num_gridpoints = agc.start_gridpoints

#     # Current left point
#     x1 = init_grid[1]
#     x1_id = agc.xscale_inv(x1)
#     y1 = apply_tuple3(yscales, func(x1_id))
#     y1_id = func(x1_id)


#     # The full set of gridpoints
#     grid = [x1_id,]
#     # Tuple with list of all values (Faster than list of tuples + reshaping)
#     #values = ntuple(i -> [y1_id[i],], length(yscales)) # Type unstable, but low performance cost, but affects output
#     values = tuple_of_vectors(y1_id)

#     # Compute some values on initial grid
#     init_grid_vals = similar.(values, agc.start_gridpoints)

#     for (i,x) in enumerate(init_grid)
#         x_id = agc.xscale_inv(x)
#         y_id = func(x_id)
#         y = apply_tuple3(yscales, y_id)
#         setindex!.(init_grid_vals, y, i)
#     end
#     @show abs.(init_grid_vals[1])
#     @show real.(diff(init_grid_vals[1]))
#     grid_norms = [norm.(init_grid_vals[i])./(init_grid[2]-init_grid[1]) for i in 1:length(init_grid_vals)]
#     grid_diffs_norms = [sqrt.(normdiff2.(init_grid_vals[i][1:end-1],init_grid_vals[i][2:end]))./(init_grid[2]-init_grid[1]) for i in 1:length(init_grid_vals)]
#     println(grid_norms)
#     println(grid_diffs_norms)
#     # println(diffs)

#     for x2 in init_grid[2:end]
#         # Scale back to identity
#         x2_id = agc.xscale_inv(x2)
#         y2_id = func(x2_id)
#         y2 = apply_tuple3(yscales, y2_id)
#         # Refine (if nessesary) section (x1,x2)
#         num_gridpoints = refine_grid!(values, grid, x1, x2, y1, y2, y1_id, y2_id, num_gridpoints, agc)
#         # We are now done with [x1,x2]
#         # Continue to next segment
#         x1 = x2
#         y1 = y2
#         y1_id = y2_id
#     end
#     return values, grid
# end
