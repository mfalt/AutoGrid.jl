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
                reltol_dd::FloatT=0.1, abstol_dd::FloatT=1e-6) where {N, FloatT<:Number, FuncT<:Function, XscaleT<:Function, XscaleInvT<:Function, FuncConfigs<:Tuple{Vararg{<:AbstractFuncConfig,N}}, LT<:Number} =
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
    num_gridpoints = agc.start_gridpoints
    func = agc.func
    funcConfigs = agc.funcConfigs
    
    # Linearly spaced grid in xscale
    init_grid = collect(LinRange(agc.xscale(agc.lims[1]), agc.xscale(agc.lims[2]), agc.start_gridpoints))

    x1, x1_id, y1, y1_id = eval_at_point(agc, init_grid[1])

    values_id = tuple_of_vectors(y1_id, funcConfigs)
    grid_id = [x1_id,]

    xend, xend_id, yend, yend_id = eval_at_point(agc, init_grid[end])
    
    buffer_x = [xend,]
    buffer_x_id = [xend_id,]
    buffer_values = tuple_of_vectors(yend, funcConfigs)
    buffer_values_id = tuple_of_vectors(yend_id, funcConfigs)
    # Fill right buffer with 2,.....,(end-1) in reverse order
    # xend is already added
    for i in (length(init_grid)-1):-1:2
        xi, xi_id, yi, yi_id = eval_at_point(agc, init_grid[i])
        push_buffers!(xi, xi_id, yi, yi_id, buffer_x, buffer_x_id, buffer_values, buffer_values_id, funcConfigs)
    end

    x2, x2_id, y2, y2_id = pop_buffers!(buffer_x, buffer_x_id, buffer_values, buffer_values_id, funcConfigs)
    x3, x3_id, y3, y3_id = pop_buffers!(buffer_x, buffer_x_id, buffer_values, buffer_values_id, funcConfigs)

    # Inside loop we assume x1,x2,x3 are set to current interval
    # x1, y1 should already be in `grid`, `values_id` etc., x2, x3 not
    while true
        stop_refine = num_gridpoints >= agc.maxgridpoints
        stop_refine && @warn "Maximum number of gridpoints reached in refine_grid! at $x2_id, no further refinement will be made. Increase maxgridpoints to get better accuracy." maxlog=1
        # TODO! use differences in x
        if !stop_refine && (max(abs(x2 - x1),abs(x3 - x2)) >= agc.min_eps) && should_refine(y1, y2, y3, x1, x2, x3, agc)
            xm_right = (x2+x3)/2
            xm_left = (x1+x2)/2
            # Rightmost point
            push_buffers!(x3, x3_id, y3, y3_id, buffer_x, buffer_x_id, buffer_values, buffer_values_id, funcConfigs)
            # Second to last point 
            xi, xi_id, yi, yi_id = eval_at_point(agc, xm_right)
            push_buffers!(xi, xi_id, yi, yi_id, buffer_x, buffer_x_id, buffer_values, buffer_values_id, funcConfigs)
            # Old middle point is now last point
            x3, x3_id, y3, y3_id = x2, x2_id, y2, y2_id
            # New middle point
            x2, x2_id, y2, y2_id = eval_at_point(agc, xm_left)
            num_gridpoints += 2
        else
            # Save completed values to output
            push!(grid_id, x2_id)
            map(save, values_id, y2_id, agc.funcConfigs)

            if isempty(buffer_x) # We are done
                push!(grid_id, x3_id)
                map(save, values_id, y3_id, agc.funcConfigs)
                break
            end
            # Prepare for next loop
            x1, x1_id, y1, y1_id = x2, x2_id, y2, y2_id
            x2, x2_id, y2, y2_id = x3, x3_id, y3, y3_id
            x3, x3_id, y3, y3_id = pop_buffers!(buffer_x, buffer_x_id, buffer_values, buffer_values_id, funcConfigs)
            # Go to the next point
        end
    end
    return values_id, grid_id
end

function eval_at_point(agc, xi)
    xi_id = agc.xscale_inv(xi)
    yi_id = agc.func(xi_id)
    yi = apply_config_scales(agc.funcConfigs, yi_id)
    return xi, xi_id, yi, yi_id
end

function push_buffers!(xi, xi_id, yi, yi_id, buffer_x, buffer_x_id, buffer_values, buffer_values_id, funcConfigs)
    push!(buffer_x, xi) # Slightly unnessesary alloc
    push!(buffer_x_id, xi_id)
    map(save, buffer_values, yi, funcConfigs)
    map(save, buffer_values_id, yi_id, funcConfigs)
    return
end

function pop_buffers!(buffer_x, buffer_x_id, buffer_values, buffer_values_id, funcConfigs)
    xi = pop!(buffer_x)
    xi_id = pop!(buffer_x_id)
    yi = pop_tuple_of_vectors!(buffer_values, funcConfigs)
    yi_id = pop_tuple_of_vectors!(buffer_values_id, funcConfigs)
    return xi, xi_id, yi, yi_id
end

apply_scale(fc::AbstractFuncConfig, y_id) = fc.yscale(y_id)
#apply_scale(fc::SavingFuncConfig, y_id) = nothing

function apply_config_scales(funcConfigs, y_id)
    map(apply_scale, funcConfigs, y_id)
end

# TODO We can scale when saving instead of recomputing scaling
# Svectors should also be faster in general
function should_refine(v1, v2, v3, x1, x2, x3, agc::AbstractAutoGridConfig)
    not_linear = !all(is_almost_linear.(v1, v2, v3, x1, x2, x3, agc.funcConfigs, Ref(agc))) # Ref so struct is scalar
    if not_linear
        # println("Not linear")
        return not_linear
    else
        # println("Linear")
    end
    is_peak = any(potential_peak.(v1, v2, v3, agc.funcConfigs, Ref(agc)))
    # println("Peak: $is_peak")
    return is_peak
end

function potential_peak(y1, y2, y3, fc::AbstractFuncConfig, agc::AbstractAutoGridConfig)
    return fc.find_peaks && dotdiff(y2,y1,y3,y2) < 0
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

@inline function is_almost_linear(y1, y2, y3, x1, x2, x3, fc::AbstractFuncConfig, agc::AbstractAutoGridConfig)
    # We assume that x2-x1 \approx  x3-x2, so we need to check that y2-y1 approx y3-y2
    # Essentially low second derivative compared to derivative, so that linear approximation is good.
    # Second argument to avoid too small steps when derivatives are small
    dd = derivderiv2(y1,y2,y3,x1,x2,x3)
    #rhs = max((agc.reltol_dd)^2*max(deriv2(y2,y1,x2,x1), deriv2(y3,y2,x3,x2)), (agc.abstol_dd)^2)
    d = max(deriv2(y2,y1,x2,x1), deriv2(y3,y2,x3,x2))
    # println("$dd < $d")
    # println("$dd < $(max(d, agc.abstol_dd))")
    # println("$y1, $y2, $y3")
    # println("$(lhs < rhs) lhs/rhs: $lhs, $rhs")
    return dd <  agc.reltol_dd*max(d, agc.abstol_dd)/(x3-x1)^2
end
@inline is_almost_linear(y1, ym, y2, x1, x2, x3, fc::SavingFuncConfig, agc::AbstractAutoGridConfig) = true

deriv2(y1,y2,x1,x2) = normdiff2(y1,y2)/abs2(x1-x2)
# This will probably never (almost) be called
@inline normdiff2(y1::Number,y2::Number) = abs2(y1-y2)

function normdiff2(y1::A,y2::A) where {T<:Number, N, A <: Union{NTuple{N,T},AbstractArray{T}}}
    length(y1) == length(y2) || throw(DimensionMismatch())
    sum = abs2(y1[1]-y2[1])
    @inbounds for i in 2:length(y1)
        sum += abs2(y1[i]-y2[i])
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


function derivderiv2(y1,y2,y3,x1,x2,x3)
    # println("x:",x1,x2,x3)
    # println("y:",y1,y2,y3)
    d1 = abs(x1-x2)
    d2 = abs(x2-x3)
    diff = normdiffdiff2(y1,y2,y3,d1,d2)
    dsquare = (d1*d2*(d1+d2))^2/4
    # println("before: $diff, dsquare: $dsquare")
    diff/dsquare
end

# This will probably never(almost) be called
@inline normdiffdiff2(x::Number,y::Number,z::Number,d1,d2) = abs2(d2*x+d1*z-(d1+d2)*y)

function normdiffdiff2(x::A,y::A,z::A,d1,d2) where {T<:Number, N, A <: Union{NTuple{N,T},AbstractArray{T}}}
    length(x) == length(y) || throw(DimensionMismatch())
    sum = abs2(d2*x[1]+d1*z[1]-(d1+d2)*y[1])
    @inbounds for i in 2:length(x)
        sum += abs2(d2*x[i]+d1*z[i]-(d1+d2)*y[i])
    end
    return sum
end
function normdiffdiff2(x,y,z,d1,d2)
    length(x) == length(y) || throw(DimensionMismatch())
    sum = normdiffdiff2(x[1],y[1],z[1],d1,d2)
    @inbounds for i in 2:length(x)
        sum += normdiffdiff2(x[i],y[i],z[i],d1,d2)
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

@generated function tuple_of_vectors(yi_id::Tuple{Vararg{<:Any,N}}, funcConfigs::Tuple{Vararg{<:AbstractFuncConfig,N}}) where {N}
    vec = []
    for i in 1:N
        # println("Generating")
        # println(funcConfigs)
        # if funcConfigs.types[i] <: AbstractFuncConfig{true}
            push!(vec, :([yi_id[$i],]))
        # else
            # push!(vec, :(nothing))
        # end
    end
    :(Core.tuple($(vec...),))
end

@generated function pop_tuple_of_vectors!(vals::Tuple{Vararg{<:Any,N}}, funcConfigs::Tuple{Vararg{<:AbstractFuncConfig,N}}) where {N}
    vec = []
    for i in 1:N
        # println("Generating")
        # println(funcConfigs)
        # if funcConfigs.types[i] <: AbstractFuncConfig{true}
            push!(vec, :(pop!(vals[$i])))
        # else
            # push!(vec, :(nothing))
        # end
    end
    :(Core.tuple($(vec...),))
end

#save(values, ym_id, funcConfig::AbstractFuncConfig{false}) = nothing
function save(values, ym_id, funcConfig::AbstractFuncConfig)
    push!(values, ym_id)
    return nothing
end

end
