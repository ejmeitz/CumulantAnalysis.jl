struct BlockAveragedData{R,S,T}
    block_sizes::AbstractVector{Int}
    N_blocks::AbstractVector{Int}
    σ_estimates::AbstractVector{R}
    SE_estimates::AbstractVector{S}
    statistic_estimate::T
end

function single_block_average(data::AbstractVector, n_blocks::Integer, statistic::Function)

    N = length(x)
    block_size = div(N, n_blocks)
    N_max = block_size * n_blocks # ignore data at end to keep block size consistent
    block_ranges = [s : min(s + block_size - 1, N) for s in 1:block_size:N_max]

    block_values = statistic.(view.(Ref(data), block_ranges))
    σ = std(block_values)
    SE_estimate = σ / sqrt(n_blocks)

    return block_size, σ, SE_estimate
end

# Calculates block averages from N_blocks = 2, to a dynamic
# block size chosen such that the smallest block contains
# `min_block_size` elements from the original array.
function block_average(data::AbstractVector, statistic::Function; min_block_size::Integer = 25)

    # Space n_blocks by powers of 2 such that smallest block
    # has at least min_block_size elements
    N = length(x)
    K = div(N, min_block_size)
    max_blocks = max(K, 2) # 2 here is minimum number of blocks
    max_exp = floor(Int, log2(max_blocks))
    n_blocks = 2 .^ (1:max_exp)

    value = statistic(data)
    results = single_block_average.(Ref(data), n_blocks, Ref(statistic))
    block_sizes = getindex.(results, 1)
    σ_estimates = getindex.(results, 2)
    SE_estimates = getindex.(results, 3)

    return BlockAveragedData(block_sizes, n_blocks, σ_estimates, SE_estimates, value)
end


# Returns the index where standard deviation is within `tol` percent
# of the the last `lookback` values averaged. If no values are found
# the last index is returned. 
function select_plateau(bad::BlockAveragedData; tol::Real=0.05,  lookback::Int=3)

    bs = bad.block_sizes
    σ = bad.σ_estimates

    n = length(σ)
    @assert length(bs) == n "Block sizes and σ estimates must match length"

    # make sure lookback ≤ n
    L = min(lookback, n)
    ref = mean(σ[n-L+1 : n])

    idx = findfirst(s -> s >= (1 - tol)*ref, σ)
    return isnothing(idx) ? lastindex(bs) : idx

end


function block_average_plot(bad::BlockAveragedData)
    #TODO
    #* x-axis log2
    # use plotlylight?
end


function do_block_averaging(X, statistic)
    bad = block_average(X, statistic)
    converged_idx = select_plateau(bad)
    #* MAKE block average plot??
    return measurement(bad.statistic_estimate, bad.SE_estimate[converged_idx])
end


"""
Applies block averaging to the data with the function specificying the statistic.
The first argument is assumed to be the data to do the blocking analysis one. All
other args are captured and passed to the end of the function. Key word arguments 
specified by name will not work!

For example,
````
@macroexpand @block_average y = mean(rand(3,3))
quote
    __ba_fn = (ba_data->mean(ba_data))
    y = do_block_averaging(rand(3, 3), __ba_fn)
end

@macroexpand y = @block_average mean(rand(3,3))
:(y = begin
    __ba_fn = (ba_data->mean(ba_data))
    do_block_averaging(rand(3, 3), __ba_fn)
end)

@macroexpand @block_average y = test(z, 2)
quote
    __ba_fn = (ba_data->test(ba_data, 2))
    y = do_block_averaging(z, __ba_fn)
end


```
"""
macro block_average(ex)

    result_variable = nothing

    # Figure out if macro put before of after the = sign
    if ex.head == :(=)
        # @block_average y = f(data; kwargs...)
        result_variable = ex.args[1]
        fn_sig = ex.args[2]
        fn_name = fn_sig.args[1]
        data = fn_sig.args[2]
        kwargs = fn_sig.args[3:end]
    elseif ex.head == (:call)
        # @block_average f(data; kwargs...)
        fn_name = ex.args[1]
        data = ex.args[2]
        kwargs = ex.args[3:end]
    else
        Base.throw("Block average macro can only be applied to function calls or assignments")
    end


    fn_sym = :__ba_fn

    # Build the "statistic" function with 1 arg. Kwargs are
    # captured since they exist in same scope
    #   __ba_fn = (ba_data) -> f(ba_data, kwargs...) 
    set_new_function = Expr(
        :(=), # head
        fn_sym, # lhs
        Expr(
            :->,
            :ba_data,
            Expr(:call, fn_name, :ba_data, kwargs...)
        ) # rhs
    )

    call_block_avg = Expr(:call, :do_block_averaging, data, fn_sym)

    body = if isnothing(result_variable)
        Expr(:block, set_new_function, call_block_avg)
    else
        Expr(:block,
             set_new_function,
             Expr(:(=), result_variable, call_block_avg)
        )
    end

    return esc(body)

end