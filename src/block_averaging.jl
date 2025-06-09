struct BlockAveragable{N,T}
    data::AbstractVector{T}
end

function BlockAveragable(args...)
    N = length(args)
    @assert N > 0 "Tried to create empty BlockAveragable??"
    @assert allequal(length(a) for a in args) "BlockAveragable data must have same length"
    return BlockAveragable{N, typeof(args[1])}([args...])
end

n_datasets(::BlockAveragable{N}) where N = N
length(ba::BlockAveragable) = length(first(ba.data))
getindex(ba::BlockAveragable, slice::UnitRange{<:Integer}) = @views BlockAveragable([d[slice] for d in ba.data]...)
getindex(ba::BlockAveragable, i::Integer) = ba.data[i]

struct BlockAveragedData{R,S,T}
    block_sizes::AbstractVector{Int}
    N_blocks::AbstractVector{Int}
    σ_estimates::AbstractVector{R}
    SE_estimates::AbstractVector{S}
    statistic_estimate::T
end

function single_block_average(ba::BlockAveragable, n_blocks::Integer, statistic::Function)

    N = length(x)
    block_size = div(N, n_blocks)
    N_max = block_size * n_blocks # ignore data at end to keep block size consistent
    block_ranges = [s : min(s + block_size - 1, N) for s in 1:block_size:N_max]

    # getindex for BlockAveragable returns a new
    # BlockAveragable which stores the slices as views
    ba_slices = getindex.(Ref(ba), block_ranges)
    block_values = [statistic(ba_slice.data...) for ba_slice in ba_slices]
    σ = std(block_values)
    SE_estimate = σ / sqrt(n_blocks)

    return block_size, σ, SE_estimate
end

# Calculates block averages from N_blocks = 2, to a dynamic
# block size chosen such that the smallest block contains
# `min_block_size` elements from the original array.
function block_average(ba::BlockAveragable, statistic::Function; min_block_size::Integer = 25)

    # Space n_blocks by powers of 2 such that smallest block
    # has at least min_block_size elements
    N = length(ba)
    K = div(N, min_block_size)
    max_blocks = max(K, 2) # 2 here is minimum number of blocks
    max_exp = floor(Int, log2(max_blocks))
    n_blocks = 2 .^ (1:max_exp)

    value = statistic(ba.data)
    results = single_block_average.(Ref(ba), n_blocks, Ref(statistic))
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


function do_block_averaging(X::BlockAveragable, statistic::Function)
    bad = block_average(X, statistic)
    converged_idx = select_plateau(bad)
    #* MAKE block average plot??
    return measurement(bad.statistic_estimate, bad.SE_estimate[converged_idx])
end


"""
@ba y = mean(rand(3,3))
@ba y = test(z, 2)
y = @ba test(z,2)
@ba nargs=2 y=mean([1,2,3],[2,3,4])

Applies block averaging to the data with the function specificying the statistic.
The first argument is assumed to be the data to do the blocking analysis one. All
other args are captured and passed to the end of the function. Key word arguments 
specified by name will not work! Only functions like f(data..., args...) will work
where nargs optionally specifies how many args to expect in data.

For example,
```
@macroexpand @ba y = mean(rand(3,3))
quote
    __ba_fn = (ba_data->mean(ba_data))
    y = do_block_averaging(rand(3, 3), __ba_fn)
end

@macroexpand y = @ba mean(rand(3,3))
:(y = begin
    __ba_fn = (ba_data->mean(ba_data))
    do_block_averaging(rand(3, 3), __ba_fn)
end)

@macroexpand @ba y = test(z, 2)
quote
    __ba_fn = (ba_data->test(ba_data, 2))
    y = do_block_averaging(z, __ba_fn)
end
```

If your function takes multiple inputs you can specify the number of 
arguments which should be sampled concurrently. For example,
```
@macroexpand @ba nargs=2 y = mean([1,2,3],[2,3,4])
quote
    __ba_fn = (ba_data->mean(ba_data))
    y = do_block_averaging([1, 2, 3], [2, 3, 4], __ba_fn)
end
```
"""
macro ba(nargs, ex)

    if nargs isa Expr && nargs.head === :(=)
        N = nargs.args[2]
        @assert N isa Integer
    else
        Base.throw(ArgumentError("If passing param to @ba must be nargs=<int>"))
    end

    result_variable = nothing

    # Figure out if macro put before of after the = sign
    if ex.head === :(=)
        # @ba y = f(data, args...)
        result_variable = ex.args[1]
        fn_sig = ex.args[2]
        fn_name = fn_sig.args[1]
        data_end_idx = 2 + N - 1
        data = fn_sig.args[2:data_end_idx]
        kwargs = fn_sig.args[data_end_idx+1:end]
    elseif ex.head === (:call)
        # @ba f(data, args...)
        fn_name = ex.args[1]
        data_end_idx = 2 + N - 1
        data = ex.args[2:data_end_idx]
        kwargs = ex.args[data_end_idx+1:end]
    else
        Base.throw("Block average macro can only be applied to function calls or assignments")
    end

    tmp_fn = gensym()
    tmp_ba = gensym()

    body = if isnothing(result_variable)
        quote
            local $tmp_fn = (d) -> $(fn_name)(d, $(kwargs...))
            local $tmp_ba = $CumulantAnalysis.BlockAveragable($(data...))
            $CumulantAnalysis.do_block_averaging($tmp_ba, $tmp_fn)
        end
    else
        quote
            local $tmp_fn = (d) -> $(fn_name)(d, $(kwargs...))
            local $tmp_ba = $CumulantAnalysis.BlockAveragable($(data...))
            $result_variable = $CumulantAnalysis.do_block_averaging($tmp_ba, $tmp_fn)
        end
    end

    return esc(body)

end

macro ba(ex)
    return esc(
        quote
            $CumulantAnalysis.@ba num=1 $ex
        end
    )
end