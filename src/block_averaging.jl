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
Base.length(ba::BlockAveragable) = length(first(ba.data))
Base.getindex(ba::BlockAveragable, slice::UnitRange{<:Integer}) = @views BlockAveragable([d[slice] for d in ba.data]...)
Base.getindex(ba::BlockAveragable, i::Integer) = ba.data[i]

struct BlockAveragedData{R,S,T}
    block_sizes::AbstractVector{Int}
    N_blocks::AbstractVector{Int}
    σ_estimates::AbstractVector{R}
    SE_estimates::AbstractVector{S}
    statistic_estimate::T
end

function single_block_average(ba::BlockAveragable, n_blocks::Integer, statistic::Function)

    N = length(ba)
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

    value = statistic(ba.data...)
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


function block_average_plot(bad::BlockAveragedData, outpath::String)
    #TODO
    #* x-axis log2
    # use plotlylight?
end


function do_block_averaging(X::BlockAveragable, statistic::Function, outpath::Union{String, Nothing})
    bad = block_average(X, statistic)
    converged_idx = select_plateau(bad)
    !isnothing(outpath) && block_average_plot(bad, outpath)
    return measurement(bad.statistic_estimate, bad.SE_estimates[converged_idx])
end


"""
@ba y = mean(rand(3,3))
@ba y = test(z, 2)
y = @ba test(z,2)
@ba nargs=2 y=mean([1,2,3],[2,3,4])
@ba nargs=2 path="/home" y=mean([1,2,3],[2,3,4])

Applies block averaging to f(data..., args...) with the function `f` specificying the statistic.
The first `nargs` arguments are the data to do the blocking analysis one. All
other args are captured and passed to the end of the function. Key word arguments 
specified by name will not work!

Optional Arguments:
- nargs: How many paramters to take from the front as data to sample from
- path: If the path argument is set, then a convergence plot will be created and saved. 
"""
macro ba(expressions...)

    arg_expressions = expressions[1:end-1]
    ex = expressions[end]
    path = nothing
    N::Int = 1

    # Parse arguments if they are passed
    for ae in arg_expressions
        if ae isa Expr && ae.head === :(=)
            arg_name = ae.args[1]
            if arg_name === :path
                path = ae.args[2]
            elseif arg_name === :nargs
                N = ae.args[2]
                @assert N isa Integer
            else
                Base.throw(ArgumentError("Unxpected arg in @ba. Expects nargs or path, got $(arg_name)"))
            end
        else
            Base.throw(ArgumentError("If passing param to @ba must be nargs=<int> and/or path=<str>"))
        end
    end

    result_variable::Union{Nothing, Symbol} = nothing
    fn_sig::Symbol = Symbol()
    fn_name::Symbol = Symbol()

    # Figure out if macro put before of after the = sign
    if ex.head === :(=)
        # @ba y = f(data, args...)
        result_variable = ex.args[1]
        fn_sig = ex.args[2]
        fn_name = fn_sig.args[1]
        data_end_idx = N + 1
        data = fn_sig.args[2:data_end_idx]
        kwargs = fn_sig.args[data_end_idx+1:end]
    elseif ex.head === (:call)
        # @ba f(data, args...)
        fn_name = ex.args[1]
        data_end_idx = N + 1
        data = ex.args[2:data_end_idx]
        kwargs = ex.args[data_end_idx+1:end]
    else
        Base.throw("Block average macro can only be applied to function calls or assignments")
    end

    # Generate code
    tmp_fn = gensym()
    tmp_ba = gensym()

    body = if isnothing(result_variable)
        quote
            local $tmp_fn = (d) -> $(fn_name)(d, $(kwargs...))
            local $tmp_ba = $CumulantAnalysis.BlockAveragable($(data...))
            $CumulantAnalysis.do_block_averaging($tmp_ba, $tmp_fn, $path)
        end
    else
        quote
            local $tmp_fn = (d) -> $(fn_name)(d, $(kwargs...))
            local $tmp_ba = $CumulantAnalysis.BlockAveragable($(data...))
            $(result_variable) = $CumulantAnalysis.do_block_averaging($tmp_ba, $tmp_fn, $path)
        end
    end

    return esc(body)

end