export  CumulantCorrections, save_errors

struct CumulantCorrections{O,T,G}
    harmonic::T
    corrections::SVector{O,T} # O is order of corrections
    ground_truth::G
    property::String
    unit_str::String
end

function CumulantCorrections(harmonic, corrections, ground_truth, property, unit_str)

    return CumulantCorrections{length(corrections), typeof(harmonic), typeof(ground_truth)}(
                                harmonic, corrections, ground_truth, property, unit_str)
end

estimate(cc::CumulantCorrections) = cc.harmonic + sum(cc.corrections)

function Base.getindex(cc::CumulantCorrections{O}, i::Integer) where O
    if i > O
        raise(ArgumentError("Cannot get correction $(i), corrections only calculated to order $(O)"))
    end
    
    if i == 0
        return cc.harmonic
    else
        return cc.corrections[i]
    end
end

CumulantAnalysis.mean(vcc::Vector{<:CumulantCorrections}) = Statistics.mean(estimate.(vcc))
CumulantAnalysis.mean(vcc::Vector{<:CumulantCorrections}, order::Integer) = Statistics.mean(getindex.(vcc, order))

CumulantAnalysis.std(vcc::Vector{<:CumulantCorrections}) = Statistics.std(estimate.(vcc))
CumulantAnalysis.std(vcc::Vector{<:CumulantCorrections}, order::Integer) = Statistics.std(getindex.(vcc, order))

se(vcc::Vector{<:CumulantCorrections}) = CumulantAnalysis.std(vcc) / sqrt(length(vcc))
se(vcc::Vector{<:CumulantCorrections}, order::Integer) = CumulantAnalysis.std(vcc, order) / sqrt(length(vcc))

harmonic_part(cc::CumulantCorrections) = cc.harmonic
correction(cc::CumulantCorrections) = sum(cc.corrections)



# Calculates mean and se for each order
function save_errors(vcc::Vector{<:CumulantCorrections{ORDER}}, outdir::String) where ORDER
    prop_name = first(vcc).property
    unit_str = first(vcc).unit_str

    outpath_mean = (ext) -> joinpath(outdir, prop_name * "_mean.$(ext)")
    outpath_seed = (ext) -> joinpath(outdir, prop_name * "_seed.$(ext)")

    mean_data = OrderedDict(prop_name*"0" => CumulantAnalysis.mean(vcc, 0))
    seed_data = OrderedDict(prop_name*"0" => harmonic_part.(vcc))

    float_fmt_str = (N) -> Printf.Format(join(fill("%15.7f", N), " "))
    str_fmt_str = (N) -> Printf.Format(join(fill("%15s", N), " "))

    for order in 1:ORDER
        mean_data[prop_name * "$(order) $(unit_str)"] = CumulantAnalysis.mean(vcc, order)
        mean_data[prop_name * "$(order)_SE"] = se(vcc, order)

        seed_data[prop_name * "$(order) $(unit_str)"] = getindex.(vcc, order)
    end

    mean_data[prop_name*"_total $(unit_str)"] = CumulantAnalysis.mean(vcc)
    mean_data[prop_name*"_total_SE"] = se(vcc)


    # Human Readable Version
    mean_header = collect(keys(mean_data))
    mean_values = collect(values(mean_data))
    N = length(mean_header)
    open(outpath_mean("txt"), "w") do f
        println(f, Printf.format(str_fmt_str(N), mean_header...))
        println(f, Printf.format(float_fmt_str(N), mean_values...))
    end

    seed_header = collect(keys(seed_data))
    seed_values = collect(values(seed_data))
    N = length(seed_header)
    N_seeds = length(vcc)
    open(outpath_seed("txt"), "w") do f
        println(f, Printf.format(str_fmt_str(N), seed_header...))
        for i in 1:N_seeds
            println(f, Printf.format(float_fmt_str(N), getindex.(seed_values, i)...))
        end
    end

    # Better for computers...
    save(outpath_mean("h5"), mean_data)
    save(outpath_seed("h5"), seed_data)

end

function save_errors(cc::BootstrapCumualantEstimate{ORDER}, outdir::String) where ORDER
    prop_name = cc.property
    unit_str = cc.unit_str

    outpath_mean = (ext) -> joinpath(outdir, prop_name * "_mean.$(ext)")
    mean_data = OrderedDict(prop_name*"0" => cc.harmonic)

    for order in 1:ORDER
        mean_data[prop_name * "$(order) $(unit_str)"] = cc.corrections[order]
        mean_data[prop_name * "$(order)_SE"] = cc.correction_SEs[order]
    end

    mean_data[prop_name*"_total $(unit_str)"] = cc.total
    mean_data[prop_name*"_total_SE"] = cc.total_SE

    float_fmt_str = (N) -> Printf.Format(join(fill("%15.7f", N), " "))
    str_fmt_str = (N) -> Printf.Format(join(fill("%15s", N), " "))

     # Human Readable Version
    mean_header = collect(keys(mean_data))
    mean_values = collect(values(mean_data))
    N = length(mean_header)
    open(outpath_mean("txt"), "w") do f
        println(f, Printf.format(str_fmt_str(N), mean_header...))
        println(f, Printf.format(float_fmt_str(N), mean_values...))
    end
    # Save to HDF5
    save(outpath_mean, mean_data)
end