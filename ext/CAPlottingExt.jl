module CAPlottingExt

using CumulantAnalysis
using Plots
using HDF5

function CumulantAnalysis.block_average_plot(hdf5_path, outpath::String)

    file = h5open(hdf5_path, "r")
    bs = read(file["block_sizes"])
    sigs = read(file["all_std_estimates"])
    converged_idx = read(file["converged_idx"])
    original_data_length = read(file["original_data_length"])
    close(file)

    xticks = [bs[i] for i in 1:3:length(bs)]

    # start your scatter
    p = scatter(
        bs, sigs,
        xscale     = :log2,
        title      = "Total Elements: $(original_data_length)",
        xticks     = (xticks, string.(xticks)),
        xlabel     = "Elements per Block",
        ylabel     = "Standard Deviation",
        label      = "Std Dev",
        color      = "#4682b4",
        legend     = :topright,
        markersize = 6,
        show       = false
    );

    # add horizontal line at the converged σ
    hline!(
        p,
        [sigs[converged_idx]],
        linestyle = :dash,
        linewidth = 2,
        color     = :red,
        label     = "Converged σ",
        show      = false
    );

    savefig(p, outpath)

    return nothing
end




end # CAPlottingExt