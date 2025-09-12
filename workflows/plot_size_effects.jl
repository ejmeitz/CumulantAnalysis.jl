using CairoMakie
using DelimitedFiles
using HDF5

order = 2

get_outpath = (M,T) -> joinpath(get_basepath(M), "T$(T)")
get_datapath = (M,T,s) -> joinpath(get_outpath(M, T), "(s)UC")

unit_map = Dict(
    "F" => "eV / atom",
    "S" => "kB / atom",
    "U" => "eV / atom",
    "Cᵥ" => "kB / atom",
)

method_colors = ["#4287f5", "#f59b3b", "#fa144d"]
true_color = "#b31041"


# material = "SW_SILICON"
# get_basepath = (M) -> "/mnt/merged/emeitz/CumulantAnalysisTest/$(M)_LJ_size_effects"
# temperatures = [100, 1300]
# F_true = [-4.2957902, -4.6756675]
# F_tol = 1e-3 # eV / atom
# sizes = [2,3,4,5,6]

# SW:
# <V> @ 100K  = -933.9152951740965
# <V> @ 1300K = -898.2383308713825

# LJ
# <V> @ 10K = -451.404009
# <V> @ 80K = -401.9778595676087


# ground_truths = Dict(
#     "F" => Dict(100 => -4.2957902, 1300 => -4.6756675), #eV / atom
#     "F_tol" => 1e-3, #eV / atom
#     "U" => Dict(100 => -4.31081576, 1300 => -3.9912508), #eV / atom
#     "U_tol" => missing, #eV / atom
#     "Cv" => Dict(100 => 3.0162458, 1300 => 3.2265268), # eV / kB * atom
#     "Cv_tol" => missing, # 1% of DP -- eV / kB * atom
#     "S" => Dict(100 => -1.7436438, 1300 => 6.10948172), # just from S = (U - F) / T
#     "S_tol" => missing, # eV / kB * atom
# )

material = "LJ_ARGON"
get_basepath = (M) -> "/mnt/merged/emeitz/CumulantAnalysisTest/$(M)_LJ_size_effects"
temperatures = [10, 80]
sizes = [3,4,5,6,7]


ground_truths = Dict(
    "F" => Dict(10 => -0.0729854, 80 => -0.0820714), #eV / atom
    "F_tol" => 1e-3, #eV / atom
    "U" => Dict(10 => -0.07517623, 80 => -0.057791034), #eV / atom
    "U_tol" => missing, #eV / atom
    "Cv" => Dict(10 => 2.9662821, 80 => 2.829351244), # eV / kB * atom
    "Cv_tol" => 0.0, # eV / kB * atom
    "S" => Dict(10 => -0.11024685, 80 => 0.15272945), # just from S = (U - F) / T
    "S_tol" => missing, # eV / kB * atom
)


function parse_totals_per_size(prop, M, T)

    total = zeros(length(sizes))
    total_SE = zeros(length(sizes))

    for (j,s) in enumerate(sizes)
        path = joinpath(get_datapath(M,T,s), "$(prop)_mean.txt")
        file = h5open(path, "r")

        total[j] = read(file, "$(prop)_total")
        total_SE[j] = read(file, "$(prop)_total_SE")

        close(file)
    end

    return total, total_SE

end

function parse_nsamples_study_per_size(M, T)

    N = nothing
    k1s = Float64[]
    k1_SEs = Float64[]
    k2s = Float64[]
    k2_SEs = Float64[]
    k3s = Float64[]
    k3_SEs = Float64[]


    for s in sizes
        data = readdlm(joinpath(get_datapath(M, T, s), "outfile.size_study"))
        N = Int.(data[:,1])
        push!(k1s, data[:,2])
        push!(k1_SEs, data[:,3])
        push!(k2s, data[:,4])
        push!(k2_SEs, data[:,5])
        push!(k3s, data[:,6])
        push!(k3_SEs, data[:,7])
    end
    return N, k1s, k1_SEs, k2s, k2_SEs, k3s, k3_SEs
end

# compares bulk to ground truth
function make_line_plot(sizes, total, total_SE, truth, tol, T, 
                            prop, mthd, study_type::String, xlabel)
    size_in_inches = (3, 2.25)
    dpi = 300
    size_in_pixels = size_in_inches .* dpi

    f = Figure(resolution = size_in_pixels);
    ax = Axis(f[1,1], xlabel = xlabel, ylabel = "$(prop) [$(unit_map[prop])]",
        ylabelsize = 40, xlabelsize = 40, yticklabelsize = 30, xticklabelsize = 30,
        xticks = sizes, xgridvisible = false, ygridvisible = false, xticksmirrored = true,
        yticksmirrored = true, xticklabelpad = 4, xtickalign=1, ytickalign = 1)

    xmin = minimum(sizes) - 1
    xmax = maximum(sizes) + 1
    xlims!(xmin, xmax)
    ylims!(minimum([total; truth]) - 0.01, maximum([total; truth]) + 0.01)

    h = hlines!(truth, xmin, xmax, color = "#b31041", linestyle = :dash, linewidth = 4);
    band!(xmin:0.1:xmax, truth - tol, truth + tol, color= "#b31041", alpha = 0.4);
    s1 = errorbars!(sizes, total, total_SE, total_SE, markersize = 30, color = :black, whiskerwidth = 10);

    axislegend(ax, [[h], [s1]], [["Ground Truth"], ["Cumulant Expansion, Order $(order)"]],
                position = :lt, labelsize = 35, orientation = :horizontal, framevisible = false, nbanks = 2, 
                labelhalign = :left, colgap = 25, patchlabelgap = 12)

    save(joinpath(get_outpath(mthd,T),"$(mthd)_$(prop)_$(study_type)_$(material)_T$(T)_ORDER$(order).svg"), f)
    save(joinpath(get_outpath(mthd,T),"$(mthd)_$(prop)_$(study_type)_$(material)_T$(T)_ORDER$(order).png"), f)
end

function make_nsamples_plot(N, total, total_SE, T, 
                            prop, unit_str, mthd, study_type::String, xlabel)
    size_in_inches = (3, 2.25)
    dpi = 300
    size_in_pixels = size_in_inches .* dpi

    f = Figure(resolution = size_in_pixels);
    ax = Axis(f[1,1], xlabel = "Number of Samples", ylabel = "$(prop) [$(unit_str)]",
        ylabelsize = 40, xlabelsize = 40, yticklabelsize = 30, xticklabelsize = 30,
        xticks = N, xgridvisible = false, ygridvisible = false, xticksmirrored = true,
        yticksmirrored = true, xticklabelpad = 4, xtickalign=1, ytickalign = 1)

    s1 = errorbars!(N, total, total_SE, total_SE, markersize = 30, color = :black, whiskerwidth = 10);

    axislegend(ax, [s1], ["Cumulant Expansion, Order $(order)"],
                position = :lt, labelsize = 35, orientation = :horizontal, framevisible = false, nbanks = 2, 
                labelhalign = :left, colgap = 25, patchlabelgap = 12)

    save(joinpath(get_outpath(mthd,T),"$(mthd)_$(prop)_$(study_type)_$(material)_T$(T)_ORDER$(order).svg"), f)
    save(joinpath(get_outpath(mthd,T),"$(mthd)_$(prop)_$(study_type)_$(material)_T$(T)_ORDER$(order).png"), f)
end

function make_method_comp_barplot(vals, val_SE, truth, tol, colors, T, prop, mthd_labels::AbstractVector{String})
    size_in_inches = (3, 2.25)
    dpi = 300
    size_in_pixels = size_in_inches .* dpi

    f = Figure(resolution = size_in_pixels);
    ax = Axis(f[1,1], xticks = (1:3, mthd_labels), ylabel = "$(prop) [$(unit_map[prop])]",
        ylabelsize = 40, xlabelsize = 40, yticklabelsize = 30, xticklabelsize = 30,
        xgridvisible = false, ygridvisible = false, xticksmirrored = true,
        yticksmirrored = true, xticklabelpad = 4, xtickalign=1, ytickalign = 1)

    xlims!(0, 4)

    barplot!(1:3, vals, color = colors)
    errorbars!(1:3, vals, val_SE, val_SE; color = :black, linewidth = 4, overdraw = true, whiskerwidth = 8)
    hlines!(truth, 0, 4, color = true_color, linestyle = :dash, linewidth = 4);
    band!(0:0.1:4, truth - tol, truth + tol, color= true_color, alpha = 0.4);

    eh = PolyElement(color = true_color, alpha = 0.4)

    axislegend(ax, [eh], ["Ground Truth"], position = :rb, labelsize = 35, orientation = :vertical, 
                framevisible = false, nbanks = 4, labelhalign = :center, 
                colgap = 25, patchlabelgap = 12)

    save(joinpath(get_outpath(T),"MTHD_COMP_$(prop)_$(material)_T$(T)_ORDER$(order).svg"), f)
    save(joinpath(get_outpath(T),"MTHD_COMP_$(prop)_$(material)_T$(T)_ORDER$(order).png"), f)

end

for mthd in ("MDV0", "HarmV0", "QuarticV0")
    for (i,T) in enumerate(temperatures)
        for prop in ("F", "S", "Cv", "U")
            tol = ground_truths["$(prop)_tol"]
            truth = ground_truths["$(prop)"][T]
            
            all_totals = []; all_total_SEs = []; 
            all_k1 = []; all_k1_SEs = [];
            all_k2 = []; all_k2_SEs = [];
            N = nothing

            mthd_labels = ["MDV0", "HarmV0", "QuarticV0"]
            for mthd in mthd_labels
                total, total_SE = parse_totals_per_size(prop, order, T)
                push!(all_totals, total); push!(all_total_SEs, total_SE)
                make_line_plot(sizes, total, total_SE, truth, tol, T, prop, mthd,
                                "SIZE_EFF", "Number of Conventional Cells (N x N x N)")

                N, k1s, k1_SEs, k2s, k2_SEs, k3s, k3_SEs = parse_nsamples_study_per_size(mthd, T)
                make_nsamples_plot(N, k1s, k1_SEs, T, "β * κ₁", "dimensionless", mthd, "")

                #! PLOT of F,S,U,Cv vs N not just cumulants
            end

            make_method_comp_barplot(all_totals, all_total_SEs, truth, tol, method_colors, T, prop, mthd_labels)
        end
    end
end





# breaks down property by term in the cumulant expansion
# function make_bar_plot(F0_T, F_terms_T, F_total_T, F_total_SE_T, F_true, order_colors, T)

#     @assert length(order_colors) == order + 1 "Expected $(order + 1) colors, got $(length(order_colors))"

#     size_in_inches = (3, 2.25)
#     dpi = 300
#     size_in_pixels = size_in_inches .* dpi

#     f = Figure(resolution = size_in_pixels);
#     ax = Axis(f[1,1], xlabel = L"Number of Conventional Cells (N x N x N)", ylabel = L"\mathcal{F} [eV / atom]",
#         ylabelsize = 40, xlabelsize = 40, yticklabelsize = 30, xticklabelsize = 30,
#         xticks = sizes, xgridvisible = false, ygridvisible = false, xticksmirrored = true,
#         yticksmirrored = true, xticklabelpad = 4, xtickalign=1, ytickalign = 1)

#     xmax = maximum(sizes) + 1
#     xlims!(0, xmax)

#     xs = [fill(j, order + 1) for j in 1:length(sizes)]
#     xs = reduce(vcat, xs)

#     ys = [[F0_T[j], F_terms_T[j, :]...] for j in 1:length(sizes)]
#     ys = reduce(vcat, ys)
#     yerr = F_total_SE_T

#     stack = [collect(1:length(order_colors)) for j in 1:length(sizes)]
#     stack = reduce(vcat, stack)

#     colors = [order_colors for _ in 1:length(sizes)]
#     colors = reduce(vcat, colors)

#     band!(0:0.1:xmax, F_true - F_tol, F_true + F_tol, color= "#b31041", alpha = 0.4)

#     barplot!(
#         xs,
#         ys,
#         stack = stack,
#         color = colors,
#     )

#     # Outline each bar with black
#     # hopefully makes it clear the 
#     # errorbar is for the whole thing
#     barplot!(
#         collect(1:length(sizes)),
#         F_total_T,
#         color = RGBAf(0,0,0,0),
#         strokewidth = 4,
#         strokecolor = :black,
#     )

#     # Error bar for total error on F
#     errorbars!(sizes, F_total_T, yerr, yerr; color = :black, linewidth = 4, overdraw = true, whiskerwidth = 8)

#     eh = PolyElement(color = "#b31041", alpha = 0.4)
#     ef = [PolyElement(color = oc) for oc in order_colors]
#     legend_markers = [eh; ef]

#     legend_labels = ["Ground Truth"; ["F$(o)" for o in 0:order]]
#     axislegend(ax, legend_markers, legend_labels, position = :cb, labelsize = 35, orientation = :vertical, 
#                 framevisible = false, nbanks = 4, labelhalign = :center, 
#                 colgap = 25, patchlabelgap = 12)

#     save(joinpath(get_outpath(T),"SIZE_EFFECTS_$(material)_T$(T)_ORDER$(order)_COMPS.svg"), f)
#     save(joinpath(get_outpath(T),"SIZE_EFFECTS_$(material)_T$(T)_ORDER$(order)_COMPS.png"), f)
# end
