using CairoMakie
using DelimitedFiles
using HDF5

# material = "SW_SILICON"
# get_outpath = (T) -> "/mnt/merged/emeitz/CumulantAnalysisTest/sTDEP_SW_size_effects/T$(T)"
# get_datapath = (T,s) -> "/mnt/merged/emeitz/CumulantAnalysisTest/sTDEP_SW_size_effects/T$(T)/$(s)UC"
# temperatures = [100, 1300]
# F_true = [-4.2957902, -4.6756675]
# F_tol = 1e-3 # eV / atom
# sizes = [2,3,4,5,6,7]

# SW:
# <V> @ 100K  = -933.9152951740965
# <V> @ 1300K = -898.2383308713825

# LJ
# <V> @ 10K = -451.404009
# <V> @ 80K = -401.9778595676087


ground_truths = Dict(
    "F" => Dict(100 => -4.2957902, 1300 => -4.6756675), #eV / atom
    "F_tol" => 1e-3, #eV / atom
    "U" => Dict(100 => -4.31081576, 1300 => -3.9912508), #eV / atom
    "U_tol" => missing, #eV / atom
    "Cv" => Dict(100 => 3.0162458, 1300 => 3.2265268), # eV / kB * atom
    "Cv_tol" => missing, # 1% of DP -- eV / kB * atom
    "S" => Dict(100 => -1.7436438, 1300 => 6.10948172), # just from S = (U - F) / T
    "S_tol" => missing, # eV / kB * atom
)

material = "LJ_ARGON"
get_outpath = (T) -> "/mnt/merged/emeitz/CumulantAnalysisTest/sTDEP_LJ_size_effects/T$(T)"
get_datapath = (T,s) -> "/mnt/merged/emeitz/CumulantAnalysisTest/sTDEP_LJ_size_effects/T$(T)/$(s)UC"
temperatures = [10, 80]

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

sizes = [3,4,5,6,7,8]

order = 2 # order to actually use in plots
# max_order = 2
# order_colors =  ["#080808", "#4287f5", "#f7952d", "#259c2f"]

F0s = zeros(length(temperatures), length(sizes))
F_terms = zeros(length(temperatures), length(sizes), order)
F_terms_SE = zeros(length(temperatures), length(sizes), order)
F_total = zeros(length(temperatures), length(sizes))
F_total_SE = zeros(length(temperatures), length(sizes))


function parse_totals(prop, order)
    # harmonic = zeros(length(temperatures), length(sizes))
    # offsets = zeros(length(temperatures), length(sizes))
    # offsets_SE = zeros(length(temperatures), length(sizes))
    # terms = zeros(length(temperatures), length(sizes), order)
    # terms_SE = zeros(length(temperatures), length(sizes), order)
    total = zeros(length(temperatures), length(sizes))
    total_SE = zeros(length(temperatures), length(sizes))

    for (i,T) in enumerate(temperatures)
        for (j,s) in enumerate(sizes)
            path = joinpath(get_datapath(T,s), "F_mean.txt")
            file = h5open(path, "r")

            # harmonic[i,j] = read(file, "$(prop)0")

            # offsets[i,j] = read(file, "$(prop)_offset")
            # offsets_SE[i,j] = read(file, "$(prop)_offset_SE")

            total[i,j] = read(file, "$(prop)_total")
            total_SE[i,j] = read(file, "$(prop)_total_SE")

            close(file)
        end
    end

    return total, total_SE

end


# compares bulk to ground truth
function make_line_plot(sizes, total, total_SE, truth, tol, T, prop)
    size_in_inches = (3, 2.25)
    dpi = 300
    size_in_pixels = size_in_inches .* dpi

    f = Figure(resolution = size_in_pixels);
    ax = Axis(f[1,1], xlabel = "Number of Conventional Cells (N x N x N)", ylabel = "F [eV / atom]",
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

    save(joinpath(get_outpath(T),"$(prop)_SZ_EFF_$(material)_T$(T)_ORDER$(order).svg"), f)
    save(joinpath(get_outpath(T),"$(prop)_SZ_EFF_$(material)_T$(T)_ORDER$(order).png"), f)
end


for (i,T) in enumerate(temperatures)
    for prop in ("F", "S", "Cv", "U")
        total, total_SE = parse_totals(prop, order)
        tol = ground_truths["$(prop)_tol"]
        truth = ground_truths["$(prop)"][T]
        make_line_plot(sizes, total, total_SE, truth, tol, T, prop)
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
