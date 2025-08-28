using CairoMakie
using DelimitedFiles

material = "SW_SILICON"
get_datapath = (T,s) -> "/mnt/merged/emeitz/CumulantAnalysisTest/sTDEP_SW_TEST/T$(ustrip(T))/$(s)UC"
temperatures = [100, 1300]
F_true = [-4.2957902, -4.6756675]
F_tol = 1e-3 # eV / atom
sizes = [2,3,4,5,6,7,8]
order = 3
order_colors =  ["#080808", "#4287f5", "#f7952d", "#259c2f"]

F0s = zeros(length(temperatures), length(sizes))
F_terms = zeros(length(temperatures), length(sizes), order)
F_terms_SE = zeros(length(temperatures), length(sizes), order)
F_total = zeros(length(temperatures), length(sizes))
F_total_SE = zeros(length(temperatures), length(sizes))

for (i,T) in enumerate(temperatures)
    for (j,s) in sizes
        F_path = joinpath(get_datapath(T,s), "F_mean.txt")
        data = readdlm(F_path, skipstart = 1)
        F0s[i,j] = data[1] 
        F_terms[i,j,:] = data[2:2:end-3]
        F_terms_SE[i,j,:] = data[3:2:end-2]
        F_total[i,j] = data[end-1]
        F_total_SE[i,j] = data[end]
    end
end

# breaks down property by term in the cumulant expansion
function make_bar_plot(F0_T, F_terms_T, F_total_SE_T, order_colors)

    @assert length(order_colors) == order + 1 "Expected $(order + 1) colors, got $(length(order_colors))"

    f = Figure(resolution = size_in_pixels);
    ax = Axis(f[1,1], xlabel = L"Number of Conventional Cells (N x N x N)", ylabel = L"\mathcal{F} [eV / atom]",
        ylabelsize = 40, xlabelsize = 40, yticklabelsize = 30, xticklabelsize = 30,
        xticks = sizes, xgridvisible = false, ygridvisible = false, xticksmirrored = true,
        yticksmirrored = true, xticklabelpad = 4, xtickalign=1, ytickalign = 1)

    # xlims!(0,1350)
    # ylims!(0.4,0.75)

    xs = [fill(j, order + 1) for j in 1:length(sizes)]
    xs = reduce(vcat, xs)

    ys = [[F0_T[j], F_terms_T[j, :]...] for j in 1:length(sizes)]
    ys = reduce(vcat, ys)
    yerr = F_total_SE_T

    stack = [collect(1:length(order_colors)) for j in 1:length(sizes)]
    stack = reduce(vcat, stack)

    colors = [order_colors for _ in 1:length(sizes)]
    colors = reduce(vcat, colors)

    barplot(
        xs,
        ys,
        stack = stack,
        color = colors,
        axis = (xticks = (1:length(sizes), ["$s" for s in sizes]))
    )

    legend_labels = ["F$(o)" for o in 0:order]
    axislegend(ax, legend_labels, position = :cb, labelsize = 35, orientation = :vertical, 
                framevisible = false, nbanks = 4, labelhalign = :center, 
                colgap = 25, patchlabelgap = 12)

    save(joinpath(base_outpath,"SIZE_EFFECTS_$(material)_T$(T)_ORDER$(order)_COMPS.svg"), f)
    save(joinpath(base_outpath,"SIZE_EFFECTS_$(material)_T$(T)_ORDER$(order)_COMPS.png"), f)
end

# compares bulk to ground truth
function make_line_plot(sizes, F_total, F_true)
    f = Figure(resolution = size_in_pixels);
    ax = Axis(f[1,1], xlabel = L"Number of Conventional Cells (N x N x N)", ylabel = L"\mathcal{F} [eV / atom]",
        ylabelsize = 40, xlabelsize = 40, yticklabelsize = 30, xticklabelsize = 30,
        xticks = sizes, xgridvisible = false, ygridvisible = false, xticksmirrored = true,
        yticksmirrored = true, xticklabelpad = 4, xtickalign=1, ytickalign = 1)

    xmin = minimum(sizes) - 1
    xmax = maximum(sizes) + 1
    xlims!(xmin, xmax)

    h = hlines!(F_true, xmin, xmax, color = "#b31041", linestyle = :dash, linewidth = 4);
    band!(xmin:0.1:xmax, F_true - F_tol, F_true + F_tol, color= "#b31041", alpha = 0.4)

    s1 = scatter!(sizes, F_total, markersize = 30, color = :black);

    axislegend(ax, [h, s1], ["Ground Truth", "Cumulant Expansion, Order $(order)"],
                position = :rb, labelsize = 35, orientation = :vertical, framevisible = false, nbanks = 4, 
                labelhalign = :center, colgap = 25, patchlabelgap = 12)

    save(joinpath(base_outpath,"SIZE_EFFECTS_$(material)_T$(T)_ORDER$(order).svg"), f)
    save(joinpath(base_outpath,"SIZE_EFFECTS_$(material)_T$(T)_ORDER$(order).png"), f)
end


for (i,T) in enumerate(temperatures)
    @views make_line_plot(sizes, F_total[i,:], F_true[i])
    @views make_bar_plot(F0s[i,:], F_terms[i, :, :], F_total_SE[i], order_colors)
end
