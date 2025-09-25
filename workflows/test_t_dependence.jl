using CairoMakie
using Statistics
using DelimitedFiles

basepath = "/mnt/merged/emeitz/CumulantAnalysisTest/LJ_TEST_MDV0"
true_energy_path = (T) -> joinpath(basepath, "T$(T)", "6UC", "outfile.energies")
get_path = (T) -> joinpath(basepath, "T$(T)", "6UC", "outfile.energies")

temps = [10, 20, 30, 40, 50, 60, 70, 80]

function parse_V2(path)

    data = readdlm(path, comments = true)

    # energies in file are meV / atom, Converts to eV / atom
    conv = 1 / 1000
    
    @views V2 = data[:, 4] .* conv
    @views V3 = data[:, 5] .* conv
    @views V4 = data[:, 6] .* conv

    return V2, V3, V4
end

V2_means = []
V4_means = []

for T in temps
    p = get_path(T)
    V = vec(readdlm())
    V2, V3, V4 = parse_V2(p)
    push!(V2_means, mean(V2))
    push!(V4_means, mean(V4))
end

kB = 8.617333262145e-5 # eV/K
V2_means ./= kB
V4_means ./=  (kB .* temps .^ 2)

###########

size_in_inches = (3, 2.25)
dpi = 300
size_in_pixels = size_in_inches .* dpi

f = Figure(size = size_in_pixels);
ax = Axis(f[1,1], xlabel = "T [K]", ylabel = "<V2>/(N*kB) (eV/kB/atom)",
    ylabelsize = 40, xticks = [20, 40, 60, 80], xlabelsize = 40, yticklabelsize = 30, xticklabelsize = 30,
    xgridvisible = false, ygridvisible = false, xticksmirrored = true,
    yticksmirrored = true, xticklabelpad = 4, xtickalign=1, ytickalign = 1,
    yautolimitmargin = (0.15f0, 0.15f0))

s = scatter!(temps, V2_means, markersize = 30)

xs = collect(0:0.1:85)
ys = 3 .* (0.5 .* xs)

lines!(xs, ys, color = :red, linewidth = 4, linestyle = :dash)

xlims!(0,85)
ylims!(0, maximum(V2_means) * 1.1)

save(joinpath(basepath,"V2_mean_vs_T.png"), f)

######

f = Figure(size = size_in_pixels);
ax = Axis(f[1,1], xlabel = "T [K]", ylabel = "<V4>/(N*kB*T^2) (eV/kB/atom)",
    ylabelsize = 40, xticks = [20, 40, 60, 80], xlabelsize = 40, yticklabelsize = 30, xticklabelsize = 30,
    xgridvisible = false, ygridvisible = false, xticksmirrored = true,
    yticksmirrored = true, xticklabelpad = 4, xtickalign=1, ytickalign = 1,
    yautolimitmargin = (0.15f0, 0.15f0))

s = scatter!(temps, V4_means, markersize = 30)

# xs = collect(0:0.1:85)
# ys = 3 .* (2 ./ xs)

# lines!(xs, ys, color = :red, linewidth = 4, linestyle = :dash)

xlims!(0,85)
ylims!(0, maximum(V4_means) * 1.1)

save(joinpath(basepath,"V4_mean_vs_T.png"), f)