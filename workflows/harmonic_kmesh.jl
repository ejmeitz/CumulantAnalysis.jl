using LatticeDynamicsToolkit
import LatticeDynamicsToolkit: Hartree_to_eV, kB_eV
using CairoMakie
# using LaTeXStrings

limit = Classical
ks = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
outpath = raw"Z:\emeitz\Data\FreeEnergies\KMeshConvergence"

# LENNARD JONES
# pot = "LJ"
# Ts = [10, 80]
# ucposcar_path = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data\LJ\infile.ucposcar"
# ifc2_path = (T) -> "Z:/emeitz/Data/ForceConstants/LJ_IFC_INTERPOLATION_NODES/IFCs/T$(T)/infile.forceconstant"
# uc = CrystalStructure(ucposcar_path)

# STILLINGER WEBER
pot = "SW"
Ts = [100, 1300]
ucposcar_path = raw"C:\Users\ejmei\repos\TDEP_IFCs.jl\data\SW\infile.ucposcar"
ifc2_path = (T) -> "Z:/emeitz/Data/ForceConstants/SW_IFC_INTERPOLATION_NODES/IFCs/T$(T)/infile.forceconstant"
uc = CrystalStructure(ucposcar_path)


Fs = zeros(length(Ts), length(ks))
Ss = zeros(length(Ts), length(ks))
Us = zeros(length(Ts), length(ks))
Cvs = zeros(length(Ts), length(ks))


for (i,T) in enumerate(Ts)
    ifc2 = read_ifc2(ifc2_path(T), ucposcar_path)

    for (j,k) in enumerate(ks)
        @info "T = $(T), k = $(k)"
        Fs[i,j], Ss[i,j], Us[i,j], Cvs[i,j] = Hartree_to_eV .* harmonic_properties(
                    T, 
                    uc,
                    ifc2,
                    [k,k,k],
                    limit;
                )
    end
end

# conert to meV or normalize by kB
Fs  .*= 1000
Ss  ./= kB_eV
Us  .*= 1000
Cvs ./= kB_eV

ylabels = Dict(
    "F"  => L"\text{Absolute Error } [{\mathrm{meV}}/\mathrm{atom}]",
    "S"  => L"\text{Absolute Error } [k_{\mathrm{B}}]",
    "U"  => L"\text{Absolute Error } [{\mathrm{meV}}/\mathrm{atom}]",
    "Cv" => L"\text{Absolute Error } [k_{\mathrm{B}}]"
)

title_map = Dict(
    "F" => "$(pot) Free Energy",
    "S" => "$(pot) Entropy",
    "U" => "$(pot) Internal Energy",
    "Cv" => "$(pot) Heat Capacity"
)

function get_tol(prop, T)
    if prop == "F"
        return 1.0 # 1 meV / atom
    elseif prop == "S"
        return 1e-3/(kB_eV*T) # unitless # ignores error from U, but calculated form F = U - TS
    elseif prop == "U"
        return 1.0 # just copied F
    elseif prop == "Cv"
        return 0.03 # 1% of Dulong-Petit in units of kB / atom
    else
        error("Unrecognized property $(prop). Cannot get tolerance.")
    end
end

all_data = Dict(
    "F"  => Fs,
    "S"  => Ss,
    "U"  => Us,
    "Cv" => Cvs
)

size_in_inches = (3, 2.25)
dpi = 300
size_in_pixels = size_in_inches .* dpi

for (key, data) in all_data
    for (j, T) in enumerate(Ts)

        if limit == Classical && key ∈ ["U", "Cv"]
            continue
        end


        f = Figure(size = size_in_pixels);
        ax = Axis(f[1,1], xlabel = "1 / k", ylabel = ylabels[key],
            ylabelsize = 40, xlabelsize = 40, yticklabelsize = 30, xticklabelsize = 30,
            xscale = log10, yscale = log10,
            xgridvisible = false, ygridvisible = false, xticksmirrored = true,
            yticksmirrored = true, xticklabelpad = 4, xtickalign=1, ytickalign = 1,
            title = "$(title_map[key]) T = $(T)K", titlesize = 40
        )

        X   = 1.0 ./ ks
        Y   = vec(data[j,:])
        fit = fit_Qinf_vs_k!(ks, Y)  # your fitting function

        p = fit.p_best

        idxs = [1,2,3,4,5,7,9]
        tick_labels = [L"\frac{1}{%$(k)}" for k in ks] 
        tick_labels = tick_labels[idxs]
        ax.xticks = (X[idxs], tick_labels)
        ax.xlabel = L"1/k"

        err      = abs.(Y .- fit.Qinf)
        err_plot = max.(err, eps())   # avoid zero on log scale

        # Figure out y-ticks 
        ymin = minimum(err_plot)
        ymax = maximum(err_plot)
        set_log10_integer_power_ticks!(ax; ymin=ymin, ymax=ymax, which=:y)

        # data points
        scatter!(X, err_plot; markersize = 30, color = :black)

        # best-fit power law on log–log of error vs 1/k:
        # err ≈ a * (1/k)^(3p)
        order   = sortperm(X)
        Xs      = X[order]
        # anchor line at densest nonzero error
        anchor_idx = findlast(>(0), err)
        anchor_idx = isnothing(anchor_idx) ? length(err) : anchor_idx
        Xref    = X[anchor_idx]
        yref    = err_plot[anchor_idx]
        yfit    = yref .* (Xs ./ Xref).^(3p)

        tol = get_tol(key, T)
        l = lines!(Xs, yfit; color = "#9803fc", linestyle = :dash, linewidth = 5)
        # h = hlines!(ax, [tol]; color = "#ffae21", linewidth = 5)

        # textbox with p
        # txt = "p = $(round(p, digits=3))"
        # # place near top-left of current data range
        # x_text = minimum(X) * 1.15
        # y_text = maximum(err_plot) / 1.7
        # text!(ax, txt; position = (x_text, y_text), align = (:left, :center), fontsize = 28, color = :black)

        axislegend(ax, [l], ["Power Law Fit"], position = :lt, labelsize = 35, 
                orientation = :horizontal, framevisible = false, nbanks = 2, labelhalign = :center, 
                colgap = 25, patchlabelgap = 12)

        save(joinpath(outpath,"$(pot)_harmonic_$(key)_kmesh_convergence_T$(T).svg"), f)
        save(joinpath(outpath,"$(pot)_harmonic_$(key)_kmesh_convergence_T$(T).png"), f)
    end
end


"""
    fit_Qinf_vs_k!(ks, Qk; use_last=5, p_grid=0.5:0.01:2.0)

Given meshes `ks` (for k×k×k) and values `Qk` at a fixed temperature,
fit Q(N) = Q_inf + a * N^(-p) with N = k^3.

We scan p over `p_grid`. For each p, we do a linear least-squares fit of
Q ≈ c0 + c1 * x with x = 1 / N^p using normal equations (X'X)^-1 X'y.
Only the densest `use_last` meshes are used (asymptotic regime).

Returns NamedTuple:
  (Qinf, a, p_best, c1, sse, dof, Qfit, x, idx)

where Qfit are fitted values at the provided ks, x = 1/N^p_best,
and idx are the indices of the points actually used in the fit.
"""
function fit_Qinf_vs_k!(ks::AbstractVector, Qk::AbstractVector;
                        use_last::Int=5, p_grid=0.5:0.01:2.0)
    @assert length(ks) == length(Qk) "ks and Qk must have the same length"
    n = length(ks)
    m = clamp(use_last, 3, n)
    idx = (n - m + 1):n               # use densest m points (assumes ks sorted ascending)
    ksub = float.(ks[idx])
    Qsub = float.(Qk[idx])

    Nsub = ksub .^ 3
    best = nothing
    best_sse = Inf
    best_store = nothing

    for p in p_grid
        x = 1.0 ./ (Nsub .^ p)         # predictor
        # Design matrix for intercept + slope
        X = hcat(ones(length(x)), x)
        # Normal equations: θ = (X'X)^-1 X' y
        XtX = transpose(X) * X
        XtY = transpose(X) * Qsub
        θ = inv(XtX) * XtY             # θ = [c0; c1]
        c0, c1 = θ
        residuals = Qsub .- (X * θ)
        sse = sum(residuals.^2)
        if sse < best_sse
            best_sse = sse
            best = (Qinf=c0, a=c1, p_best=p)
            best_store = (x_fit = 1.0 ./ (float.(ks) .^ 3) .^ p,  # x over all ks
                          Qfit = c0 .+ c1 .* (1.0 ./ (float.(ks) .^ 3) .^ p))
        end
    end

    dof = m - 2
    return (Qinf = best.Qinf,
            a    = best.a,
            p_best = best.p_best,
            c1   = best.a,
            sse  = best_sse,
            dof  = dof,
            Qfit = best_store.Qfit,
            x    = best_store.x_fit,
            idx  = collect(idx))
end

# Force log10 ticks at integer powers within [ymin, ymax]
function set_log10_integer_power_ticks!(ax; ymin::Float64, ymax::Float64, which::Symbol=:y)
    @assert ymin > 0 "ymin must be > 0 on a log scale"
    lo = floor(Int, log10(ymin))
    hi = ceil(Int,  log10(ymax))
    vals   = 10.0 .^ (lo:hi)
    labels = [L"10^{%$(p)}" for p in lo:hi]
    if which === :y
        ax.yticks = (vals, labels)
    elseif which === :x
        ax.xticks = (vals, labels)
    else
        error("which must be :x or :y")
    end
end