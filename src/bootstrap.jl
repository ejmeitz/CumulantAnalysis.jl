function calculate_cumulants(V, V₂, V₃, V₄, V_ref, T, ce::AnalyticalEstimator, use_hot::Bool)
    c0 = CumulantData(V, V₂, V₃, V₄, V_ref, T, Val{0}(), ce, use_hot = use_hot)
    return constant_corrections(c0, T)
end

# F0/S0/U0/Cv0 should be per atom
function bootstrap_corrections(
        V, V₂, V₃, V₄, V_ref, T,
        F₀, S₀, U₀, Cᵥ₀,
        ac, # analytical corrections
        ce::AnalyticalEstimator, 
        Nat::Int, ::Type{L},
        use_hot
    ) where {L <: Limit}

    F_const, S_const, U_const, Cv_const = 
        calculate_cumulants(V, V₂, V₃, V₄, V_ref, T, ce, use_hot)

    kBNat = CumulantAnalysis.kB * Nat

    # non-dimensionalize
    S₀ /= CumulantAnalysis.kB
    Cᵥ₀ /= CumulantAnalysis.kB

    F_total_point = F₀ + (F_const/Nat) + ac.F3 + ac.F4
    S_total_point = S₀ + (S_const/kBNat) + ac.S3 + ac.S4
    U_total_point = U₀ + (U_const/Nat) + ac.U3 + ac.U4
    Cᵥ_total_point = Cᵥ₀ + (Cv_const/kBNat) + ac.Cv3 + ac.Cv4

    # Estimate standard error by bootstrapping
    is = zeros(Int, length(V))
    ΔFs = zeros(ce.n_boot); ΔSs = zeros(ce.n_boot)
    ΔUs = zeros(ce.n_boot); ΔCᵥs = zeros(ce.n_boot)

    p = Progress(ce.n_boot, "Bootstrapping Corrections")
    for i in 1:ce.n_boot
        sample!(1:length(V), is; replace = true)
        ΔFs[i], ΔSs[i], ΔUs[i], ΔCᵥs[i] =
             calculate_cumulants(V[is], V₂[is], V₃[is], V₄[is], V_ref[is], T, ce, use_hot)
        next!(p)
    end
    finish!(p)

    # Sampling error is only from the constant corrections
    F_SE = std(ΔFs) / Nat; S_SE = std(ΔSs) / kBNat
    U_SE = std(ΔUs) / Nat; Cv_SE = std(ΔCᵥs) / kBNat

    F = BootstrapCumualantEstimate(
        F₀, SVector(F_const / Nat, ac.F4, ac.F3), SVector(F_SE, 0.0, 0.0),
        F_total_point, F_SE, "F", "[eV/atom]"
    )

    S = BootstrapCumualantEstimate(
        S₀, SVector(S_const / kBNat, ac.S4, ac.S3),
        SVector(S_SE, 0.0, 0.0), S_total_point,
        S_SE, "S", "[kB / atom]"
    )

    U = BootstrapCumualantEstimate(
        U₀, SVector(U_const / Nat, ac.U4, ac.U3), SVector(U_SE, 0.0, 0.0),
        U_total_point, U_SE, "U", "[eV/atom]"
    )

    Cᵥ = BootstrapCumualantEstimate(
        Cᵥ₀, SVector(Cv_const / kBNat, ac.Cv4, ac.Cv3),
        SVector(Cv_SE, 0.0, 0.0) ./ kBNat, Cᵥ_total_point,
        Cv_SE, "Cv", "[kB / atom]"
    )

    return F, S, U, Cᵥ

end

# Potentially useful for gauging convergence of different approaches
# Bootstrap estimates error on kappa and its derivatives for all orders
function do_size_study(ce::AnalyticalEstimator, outpath, V, V₂, V₃, V₄, V_ref, T, n_atoms, use_hot)

    min_samples = (length(V) < 500) ? 10 : 100

    lg_pts = range(log10(min_samples), log10(length(V)), length = 12)
    Ns = round.(Int, 10 .^ lg_pts)

    κs = zeros(length(Ns), ce.n_boot)
    ∂κs = zeros(length(Ns), ce.n_boot)
    ∂²κs = zeros(length(Ns), ce.n_boot)

    p = Progress(length(Ns) * ce.n_boot, "Sampling Study")

    κ_point = zeros(length(Ns))
    ∂κ_point = zeros(length(Ns))
    ∂²κ_point = zeros(length(Ns))
    
    for (i,N) in enumerate(Ns)
        # pre-allocate things
        idxs = zeros(Int, N)

        # Samples are IID so our "fake" smaller dataset we'll
        # just take as the first N samples 
        V_sub  = @views V[1:N]
        V₂_sub = @views V₂[1:N]
        V₃_sub = @views V₃[1:N]
        V₄_sub = @views V₄[1:N]
        V_ref_sub = @views V_ref[1:N]

        # Get Point Estimate of Mean
        c0 = CumulantData(V_sub, V₂_sub, V₃_sub, V₄_sub, V_ref_sub, T, Val{0}(), ce, use_hot = use_hot)
      
        κ_point[i] = c0.κ
        ∂κ_point[i] = c0.∂κ_∂T
        ∂²κ_point[i] = c0.∂²κ_∂T²

        # Do bootstrap to estimate standard error
        for j in 1:ce.n_boot

            sample!(1:N, idxs; replace = true)

            V_samples = V_sub[idxs]
            V₂_samples = V₂_sub[idxs]
            V₃_samples = V₃_sub[idxs]
            V₄_samples = V₄_sub[idxs]
            V_ref_samples = V_ref_sub[idxs]

            c0 = CumulantData(V_samples, V₂_samples, V₃_samples, V₄_samples, V_ref_samples, T, Val{0}(), ce, use_hot = use_hot)
            κs[i, j] = c0.κ
            ∂κs[i, j] = c0.∂κ_∂T
            ∂²κs[i, j] = c0.∂²κ_∂T²

            next!(p)
        end
    end
    finish!(p)

    #! TODO NON-DIMENSONALIZE
    β = 1 / (kB*T)

    κ_estimates = κ_point
    ∂κ_estimates = ∂κ_point
    ∂²κ_estimates = ∂²κ_point

    κ_SEs = std(κs; dims = 2)
    ∂κ_SEs = std(∂κs; dims = 2)
    ∂²κ_SEs = std(∂²κs; dims = 2)

    data_fmt_str = (N) -> Printf.Format("%7d"*join(fill("%15.8f", N), " "))
    d_fmt = data_fmt_str(6)
    str_fmt_str = (N) -> Printf.Format("%7s"*join(fill("%15s", N-1), " "))

    header = ["N" "k" "k_SE" "dk_dT" "dk_dT_SE" "d2k_dT2" "d2k_dT2_SE"]

    open(joinpath(outpath, "outfile.nsamples_study_order0"), "w") do f
        println(f, "# Standard Error estimated from $(ce.n_boot) bootstraps of size N from origianl dataset which had $(length(V)) samples")
        println(f, "# Temperature $(T), N_atoms $(n_atoms)")
        println(f, Printf.format(str_fmt_str(length(header)), header...))
        for i in eachindex(Ns)
            println(f, Printf.format(d_fmt, Ns[i], κ_estimates[i], κ_SEs[i, 1],
                                                    ∂κ_estimates[i], ∂κ_SEs[i, 1], 
                                                    ∂²κ_estimates[i], ∂²κ_SEs[i, 1]))
        end
    end

end

# function calculate_cumulants(V, V₂, V₃, V₄, T, n_atoms, ce::SamplingCumulantEstimator{O}) where O

#     ΔF = zeros(O+1); ΔS = zeros(O+1)
#     ΔU = zeros(O+1); ΔCᵥ = zeros(O+1)

#     c0 = CumulantData(V, V₂, V₃, V₄, T, n_atoms, Val{0}(), ce)
#     ΔF[1], ΔS[1], ΔU[1], ΔCᵥ[1] = constant_corrections(c0, T)

#     if O >= 1
#         c1 = CumulantData(V, V₂, V₃, V₄, T, n_atoms, Val{1}(), ce)
#         ΔF[2], ΔS[2], ΔU[2], ΔCᵥ[2] = first_order_corrections(c1, T) 
#     end

#     if O >= 2
#         c2 = CumulantData(V, V₂, V₃, V₄, T, n_atoms, c1, Val{2}(), ce)
#         ΔF[3], ΔS[3], ΔU[3], ΔCᵥ[3] = second_order_corrections(c2, T, true)
#     end

#     return ΔF, ΔS, ΔU, ΔCᵥ

# end


# # F0/S0/U0/Cv0 should be per atom
# function bootstrap_corrections(
#         V, V₂, V₃, V₄, T, outpath,
#         F₀, S₀, U₀, Cᵥ₀,
#         ce::SamplingCumulantEstimator{O}, 
#         Nat::Int, ::Type{L}
#     ) where {O, L <: Limit}


#     # Get point estimate of corrections and control variate coefficients
#     ΔF, ΔS, ΔU, ΔCᵥ = calculate_cumulants(V, V₂, V₃, V₄, T, Nat, ce)
#     F_total_point = sum(ΔF) + (F₀*Nat)
#     S_total_point = sum(ΔS) + (S₀*Nat)
#     U_total_point = sum(ΔU) + (U₀*Nat)
#     Cᵥ_total_point = sum(ΔCᵥ) + (Cᵥ₀*Nat)

#     # Estimate standard error by bootstrapping
#     is = zeros(Int, length(V))
#     ΔFs = zeros(O+1, ce.n_boot); ΔSs = zeros(O+1, ce.n_boot)
#     ΔUs = zeros(O+1, ce.n_boot); ΔCᵥs = zeros(O+1, ce.n_boot)

#     # Re-use control variate coefficients from point estimates
#     p = Progress(ce.n_boot, "Bootstrapping Corrections")
#     for i in 1:ce.n_boot
#         sample!(1:length(V), is; replace = true)
#         ΔFs[:,i], ΔSs[:,i], ΔUs[:,i], ΔCᵥs[:,i] =
#              calculate_cumulants(V[is], V₂[is], V₃[is], V₄[is], T, Nat, ce)
#         next!(p)
#     end
#     finish!(p)

#     F_totals = sum(ΔFs, dims = 1) .+ (F₀*Nat)
#     S_totals = sum(ΔSs, dims = 1) .+ (S₀*Nat)
#     U_totals = sum(ΔUs, dims = 1) .+ (U₀*Nat)
#     Cᵥ_totals = sum(ΔCᵥs, dims = 1) .+ (Cᵥ₀*Nat)
    
#     F_SEs = std(ΔFs, dims = 2); S_SEs = std(ΔSs, dims = 2)
#     U_SEs = std(ΔUs, dims = 2); Cᵥ_SEs = std(ΔCᵥs, dims = 2)

#     kBNat = CumulantAnalysis.kB * Nat

#     F = BootstrapCumualantEstimate(
#         F₀, SVector(ΔF...) ./ Nat, SVector(F_SEs...) ./ Nat,
#         F_total_point / Nat, std(F_totals) / Nat, "F", "[eV/atom]"
#     )

#     S = BootstrapCumualantEstimate(
#         S₀ / CumulantAnalysis.kB, SVector(ΔS...) ./ kBNat, SVector(S_SEs...) ./ kBNat,
#         S_total_point / kBNat, std(S_totals) / kBNat, "S", "[kB / atom]"
#     )

#     U = BootstrapCumualantEstimate(
#         U₀, SVector(ΔU...) ./ Nat, SVector(U_SEs...) ./ Nat,
#         U_total_point / Nat, std(U_totals) / Nat, "U", "[eV/atom]"
#     )

#     Cᵥ = BootstrapCumualantEstimate(
#         Cᵥ₀ / CumulantAnalysis.kB, SVector(ΔCᵥ...) ./ kBNat, SVector(Cᵥ_SEs...) ./ kBNat,
#         Cᵥ_total_point / kBNat, std(Cᵥ_totals) / kBNat, "Cv", "[kB / atom]"
#     )

#     return F, S, U, Cᵥ

# end


# # Potentially useful for gauging convergence of different approaches
# # Bootstrap estimates error on kappa and its derivatives for all orders
# function do_size_study(ce::SamplingCumulantEstimator{O}, outpath, V, V₂, V₃, V₄, T, n_atoms) where O

#     min_samples = (length(V) < 500) ? 10 : 100

#     lg_pts = range(log10(min_samples), log10(length(V)), length = 12)
#     Ns = round.(Int, 10 .^ lg_pts)

#     κs = zeros(length(Ns), O+1, ce.n_boot)
#     ∂κs = zeros(length(Ns), O+1, ce.n_boot)
#     ∂²κs = zeros(length(Ns), O+1, ce.n_boot)

#     p = Progress(length(Ns) * ce.n_boot, "Sampling Study")

#     κ_point = zeros(length(Ns), O+1)
#     ∂κ_point = zeros(length(Ns), O+1)
#     ∂²κ_point = zeros(length(Ns), O+1)
    
#     for (i,N) in enumerate(Ns)
#         # pre-allocate things
#         idxs = zeros(Int, N)

#         # Samples are IID so our "fake" smaller dataset we'll
#         # just take as the first N samples 
#         V_sub  = @views V[1:N]
#         V₂_sub = @views V₂[1:N]
#         V₃_sub = @views V₃[1:N]
#         V₄_sub = @views V₄[1:N]

#         # Get Point Estimate of Mean
#         c0 = CumulantData(V_sub, V₂_sub, V₃_sub, V₄_sub, T, n_atoms, Val{0}(), ce)
#         c1 = CumulantData(V_sub, V₂_sub, V₃_sub, V₄_sub, T, n_atoms, Val{1}(), ce)
#         c2 = CumulantData(V_sub, V₂_sub, V₃_sub, V₄_sub, T, n_atoms, c1, Val{2}(), ce)
#         cds = (c0, c1, c2)
#         for co in 0:O         
#             κ_point[i, co + 1] = cds[co + 1].κ
#             ∂κ_point[i, co + 1] = cds[co + 1].∂κ_∂T
#             ∂²κ_point[i, co + 1] = cds[co + 1].∂²κ_∂T²
#         end

#         # Do bootstrap to estimate standard error
#         for j in 1:ce.n_boot

#             sample!(1:N, idxs; replace = true)

#             V_samples = V_sub[idxs]
#             V₂_samples = V₂_sub[idxs]
#             V₃_samples = V₃_sub[idxs]
#             V₄_samples = V₄_sub[idxs]

#             c0 = CumulantData(V_samples, V₂_samples, V₃_samples, V₄_samples, T, n_atoms, Val{0}(), ce)
#             c1 = CumulantData(V_samples, V₂_samples, V₃_samples, V₄_samples, T, n_atoms, Val{1}(), ce)
#             c2 = CumulantData(V_samples, V₂_samples, V₃_samples, V₄_samples, T, n_atoms, c1, Val{2}(), ce)

#             cds = (c0, c1, c2)

#             for co in 0:O         
#                 κs[i, co + 1, j] = cds[co + 1].κ
#                 ∂κs[i, co + 1, j] = cds[co + 1].∂κ_∂T
#                 ∂²κs[i, co + 1, j] = cds[co + 1].∂²κ_∂T²
#             end
#             next!(p)
#         end
#     end
#     finish!(p)

#     #! TODO NON-DIMENSONALIZE
#     β = 1 / (kB*T)

#     κ_estimates = κ_point
#     ∂κ_estimates = ∂κ_point
#     ∂²κ_estimates = ∂²κ_point

#     κ_SEs = std(κs; dims = 3)
#     ∂κ_SEs = std(∂κs; dims = 3)
#     ∂²κ_SEs = std(∂²κs; dims = 3)

#     data_fmt_str = (N) -> Printf.Format("%7d"*join(fill("%15.8f", N), " "))
#     d_fmt = data_fmt_str(6)
#     str_fmt_str = (N) -> Printf.Format("%7s"*join(fill("%15s", N-1), " "))

#     for co in 0:O
#         header = ["N" "k" "k_SE" "dk_dT" "dk_dT_SE" "d2k_dT2" "d2k_dT2_SE"]

#         open(joinpath(outpath, "outfile.nsamples_study_order$(co)"), "w") do f
#             println(f, "# Standard Error estimated from $(ce.n_boot) bootstraps of size N from origianl dataset which had $(length(V)) samples")
#             println(f, "# Temperature $(T), N_atoms $(n_atoms)")
#             println(f, Printf.format(str_fmt_str(length(header)), header...))
#             for i in eachindex(Ns)
#                 println(f, Printf.format(d_fmt, Ns[i], κ_estimates[i, co+1], κ_SEs[i, co+1],
#                                                        ∂κ_estimates[i, co+1], ∂κ_SEs[i, co+1], 
#                                                        ∂²κ_estimates[i, co+1], ∂²κ_SEs[i, co+1]))
#             end
#         end
#     end

# end