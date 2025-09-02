

function parse_energies(path, n_atoms)

    data = readdlm(path, comments = true)
    
    @views E_polar = data[:, 2]
    @views E_pair = data[:, 3]
    @views E_triplet = data[:, 4]
    @views E_quartet = data[:, 5]
    
    # energies in file are meV / atom
    # conver to eV

    return E_polar, E_pair, E_triplet, E_quartet
end

function estimate_cumulants(V2, V3, V4, β)
    κ₁ = mean(V4)
    κ₂ = (-β/2) * mean(V3.^2)
end

function cumulants_from_effective_hamiltonian(
        outpath::String;
        ucposcar_path::String = joinpath(outpath, "infile.ucposcar"),
        ifc2_path::String = joinpath(outpath, "infile.ucposcar")
    )

    isfile(ucposcar_path) || throw(ArgumentError("ucposcar_path is not a file: $(ucposcar_path)"))
    isfile(ifc2_path) || throw(ArgumentError("Force constant path is not a file: $(ifc2_path)"))

    cp(ifc_path, joinpath(outpath, "infile.forceconstant"); force = true)
    cp(ucposcar_path, joinpath(outpath, "infile.ucposcar"); force = true)

    #! TODO RUN EFFECTIVE HAMILTONIAN

    # energies in meV / atom
    Vp, V2, V3, V4 = parse_energies(joinpath(outpath, "outfile.energies"))

    if !(sum(Vp) ≈ 0.0)
        error("Got non-zero polar term. I don't know what to do with that.")
    end

    V = V2 .+ V3 .+ V4



end