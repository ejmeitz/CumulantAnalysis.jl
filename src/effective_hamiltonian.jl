

function parse_energies(path)

    data = readdlm(path, comments = true)

    # Parse n_atoms and T from header
    f = open(path, "r")
    readline(f) # skip
    T = parse(Float64, split(strip(readline(f)))[end])
    n_atoms = parse(Int, split(strip(readline(f)))[end])
    close(f)

    # energies in file are meV / atom, Converts to eV
    conv = n_atoms / 1000
    
    @views E_polar = data[:, 2] .* conv
    @views E_pair = data[:, 3] .* conv
    @views E_triplet = data[:, 4] .* conv
    @views E_quartet = data[:, 5] .* conv

    return T, n_atoms, E_polar, E_pair , E_triplet, E_quartet
end


function cumulants_from_effective_hamiltonian(
        T::Real, # kelvin
        ::Type{L}, 
        outpath::String;
        order::Int = 3,
        energies_file::Union{String, Nothing} = nothing,
        nconf::Int = 50_000,
        ucposcar_path::String = joinpath(outpath, "infile.ucposcar"),
        ssposcar_path::String = joinpath(outpath, "infile.ssposcar"),
        ifc_path::String = joinpath(outpath, "infile.forceconstant"),
        n_boot::Int = 100,
        boot_size::Int = floor(Int, nconf / 5)
    ) where {L <: Limit}

    isfile(ucposcar_path) || throw(ArgumentError("ucposcar path is not a file: $(ucposcar_path)"))
    isfile(ifc_path) || throw(ArgumentError("Force constant path is not a file: $(ifc_path)"))
    isfile(ssposcar_path) || throw(ArgumentError("ssposcar path is not a file: $(ssposcar_path)"))


    new_ifc_path = joinpath(outpath, "infile.forceconstant")
    new_uc_path = joinpath(outpath, "infile.ucposcar")
    new_ss_path = joinpath(outpath, "infile.ssposcar")
    isfile(new_ifc_path) || cp(ifc_path, new_ifc_path; force = true)
    isfile(new_uc_path) || cp(ucposcar_path, new_uc_path; force = true)
    isfile(new_ss_path) || cp(ssposcar_path, new_ss_path; force = true)

    if isnothing(energies_file)
        #! TODO RUN EFFECTIVE HAMILTONIAN
        error("Not implemented yet")
        energies_file = joinpath(outpath, "outfile.energies")
    end

    # energies in meV / atom
    T_file, n_atoms, Vp, V2, V3, V4 = parse_energies(energies_file)

    if T_file != T
        @warn "You said tempearture was $T, but parsed temperature as $(T_file) from outfile.energies"
    end

    if !(sum(Vp) ≈ 0.0)
        error("Got non-zero polar term. I don't know what to do with that.")
    end

    ΔV = V3 .+ V4

    #! We sampled w.r.t. harmonic so derivatives take V2 as the "reference"
    return bootstrap_corrections(T, V2, ΔV, n_boot, boot_size, outpath, n_atoms, order, L)

end