

function load_tdep_ifc2()

end

function get_modes(ifc2::Matrix{T}, masses::AbstractVector{M}; D::Int = 3)

    dynmat = zeros(size(ifc2))
    N_atoms = Int(size(ifc2,1) / D)

    Threads.@threads for i in range(1, N_atoms)
        for j in range(1, N_atoms)
            for α in range(1,D)
                for β in range(1,D)
                    ii = D*(i-1) + α
                    jj = D*(j-1) + β
                    dynmat[ii,jj] /=  ustrip(sqrt(masses[i]*masses[j]))
                end
            end
        end
    end

    dynmat_units = ifc2.units / unit(masses[1])

    return get_modes(dynmat * dynmat_units, D)

end

function get_modes(dynmat::Matrix{T}, D::Int = 3)

    dynmat_unit = unit(dynmat[1,1])
    dynmat_unitless = ustrip(dynmat)

    eig_stuff = eigen(Hermitian(dynmat_unitless))
    freqs_sq = eig_stuff.values
    idx_rt = sortperm(abs.(freqs_sq))
    # Zero out rigid translation modes
    freqs_sq[idx_rt[1:D]] .= 0.0

    return freqs_sq*dynmat_unit, eig_stuff.vectors
end

# Assumes output has been re-mapped to the supercell
function parse_TDEP_second_order(ifc_path::String, N_modes, energy_units::Symbol = :NONE)

    Φ = zeros(N_modes, N_modes)

    open(ifc_path, "r") do f
    
        #First 2 lines are header
        N_atoms = parse(Int64,split(strip(readline(f)))[1])
        @assert 3*N_atoms == N_modes "Cannot handle non-monatomic case"
        r_cut = parse(Float64,split(strip(readline(f)))[1])
    
        #Next lines gives num neighbors of atom in unit cell we are looking at
        Φ_block = zeros(3,3)
        for i in 1:N_atoms
            n_neighbors, base_atom_idx_uc = parse.(Int64,split(strip(readline(f)))[[1,end-1]])
            
            for j in 1:n_neighbors
                other_atom_idx_uc = parse(Int64,split(strip(readline(f)))[1])
                lattice_vec = parse.(Float64,split(strip(readline(f)))) #Dont think this is useful?

                Φ_block[1,:] .= parse.(Float64,split(strip(readline(f)))) 
                Φ_block[2,:] .= parse.(Float64,split(strip(readline(f))))
                Φ_block[3,:] .= parse.(Float64,split(strip(readline(f))))

                Φ[3*(base_atom_idx_uc - 1) + 1 : 3*base_atom_idx_uc,
                  3*(other_atom_idx_uc - 1) + 1 : 3*other_atom_idx_uc] .= Φ_block

            end
        end
    end

    if energy_units == :REAL
        Φ .*= 23.060541945
        units = u"kcal * mol^-1 * Å^-2"
    elseif energy_units == :METAL
        units = u"eV * Å^-2"
    elseif energy_units == :NONE
        units = NoUnits
    else
        throw(ArgumentError("Unsupported unit type: $(energy_units)"))
    end

    return DenseForceConstants(Φ, units, 0.0)

end