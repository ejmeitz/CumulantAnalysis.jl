abstract type FCData{O,T} end
abstract type AtomFC{O,T,N} end

# mimics the lo_fc2_pair type
struct FC2Data{T} <: FCData{2,T}
    idxs::SVector{2, Int}
    lvs::SVector{2, SVector{3, T}} # lattice vectors for the unit cell these atoms belong to
    r::SVector{3, T} # vector between atom 1 and 2
    ifcs::SMatrix{3, 3, T}
end 

# mimics the lo_fc3_triplet type
struct FC3Data{T} <: FCData{3,T}
    idxs::SVector{3, Int}
    lvs::SVector{3, SVector{3, T}}
    rv1::SVector{3, T}
    rv2::SVector{3, T}
    rv3::SVector{3, T}
    ifcs::SArray{Tuple{3,3,3}, T}
end

struct FC4Data{T} <: FCData{4,T}
    idxs::SVector{4, Int}
    lvs::SVector{4, SVector{3, T}}
    rv1::SVector{3, T}
    rv2::SVector{3, T}
    rv3::SVector{3, T}
    rv4::SVector{4,T}
    ifcs::SArray{Tuple{3,3,3,3}, T}
end

# mimics the lo_fc2_atom type
struct AtomFC2{T,N} <: AtomFC{2,T,N}
    pairs::SVector{N, FC2Data{T}}
end

# mimics the lo_fc3_atom type
struct AtomFC3{T,N} <: AtomFC{3,T,N}
    triplets::SVector{N, FC3Data{T}}
end

# mimics the lo_fc4_atom type
struct AtomFC4{T,N} <: AtomFC{4,T,N}
    quartets::SVector{N, FC4Data{T}}
end

n_neighbors(::AtomFC{O,T,N}) where {O,T,N} = N

struct IFCS{O, T}
    na_uc::Int # number of atoms in unit cell
    r_cut::T
    atoms::AbstractVector{<:AtomFC{O,T}} # this is super ragged...
end

# Helper functions
readline_skip_text!(io, T) = parse(T, first(split(strip(readline(io)))))

function read_vec3!(io, T)
    xs = split(strip(readline(io)))
    @assert length(xs) == 3 "Expected 3 components for a 3-vector."
    SVector{3,T}(parse.(T, xs))
end

"""
    read_fc2_pairs(path, R_frac, A)

Read a 2nd-order outfile.forceconstant from TDEP. This follows the logic implemented in TDEP.
Does NOT handle polar force constants.

Inputs
- `path::AbstractString`: file path
- `r_frac_uc::AbstractVector{SVector{3,T}}`: 3×na fractional coords of atoms in the unitcell
- `A::AbstractMatrix{T}`: 3×3 lattice (columns are lattice vectors), Cartesian = A * fractional
- `chop_tol::T`: small zeroing threshold (like `lo_chop`)

"""
function read_fc2_pairs(path::AbstractString,
                        r_frac_uc::AbstractVector{SVector{3,T}},
                        A::AbstractMatrix{T};
                        chop_tol::T = T(1e-13)) where {T<:Real}

    
    function read_mat3_rows!(io) 
        r1 = read_vec3!(io, T)
        r2 = read_vec3!(io, T)
        r3 = read_vec3!(io, T)
        M = hcat(r1, r2, r3)'   # r1 is first row, etc.
        return SMatrix{3,3,T}(M)     
    end

    chop3(v::SVector{3,T}) = SVector{3,T}(ntuple(i->(abs(v[i]) < chop_tol ? zero(T) : v[i]), 3))

    @assert size(A) == (3,3) "Lattice A must be 3×3."

    max_rcut = zero(T)

    data = AtomFC2{T}[]

    open(path, "r") do io
        na = readline_skip_text!(io, Int) # number of atoms in unit cell

        # cross-check given positions
        @assert length(r_frac_uc) == na "r_frac_uc has length $(length(r_frac_uc)) but file says there are $(na) atoms in the unitcell."

        # recomputed based on max from whole file
        _cutoff_ignored = readline_skip_text!(io, Float64)

        for a1 in 1:na
            n_neighbors = readline_skip_text!(io, Int)
            pair_data = Vector{FC2Data{T}}(undef, n_neighbors) 
            for i in 1:n_neighbors
                a2 = readline_skip_text!(io, Int) # unitcell index of neighbor
                lv2_frac = read_vec3!(io, T)
                ifcs = read_mat3_rows!(io)

                # Fortran routine uses lv1 == 0
                lv1_frac = SVector{3,T}(0,0,0)

                # frac positions of the two atoms in their image cells
                v1_frac = lv1_frac + SVector{3,T}(r_frac_uc[a1])
                v2_frac = lv2_frac + SVector{3,T}(r_frac_uc[a2])

                # Convert to Cartesian
                lv1_cart = SVector{3,T}(A * lv1_frac)
                lv2_cart = SVector{3,T}(A * lv2_frac)
                r_cart   = SVector{3,T}(A * (v2_frac - v1_frac))

                lv1_cart = chop3(lv1_cart)
                lv2_cart = chop3(lv2_cart)
                r_cart   = chop3(r_cart)

                max_rcut = max(max_rcut, norm(r_cart))

                pair_data[i] = FC2Data{T}(
                    SVector{2,Int}(a1, a2),
                    SVector{2,SVector{3,T}}(lv1_cart, lv2_cart),
                    r_cart,
                    ifcs
                )
            end
            
            # Would be nice to not have to copy into SVector here
            push!(data, AtomFC2{T, n_neighbors}(SVector{n_neighbors}(pair_data)))
        end

        # technically theres more polar stuff, but ignore that for now

        r_cut = max_rcut + sqrt(eps(T))
        println(typeof(data))
        return IFCS{2, T}(na, r_cut, data)

    end
end

function read_ifc3(path)

end

function read_ifc4(path)

end

function remap()

end