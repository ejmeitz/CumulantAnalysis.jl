export sTDEPEstimator, LAMMPSEstimator, MollyEstimator, estimate, Quantum, Classical

abstract type Limit end
struct Quantum <: Limit end
struct Classical <: Limit end


abstract type ThermoEstimator{O <: Integer, L <: Limit}end
abstract type MDEstimator{O} <: ThermoEstimator{O, Classical} end

order(::ThermoEstimator{O}) where O = O

"""
"""
struct sTDEPEstimator{O,L,T,C} <: ThermoEstimator{O,L}
    cc::C # C will be CanonicalConfiguratoin, but that is not availbale here
    temperature::T
end

limit(::sTDEPEstimator{O,L}) where {O,L} = L
stochastic(::sTDEPEstimator) = true

"""
"""
struct LAMMPSEstimator{O,T} <: MDEstimator{O}
    temperature::T
end


"""
"""
struct MollyEstimator{O,T} <: MDEstimator{O}
    temperature::T
    n_configs::Int
end

limit(::MDEstimator) = Classical()
stochastic(::MDEstimator) = false
