
"""
Infostruct
"""
mutable struct Infostruct
    nK::Int64
    nA::Int64
    nMat::Int64
    nz::Int64
    nm::Int64
    nzforce::Int64
    Ffreq::Int64
    Fmodes::Vector{Int64}
    Fmult::Vector{Float64}
    nrom:: Int64
    α::Float64
    β::Float64
    neig::Int64
    Φ::Vector{Int64}
    style::Char
    max_order::Int64
    max_orderNA::Int64
    output_dir::String
    omega_mul::Float64
    Infostruct() = new()
end
info=Infostruct()


"""
Overview: set of constants used over the project.
This helps minimizing literals throughout the work.
If you need to insert constants used througout the work,
then please add them here to help minimizing unmeaningful 
numbers in the work.
Furthermore, add 
"""

"""
resonances tolerance
"""
const ϵ_tol = 1e-3

"""
Dimension of the treated problems (3D)
"""
const dim = 3

"""
Number of nodes of a vertex (1)
"""
const VE1n = 1

"""
Number of nodes of a linear line element (2)
"""
const BE2n = 2

"""
Number of nodes of a quadratic line element (3)
"""
const BE3n = 3

"""
Number of nodes of a linear triangular element (3)
"""
const TR3n = 3

"""
Number of nodes of a quadratic triangular element (6)
"""
const TR6n = 6

"""
Number of nodes of a linear tetrahedral element (4)
"""
const TE4n = 4

"""
Number of nodes of a quadratic tetrahedral element (10)
"""
const T10n = 10

"""
Number of nodes of a linear wedge element (6)
"""
const PE6n = 6

"""
Number of nodes of a quadratic 15-nodes wedge element (15)
"""
const P15n = 15

"""
Number of nodes of a quadratic 18-nodes wedge element (18)
"""
const P18n = 18

"""
Number of nodes of a linear square element (4)
"""
const QU4n = 4

"""
Number of nodes of a quadratic serendipity square element (8)
"""
const QU8n = 8

"""
Number of nodes of a quadratic Lagrange square element (9)
"""
const QU9n = 9

"""
Number of nodes of a linear hexahedral element (8)
"""
const HE8n = 8

"""
Number of nodes of a quadratic serendipity hexahedral element (20)
"""
const H20n = 20

"""
Number of nodes of a quadratic Lagrange hexahedral element (27)
"""
const H27n = 27

"""
Number of supported boundary element
"""
const Snse = 5

"""
Number of supported volume element
"""
const Vnse = 8

"""
Maximum number of elements connected one to the other
"""
const ncv_max = 30


"""
Maximum number nodes for an element
"""
const nne_max = 27

