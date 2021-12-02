module WingBiology

    using LinearAlgebra
    using ForwardDiff
    using SparseArrays
    using SparseDiffTools
    using LinearMaps
    using Arpack

    using PyCall

    include("FuselageAero.jl")
    include("WingAero.jl")
    include("Engine.jl")

    include("StructuralAnalysis.jl")

    include("DataTypes.jl")
    include("Plotting.jl")

    include("StateSpace.jl")
    include("Newton.jl")
    include("Arnoldi.jl")
    include("AssumedModes.jl")

end # module
