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

    _linv_reg(A::AbstractMatrix, λ::Real = 1e-3) = let D = A' * A
        inv(
            D .+ Matrix{Float64}(I * λ, size(D)...)
        ) * A'
    end

    include("Arnoldi.jl")
    include("AssumedModes.jl")

end # module
