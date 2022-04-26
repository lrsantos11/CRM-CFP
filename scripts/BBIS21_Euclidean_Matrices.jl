##################################################################
## Basic Functions for Euclidean Matrices Intersectio tests
##################################################################
using Combinatorics


function ProjectionPSD(A::AbstractMatrix)
    return ProjectIndicator(IndPSD(),A)
end


function createEuclideanMatrix(Points::Vector{Vector{T}}) where T
    num_points = length(Points)
    D = Matrix{T}(undef,num_points,num_points)
    for point in eachrow(Points)
        D[:,point.indices[1]] = norm.(point .- Points).^2
    end
    return D
end
function MakePoints(N::Int)
    return (p = (a,c) -> c < N-1 ? a : [Float64[x;y] for x = a for y = p(a,c-1)])(0:N-1,N)
end

N=3
Points = MakePoints(N)
E = createEuclideanMatrix(Points)
D = deepcopy(E)
D[D .>  4] .= 0

function ProjectionCompletingD(D::AbstractMatrix,
                               X::AbstractMatrix;
                               ε::Number = 0.0)
    !issymmetric(A) && error("Matrix A is not symmetric")
    num_rows, num_cols = size(D)
    proj = deepcopy(D)
    for i in 1:num_rows, j in 1:num_cols
        if D[i,j] ≈ ε && i != j
            proj[i,j] = max(A[i,j],0.0)
        end
    end
    return proj
end
Random.seed!(0)
A = randn(size(D))
A = A'*A

proj = ProjectionCompletingD(D,A)
proj = ProjectionPSD(proj)


