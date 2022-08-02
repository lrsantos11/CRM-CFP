
using QPSReader
using Pkg.Artifacts

"Return the path to the Netlib linear already  optimization test set."
function fetch_netlib_presolved()
    artifact_toml = joinpath(@__DIR__, "..", "Artifacts.toml")
    ensure_artifact_installed("netliblp_presolved", artifact_toml)
    netliblp_presolved_hash = artifact_hash("netliblp_presolved", artifact_toml)
    @assert artifact_exists(netliblp_presolved_hash)
    return joinpath(artifact_path(netliblp_presolved_hash), "MPSpre")
end

const NETLIB_PRE_DIR = fetch_netlib_presolved()
const NETLIB_DIR = fetch_netlib()

##

function LPtoSTDFormat(c, A, l, u, xlb, xub)
    nrow, ncol = size(A)
    b = Array{Float64,1}()
    # Transform l <= Ax <= u into Ax = b
    delrows = Array{Int64,1}()
    for row = 1:nrow
        if l[row] > u[row]
            throw(error("Problem is infeasible."))
        elseif l[row] == -Inf && u[row] == Inf #Constraint is always feasible
            push!(delrows, row)
            push!(b, Inf)      #Creates dummy b[row] just to delete at the end.        
        elseif l[row] == u[row] # Constraint is l = a'x = u        
            push!(b, l[row])
        elseif l[row] > -Inf && u[row] == Inf #Constraint is  a'x >= l
            ncol += 1
            A = [A spzeros(nrow)]
            A[row, end] = -1.0 # a'x - xs = l
            push!(b, l[row]) # b <- l 
            push!(c, 0.0) # no cost
            push!(xlb, 0.0) #xs >= 0
            push!(xub, Inf) #xs <= Inf
        elseif l[row] == -Inf && u[row] < Inf # Constraint is  a'x <= u
            ncol += 1
            A = [A spzeros(nrow)]
            A[row, end] = 1.0 # a'x + xs = u
            push!(b, u[row]) # b <- u
            push!(c, 0.0) # no cost
            push!(xlb, 0.0) #xs >= 0
            push!(xub, Inf) #xs <= Inf            
        elseif l[row] > -Inf && u[row] < Inf # Constraint is  l <a'x < u.
            A = [A spzeros(nrow)]
            A[row, end] = -1.0 # a'x = xs
            push!(b, 0.0) # b <-
            push!(c, 0.0)
            push!(xlb, l[row]) # adds variable
            push!(xub, u[row])
            ncol += 1
        end
    end
    A = A[setdiff(1:end, delrows), :]
    b = deleteat!(b, delrows)
    # Transform xlb <= x <= xub  into 0 <= x <= xub
    delcols = Array{Int64,1}()
    for col = 1:ncol
        if xlb[col] > xub[col]
            throw(error("Problem is infeasible."))
        elseif xlb[col] == -Inf && xub[col] == Inf #Free variable
            A = [A -A[:, col]] #x_i = xp - xm
            push!(c, -c[col]) # adds cost for xm
            xlb[col] = 0.0
            push!(xlb, 0.0) #xm >= 0
            push!(xub, Inf) #xs <= Inf
        elseif xlb[col] == xub[col] # Constraint is l = x = u        
            b -= xlb[col] * A[:, col]
            push!(delcols, col)
        elseif xlb[col] > -Inf && xub[col] <= Inf
            b -= xlb[col] * A[:, col]
            xub[col] -= xlb[col]
            xlb[col] = 0.0
        elseif xlb[col] == -Inf && xub[col] < Inf
            c[col] = -c[col]
            A[:, col] = -A[:, col]
            b = b - xub[col] * A[:, col]
            xlb[col] = 0.0
            xub[col] = Inf
        end
    end
    A = A[:, setdiff(1:end, delcols)]
    c = c[setdiff(1:end, delcols)]
    xlb = xlb[setdiff(1:end, delcols)]
    xub = xub[setdiff(1:end, delcols)]
    nrow, ncol = size(A)
    return nrow, ncol, c, A, b, xlb, xub
end


##

function LPDataFromMPS(mpsfile::String; stdformat::Bool = true)
    lp_model = readqps(mpsfile)
    cc = lp_model.c
    A = sparse(lp_model.arows, lp_model.acols, lp_model.avals)
    l = lp_model.lcon
    u = lp_model.ucon
    xlb = lp_model.lvar
    xub = lp_model.uvar
    if stdformat
        return LPtoSTDFormat(cc, A, l, u, xlb, xub)
    else
        nrow, ncol = size(A)
        return nrow, ncol, cc, A, l, u, xlb, xub
    end
end

