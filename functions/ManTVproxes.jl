using Manopt

function ManTVproxDual(N::TangentBundle,X::TBPoint,Ξ::TBTVector,anisotropic=true)
    M = N.manifold
    m = getBase(X)
    ξ = getTangent(Ξ)
    d = M.powerSize
    dims = length(d)
    N₁ = Power(M.manifold,(d[end],))

    if dims==1
        R = CartesianIndices(d)
        Kmax = 1
    else
        R = CartesianIndices(d[1:end-1])
        Kmax = dims-1
    end

    η = zeroTVector(M,m)

    for i in R
        if anisotropic || dims==1
            for k in 1:Kmax
                g = norm(M.manifold,m[i,k],ξ[i,k])
                η[i,k] = ξ[i,k]/max(1,g)
            end
        else
            g = norm(N₁,m[i,:],ξ[i,:])
            for k in 1:Kmax
                η[i,k] = ξ[i,k]/max(1,g)
            end
        end
    end
    return TBTVector(getBase(Ξ),η)
end
