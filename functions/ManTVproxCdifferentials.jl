using Manopt

function ManTVproxDualCdifferential(N::TangentBundle,n::TBPoint,
    Η::TBTVector,Ξ::TBTVector,anisotropic=true)

    M = N.manifold
    m = getBase(n)
    d = M.powerSize
    ldims = length(d)
    N₁ = Power(M.manifold,(d[end],))

    if ldims==1
        Kmax = 1
        R = CartesianIndices(d)
    else
        Kmax = ldims-1
        R = CartesianIndices(d[1:Kmax])
    end
    maxIt = [last(R).I...]

    η = getTangent(Η)
    ξ = getTangent(Ξ)

    J = zeroTVector(M,m)
    for j in R
        if anisotropic || ldims==1
            for k in 1:Kmax
                mⱼ = m[j,k]
                ηⱼ = η[j,k]
                g = norm(M.manifold,mⱼ,ηⱼ)
                if !(j[k]==maxIt[k])
                    ξⱼ = ξ[j,k]
                    if g <=1
                        J[j,k] +=  ξⱼ
                    else
                        J[j,k] +=  1/g * (ξⱼ - 1/g^2 * dot(M.manifold,mⱼ,ξⱼ,ηⱼ)*ηⱼ)
                    end
                else
                    # Taking care of the boundary equations
                    # J[j,k] = zeroTVector(M.manifold,mⱼ)
                end
            end
        else
            g = norm(N₁,m[j,:],η[j,:])
            for k in 1:Kmax
                mⱼ = m[j,k]
                ηⱼ = η[j,k]
                if !(j[k]==maxIt[k])

                    ξⱼ = ξ[j,k]
                    if g <=1
                        J[j,k] +=  ξⱼ
                    else
                        for κ in 1:Kmax
                            if κ != k
                                J[j,κ] += - 1/g^3 * dot(M.manifold,mⱼ,ξⱼ,ηⱼ)*η[j,κ]
                            else
                                J[j,k] += 1/g * (ξⱼ - 1/g^2 * dot(M.manifold,mⱼ,ξⱼ,ηⱼ)*ηⱼ)
                            end
                        end
                    end
                else
                    # Taking care of the boundary equations
                    # J[j,k] = zeroTVector(M.manifold,mⱼ)
                end
            end

        end

    end

    return J

end

function ManTVproxDualCdifferential(M::Manifold,m::MPoint,d::Tuple,
    jk::Tuple{CartesianIndex,Int},η::TVector,ξ::TVector,
    anisotropic=true)

        ldims = length(d)

        j,k = jk

        N = Power(M,(ldims,))
        nn = PowPoint(fill(m,(ldims,)))

        if anisotropic || ldims==1
            ηⱼₖ = η[k]
            g = norm(M,m,ηⱼₖ)
            if !(j[k]==d[k])
                ξⱼₖ = ξ
                if g <=1
                    Jⱼₖ =  ξⱼₖ
                else
                    Jⱼₖ =  1/g * (ξⱼₖ - 1/g^2 * dot(M,m,ξⱼₖ,ηⱼₖ)*ηⱼₖ)
                end
            else
                # Taking care of the boundary equations
                Jⱼₖ = zeroTVector(M,m)
            end

        else
            g = norm(N,nn,η)
            ηⱼₖ = η[k]
            if !(j[k]==d[k])
                ξⱼₖ = ξ
                if g <=1
                    Jⱼₖ =  ξⱼₖ
                else
                    Jⱼₖ = Tuple{Int,TVector}[]
                    for κ in 1:ldims
                        if κ != k
                            push!(Jⱼₖ, tuple(κ,- 1/g^3 * dot(M,m,ξⱼₖ,ηⱼₖ)*η[κ]))
                        else
                            push!(Jⱼₖ,tuple(k, 1/g * (ξⱼₖ - 1/g^2 * dot(M,m,ξⱼₖ,ηⱼₖ)*ηⱼₖ)))
                        end
                    end
                end
            else
                # Taking care of the boundary equations
                Jⱼₖ = zeroTVector(M,m)
            end
        end

        return Jⱼₖ

end
