include(joinpath("..","jltools","ExpTools.jl"))
include(joinpath("..","jltools","Dbg.jl"))
include(joinpath("..","plans","problem.jl"))
include("lRSSN.jl")

using Manopt, SparseArrays, LinearAlgebra

function lRSSN_ML(pr::P,x0::PowPoint;y0::Union{TVector,Missing}=missing,
    return_dual::Bool=false,tol=1e-12,max_it::Int=50,
    save_iterates::Bool=false) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n
    dims = manifoldDimension(M)
    d = M.powerSize
    ldims = length(d)
    Mdims = manifoldDimension(M.manifold)
    NN = N.manifold
    dd = NN.powerSize
    dualdims = manifoldDimension(NN)
    Ndims = manifoldDimension(N.manifold.manifold)

    Λ = pr.operator

    # Setup algorithm
    p = x0
    if !ismissing(y0)
        ξ = y0
    else
        ξ = zeroTVector(N,n)
    end
    pξ = tuple(p,ξ)

    change = 1
    cost = getcost(pr,p)
    if save_iterates
        r = Tuple{Int32,Float64,Float64,PowPoint,TBTVector,Float64}[]
    else
        r = Tuple{Int32,Float64,Float64,Float64}[]
    end
    dbg(1,"# 0 | Change: - | RelError: 1. | F(p) = $(cost)")
    k = 1
    X = construct_lX(pr,pξ)
    nX0 = norm(X)
    relerror = 1.

    if dbglevel_atleast(3)
        maxit = min(20,max_it)
    else
        maxit = max_it
    end

    while k<=maxit && relerror >tol
        dbg(4,"k = $(k)")
        dbg(4,"p = $(p)")
        dbg(4,"ξ = $(ξ)")

        DpX1 = construct_lDpX1(pr,pξ)
        dbg(2,"DpX1 constructed")
        DξX1 = construct_lDξX1(pr,pξ)
        dbg(2,"DξX1 constructed")
        DpX2 = construct_lDpX2(pr,pξ)
        dbg(2,"DpX2 constructed")
        DξX2 = construct_lDξX2(pr,pξ)
        dbg(2,"DξX2 constructed")


        G = [DpX1 DξX1;
        DpX2 DξX2]
        dbg(2,"G constructed")
        dbg(5,"G = $(Array(G))")
        dbg(4,"cond(G) = $(cond(Array(G)))")
        # eig1 = eigen(Array(G))
        # dbg(1,"eigenvalues = $(eig1.values)")

        X = construct_lX(pr,pξ)
        dbg(4,"X = $(X)")
        dbg(2,"X constructed")
        nX = norm(X)
        dbg(2,"about to solve system")
        uv = G\-X
        dbg(2,"solved system")
        dbg(4,"uv = $(uv)")

        u = uv[1:dims]
        v = uv[dims+1:end]

        # constuct Tangent vectors
        nn = getBase(n)
        Λpp = getBase(Λ(p))

        U = zeroTVector(M,p)
        for (jj,j) ∈ enumerate(CartesianIndices(d))
            Θ,κ = tangentONB(M.manifold,p[j],m[j])
            ii = (jj-1)*Mdims
            uu = u[ii+1:ii+Mdims]
            U[j] = sum(uu.*Θ)
        end
        VT = zeroTVector(N.manifold,nn)
        for (jj,j) ∈ enumerate(CartesianIndices(dd))
            Ξ,κ = tangentONB(N.manifold.manifold,nn[j],Λpp[j])
            ii = (jj-1)*Ndims
            vv = v[ii+1:ii+Ndims]
            VT[j] = sum(vv.*Ξ)
        end
        V = TBTVector(getBase(ξ),VT)

        p = exp(M,p,U)
        ξ = ξ + V
        pξ = tuple(p,ξ)
        dbg(4,"p = $(p)")
        dbg(4,"ξ = $(ξ)")

        # checkout new result
        costn = getcost(pr,p)
        change = abs(costn-cost)
        cost = costn
        X = construct_lX(pr,pξ)
        nX = norm(X)
        relerror = nX/nX0

        dbg(1,"# $(k) | Change: $(change) | RelError: $(relerror) | F(p) = $(cost)")
        if save_iterates
            rn = (k, cost, change, p, ξ, relerror)
        else
            rn = (k, cost, change, relerror)
        end
        push!(r,rn)

        k = k+1
    end

    if return_dual
        return p,ξ,r
    else
        return p,r
    end

end

function construct_lX(pr::P,pξ::Tuple{<:MPoint,<:TVector}) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

    dims = manifoldDimension(M)
    d = M.powerSize
    ldims = length(d)
    Mdims = manifoldDimension(M.manifold)
    NN = N.manifold
    dd = NN.powerSize
    dualdims = manifoldDimension(NN)
    Ndims = manifoldDimension(N.manifold.manifold)

    Λ = pr.operator
    DΛ = pr.operatorD
    DΛadj = pr.operatorAdjD
    Λm = Λ(m)

    σ = pr.σ
    τ = pr.τ

    proxf = pr.primalProx
    proxg = pr.dualProx

    p,ξ = pξ

    ξ1 = log(M,m,p)
    ξ2 = DΛ(m,ξ1)
    ξ3 = parallelTransport(N,Λm,n,ξ2)
    ξ4 = ξ + τ*ξ3
    ξn = proxg(ξ4,τ)

    η1 = parallelTransport(N,n,Λm,ξ) # ∈ TΛ(m)N
    η2 = - σ * DΛadj(m,η1) # officially ∈ T*mM, but embedded in TmM
    η3 = parallelTransport(M,m,p,η2) # ∈ TpM
    p1 = exp(M,p,η3)
    pn = proxf(p1,σ)

    X1 = - log(M,p,pn)
    X2 = getTangent(ξ) - getTangent(ξn)

    X = Float64[]
    for j ∈ CartesianIndices(d)
        Θ,κ = tangentONB(M.manifold,p[j],m[j])
        for i in 1:Mdims
            Xᵢ = dot(M.manifold,p[j],Θ[i],X1[j])
            push!(X,Xᵢ)
        end
    end

    nn = getBase(n)
    Λpp = getBase(Λ(p))
    for j ∈ CartesianIndices(dd)
        Ξ,κ = tangentONB(N.manifold.manifold,nn[j],Λpp[j])
        for i in 1:Ndims
            Xᵢ = dot(N.manifold.manifold,nn[j],Ξ[i],X2[j])
            push!(X,Xᵢ)
        end
    end

    return X

end

function construct_lDpX1(pr::P,pξ::Tuple{<:MPoint,<:TVector}) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

    d = M.powerSize
    dims = manifoldDimension(M)
    Mdims = manifoldDimension(M.manifold)

    Λ = pr.operator
    DΛadj = pr.operatorAdjD
    Λm = Λ(m)

    σ = pr.σ

    proxf = pr.primalProx
    Dcproxf = pr.primalProxCdiff

    p,ξ = pξ

    η1 = parallelTransport(N,n,Λm,ξ) # ∈ TΛ(m)N
    qξ = - σ * DΛadj(m,η1)
    qₚ = geodesic(M,m,p,1/2)
    qb = exp(M,m,qξ)
    q₅ = 2* log(M,qb,qₚ)
    q₄ = exp(M,qb,q₅)
    q₃ = -log(M,p,q₄)
    q₂ = exp(M,p,q₃)
    q₁ = proxf(q₂,σ)

    DpX1 = spzeros(dims,dims)
    for (jj,j) ∈ enumerate(CartesianIndices(d))
        Θ,κ = tangentONB(M.manifold,p[j],m[j])
        ii = (jj-1)*Mdims
        for k in 1:Mdims
            Θⱼ = Θ[k]

            Gⱼ = DyGeo(M.manifold,m[j],p[j],1/2,Θⱼ)
            Fⱼ = 2*DyLog(M.manifold,qb[j],qₚ[j],Gⱼ)
            Eⱼ = DξExp(M.manifold,qb[j],q₅[j],Fⱼ)
            D₂ⱼ = - DyLog(M.manifold,p[j],q₄[j],Eⱼ)
            D₁ⱼ = - DxLog(M.manifold,p[j],q₄[j],Θⱼ)
            Dⱼ = D₁ⱼ + D₂ⱼ
            C₂ⱼ = DξExp(M.manifold,p[j],q₃[j],Dⱼ)

            C₁ⱼ = DxExp(M.manifold,p[j],q₃[j],Θⱼ)
            Cⱼ = C₁ⱼ + C₂ⱼ
            Bⱼ = Dcproxf(j,q₂[j],σ,Cⱼ)
            A₂ⱼ = -DyLog(M.manifold,p[j],q₁[j],Bⱼ)
            A₁ⱼ = -DxLog(M.manifold,p[j],q₁[j],Θⱼ)
            Aⱼ =  A₁ⱼ + A₂ⱼ

            for i in 1:Mdims
                Θᵢ = Θ[i]
                DpX1[ii+i,ii+k] += dot(M.manifold,p[j],Θᵢ,Aⱼ)
            end

        end
    end

    return DpX1

end

function construct_lDξX1(pr::P,pξ::Tuple{<:MPoint,<:TVector}) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

    dims = manifoldDimension(M)
    d = M.powerSize
    ldims = length(d)
    Mdims = manifoldDimension(M.manifold)
        NN = N.manifold
        dualdims = manifoldDimension(NN)
        dd = NN.powerSize
        Ndims = manifoldDimension(N.manifold.manifold)

    Λ = pr.operator
    DΛadj = pr.operatorAdjD
    Λm = Λ(m)
    Λmm = getBase(Λm)

    σ = pr.σ

    proxf = pr.primalProx
    Dcproxf = pr.primalProxCdiff

    p,ξ = pξ

    nn = getBase(n)
    Λpp = getBase(Λ(p))

    η1 = parallelTransport(N,n,Λm,ξ) # ∈ TΛ(m)N
    η2 = - σ * DΛadj(m,η1)
    q₃ = parallelTransport(M,m,p,η2)
    q₂ = exp(M,p,q₃)
    q₁ = proxf(q₂,σ)

    R = CartesianIndices(d)
    maxInd = [last(R).I...]

    D = vcat([1],[prod(d[1:i]) for i in 1:length(d)-1])

    DξX1 = spzeros(dims,dualdims)
    for (jj,j) ∈ enumerate(CartesianIndices(dd))
        Ξ,κ = tangentONB(N.manifold.manifold,nn[j],Λpp[j])
        ii = (jj-1)*Ndims
        jjp₁ = ((jj-1)%prod(d))+1
        ii₁ = (jjp₁-1)*Mdims
        jp₁ = R[jjp₁]
        K = floor((jj-1)/prod(d))+1

        Jp₁ = [jp₁.I...] # array of index
        Jp₂ = Jp₁ .+ 1 .* (1:ldims .== K) #i + e_k is j
        if all( Jp₂ .<= maxInd )
            for k in 1:Ndims
                Ξⱼ = Ξ[k]

                jp₂ = CartesianIndex(Jp₂...) # neigbbor index as Cartesian Index
                jjp₂ = sum(D.*(Jp₂ .- 1 .*(1:ldims .!=1)))
                ii₂ = (jjp₂-1)*Mdims
                h1 = parallelTransport(N.manifold.manifold,nn[j],Λmm[j],Ξⱼ) # ∈ TΛ(m)N

                h2₁ = -σ*AdjDxLog(M.manifold,m[jp₁],m[jp₂],h1)
                Hⱼ₁ = parallelTransport(M.manifold,m[jp₁],p[jp₁],h2₁)
                C₂ⱼ₁ = DξExp(M.manifold,p[jp₁],q₃[jp₁],Hⱼ₁)
                Bⱼ₁ = Dcproxf(jp₁,q₂[jp₁],σ,C₂ⱼ₁)
                A₂ⱼ₁ = -DyLog(M.manifold,p[jp₁],q₁[jp₁],Bⱼ₁)

                h2₂ = -σ*AdjDyLog(M.manifold,m[jp₁],m[jp₂],h1)
                Hⱼ₂ = parallelTransport(M.manifold,m[jp₂],p[jp₂],h2₂)
                C₂ⱼ₂ = DξExp(M.manifold,p[jp₂],q₃[jp₂],Hⱼ₂)
                Bⱼ₂ = Dcproxf(jp₂,q₂[jp₂],σ,C₂ⱼ₂)
                A₂ⱼ₂ = -DyLog(M.manifold,p[jp₂],q₁[jp₂],Bⱼ₂)

                Θ₁,κ = tangentONB(M.manifold,p[jp₁],m[jp₁])
                Θ₂,κ = tangentONB(M.manifold,p[jp₂],m[jp₂])
                for i in 1:Mdims
                    Θ₁ᵢ = Θ₁[i]
                    DξX1[ii₁+i,ii+k] += dot(M.manifold,p[jp₁],Θ₁ᵢ,A₂ⱼ₁)
                    Θ₂ᵢ = Θ₂[i]
                    DξX1[ii₂+i,ii+k] += dot(M.manifold,p[jp₂],Θ₂ᵢ,A₂ⱼ₂)
                end
            end
        end
    end

    return DξX1

end

function construct_lDpX2(pr::P,pξ::Tuple{<:MPoint,<:TVector}) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

    dims = manifoldDimension(M)
    Mdims = manifoldDimension(M.manifold)
    d = M.powerSize
    ldims = length(d)
    MM = Power(M.manifold,(ldims,))

        NN = N.manifold
        dualdims = manifoldDimension(NN)
        Ndims = manifoldDimension(NN.manifold)

    Λ = pr.operator
    DΛ = pr.operatorD
    Λm = Λ(m)

    τ = pr.τ

    proxg = pr.dualProx
    Dcproxg = pr.dualProxCdiff

    p,ξ = pξ

    nn = getBase(n)
    Λpp = getBase(Λ(p))

    ηₚ = log(M,m,p)
    η₃ = DΛ(m,ηₚ)
    η₂ = parallelTransport(N,Λm,n,η₃)
    η₁ = getTangent(ξ + τ*η₂)

    R = CartesianIndices(d)
    maxInd = [last(R).I...]
    minInd = [first(R).I...]
    D = vcat([1],[prod(d[1:i]) for i in 1:length(d)-1])

    DpX2 = spzeros(dualdims,dims)
    for (jj,j) ∈ enumerate(CartesianIndices(d))
        Θ,κ = tangentONB(M.manifold,p[j],m[j])
        JJ = (jj-1)*Mdims
        J = [j.I...] # array of index
        for K in 1:ldims
            J₁ = J .+ 1 .* (1:ldims .== K) #i + e_k is j
            if all( J₁ .<= maxInd )
                for k in 1:Mdims
                    Θⱼ = Θ[k]

                    j₁ = CartesianIndex(J₁...)
                    ii = Ndims*((K-1)*prod(d) + sum(D.*(J .- 1 .*(1:ldims .!=1)))-1)

                    Mⱼ = DyLog(M.manifold,m[j],p[j],Θⱼ)
                    Lⱼ = DxLog(M.manifold,m[j],m[j₁],Mⱼ)
                    Kⱼ = τ * parallelTransport(N.manifold.manifold,getBase(Λm)[j,K],getBase(n)[j,K],Lⱼ)
                    Jⱼ = Dcproxg(j,K,η₁[j,:],τ,Kⱼ)

                    Ξ,κ = tangentONB(N.manifold.manifold,nn[j,K],Λpp[j,K])
                    for i in 1:Ndims
                        Ξᵢ = Ξ[i]
                        if Jⱼ isa TVector
                            DpX2[ii+i,JJ+k] += -dot(N.manifold.manifold,nn[j,K],Ξᵢ,Jⱼ)
                        else
                            for ll in 1:ldims
                                KK,JJⱼ = Jⱼ[ll]
                                II = Ndims*((KK-1)*prod(d) + sum(D.*(J .- 1 .*(1:ldims .!=1)))-1)
                                DpX2[II+i,JJ+k] += -dot(N.manifold.manifold,nn[j,K],Ξᵢ,JJⱼ)
                            end
                        end
                    end
                end
            end

            J₂ = J .- 1 .* (1:ldims .== K) #i - e_k is j
            if all( J₂ .>= minInd )
                j₂ = CartesianIndex(J₂...)
                for k in 1:Mdims
                    Θⱼ = Θ[k]
                    ii = Ndims*((K-1)*prod(d) + sum(D.*(J₂ .- 1 .*(1:ldims .!=1)))-1)

                    Mⱼ = DyLog(M.manifold,m[j],p[j],Θⱼ)
                    Lⱼ = DyLog(M.manifold,m[j₂],m[j],Mⱼ)
                    Kⱼ = τ * parallelTransport(N.manifold.manifold,getBase(Λm)[j₂,K],getBase(n)[j₂,K],Lⱼ)
                    Jⱼ = Dcproxg(j₂,K,η₁[j₂,:],τ,Kⱼ)

                    Ξ,κ = tangentONB(N.manifold.manifold,nn[j₂,K],Λpp[j₂,K])
                    for i in 1:Ndims
                        Ξᵢ = Ξ[i]
                        if Jⱼ isa TVector
                            DpX2[ii+i,JJ+k] += -dot(N.manifold.manifold,nn[j₂,K],Ξᵢ,Jⱼ)
                        else
                            for ll in 1:ldims
                                KK,JJⱼ = Jⱼ[ll]
                                II = Ndims*((KK-1)*prod(d) + sum(D.*(J₂ .- 1 .*(1:ldims .!=1)))-1)
                                DpX2[II+i,JJ+k] += -dot(N.manifold.manifold,nn[j₂,K],Ξᵢ,JJⱼ)
                            end
                        end
                    end
                end
            end
        end
    end

    return DpX2

end

function construct_lDξX2(pr::P,pξ::Tuple{<:MPoint,<:TVector}) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

    dims = manifoldDimension(M)
    Mdims = manifoldDimension(M.manifold)
    d = M.powerSize
    ldims = length(d)
        NN = N.manifold
        dd = NN.powerSize
        dualdims = manifoldDimension(NN)
        Ndims = manifoldDimension(NN.manifold)


    Λ = pr.operator
    DΛ = pr.operatorD
    Λm = Λ(m)

    τ = pr.τ

    proxg = pr.dualProx
    Dcproxg = pr.dualProxCdiff

    p,ξ = pξ

    ηₚ = log(M,m,p)
    η₃ = DΛ(m,ηₚ)
    η₂ = parallelTransport(N,Λm,n,η₃)
    η₁ = getTangent(ξ + τ*η₂)

    MM = Power(M.manifold,(ldims,))
    nn = getBase(n)
    Λpp = getBase(Λ(p))

    D = vcat([1],[prod(d[1:i]) for i in 1:length(d)-1])

    DξX2 = spzeros(dualdims,dualdims)
    for (jj,j) ∈ enumerate(CartesianIndices(dd))
        Ξ,κ = tangentONB(N.manifold.manifold,nn[j],Λpp[j])
        JJ = (jj-1)*Mdims
        if length(j)>1
            j_ = CartesianIndex(j.I[1:end-1])
            K = j.I[end]
        else
            j_ = j
            K = 1
        end
        J = [j_.I...]
        for k in 1:Ndims
            Ξⱼ = Ξ[k]
            Jⱼ = Dcproxg(j_,K,η₁[j_,:],τ,Ξⱼ)
            if Jⱼ isa TVector
                Iⱼ = Ξⱼ - Jⱼ
            else
                Iⱼ = Tuple{Int,TVector}[]
                for ll in 1:ldims
                    if Jⱼ[ll][1] != K
                        push!(Iⱼ,tuple(Jⱼ[ll][1],- Jⱼ[ll][2]))
                    else
                        push!(Iⱼ,tuple(Jⱼ[ll][1],Ξⱼ - Jⱼ[ll][2]))
                    end
                end
            end

            for i in 1:Ndims
                Ξᵢ = Ξ[i]
                if Iⱼ isa TVector
                    DξX2[JJ+i,JJ+k] += dot(N.manifold.manifold,nn[j],Ξᵢ,Iⱼ)
                else
                    for ll in 1:ldims
                        KK,IIⱼ = Iⱼ[ll]
                        II = Ndims*((KK-1)*prod(d) + sum(D.*(J .- 1 .*(1:ldims .!=1)))-1)
                        DξX2[II+i,JJ+k] += dot(N.manifold.manifold,nn[j],Ξᵢ,IIⱼ)
                    end
                end
            end

        end

    end

    return DξX2

end
