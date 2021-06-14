include(joinpath("..","jltools","ExpTools.jl"))
include(joinpath("..","jltools","Dbg.jl"))
include(joinpath("..","plans","problem.jl"))
include("lRSSN-ML.jl")

using Manopt, SparseArrays, LinearAlgebra, Random

struct ilrssnProblem <: tProblem
    M::Power
    m::PowPoint
    N::Union{Power,TangentBundle}
    n::Union{PowPoint,TBPoint}
    costFunction::Function
    primalProx::Function
    dualProx::Function
    operator::Function
    operatorD::Function
    operatorAdjD::Function
    primalProxCdiff::Function
    dualProxCdiff::Function
    σ::Float64
    τ::Float64
    αk::Function

    ilrssnProblem(M::Power,m::PowPoint,
    N::Union{Power,TangentBundle},n::Union{PowPoint,TBPoint},
    cF::Function,pP::Function,dP::Function,
    Λ::Function,DΛ::Function,DΛadj::Function,
    pPCd::Function,dPCd::Function,α::Function)  =
    new(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,pPCd,dPCd,1/2,1/2,α)

    ilrssnProblem(M::Power,m::PowPoint,
    N::Union{Power,TangentBundle},n::Union{PowPoint,TBPoint},
    cF::Function,pP::Function,dP::Function,
    Λ::Function,DΛ::Function,DΛadj::Function,
    pPCd::Function,dPCd::Function,σ::Float64,τ::Float64,α::Function)  =
    new(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,pPCd,dPCd,σ,τ,α)
end

function ilRSSN_ML(pr::P,x0::PowPoint;y0::Union{TVector,Missing}=missing,
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
    α = pr.αk

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
        # TODO add noise
        noisyX = X + α(nX,k)*randn(length(X))
        dbg(2,"about to solve system")
        uv = G\-noisyX
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
