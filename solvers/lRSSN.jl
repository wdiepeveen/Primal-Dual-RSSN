include(joinpath("..","jltools","ExpTools.jl"))
include(joinpath("..","jltools","Dbg.jl"))
include(joinpath("..","plans","problem.jl"))

using Manopt, SparseArrays, LinearAlgebra

struct lrssnProblem <: tProblem
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

    lrssnProblem(M::Power,m::PowPoint,
    N::Union{Power,TangentBundle},n::Union{PowPoint,TBPoint},
    cF::Function,pP::Function,dP::Function,
    Λ::Function,DΛ::Function,DΛadj::Function,
    pPCd::Function,dPCd::Function)  =
    new(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,pPCd,dPCd,1/2,1/2)

    lrssnProblem(M::Power,m::PowPoint,
    N::Union{Power,TangentBundle},n::Union{PowPoint,TBPoint},
    cF::Function,pP::Function,dP::Function,
    Λ::Function,DΛ::Function,DΛadj::Function,
    pPCd::Function,dPCd::Function,σ::Float64,τ::Float64)  =
    new(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,pPCd,dPCd,σ,τ)
end

function lRSSN(pr::P,x0::PowPoint;y0::Union{TVector,Missing}=missing,
    tol::Float64=1e-12,return_dual::Bool=false,max_it::Int=20,
    save_iterates::Bool=false) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n
    d = manifoldDimension(M)

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
    k = 1
    dbg(1,"# $(k) | Change: $(change) | RelError: 1. | F(p) = $(cost)")

    Θ,κ = tangentONB(M,p,m) # must have base p
    dbg(2,"primal ONB constructed")
    if N isa TangentBundle
        Ξ,κ = tangentONB(N.manifold,getBase(n),getBase(Λ(p))) # must have base n
    else
        Ξ,κ = tangentONB(N,n,Λ(p)) # must have base n
    end

    X = construct_X(pr,pξ,Θ,Ξ)
    nX0 = norm(X)
    relerror = 1.

    if dbglevel_atleast(3)
        maxit = 3
    else
        maxit = max_it
    end

    while k<=maxit && relerror >tol
        dbg(4,"k = $(k)")
        dbg(4,"p = $(p)")
        dbg(4,"ξ = $(ξ)")

        Θ,κ = tangentONB(M,p,m) # must have base p
        dbg(3,"primal ONB constructed")
        if N isa TangentBundle
            Ξ,κ = tangentONB(N.manifold,getBase(n),getBase(Λ(p))) # must have base n
        else
            Ξ,κ = tangentONB(N,n,Λ(p)) # must have base n
        end
        dbg(3,"dual ONB constructed")

        DpX1 = construct_DpX1(pr,pξ,Θ)
        dbg(3,"DpX1 constructed")
        DξX1 = construct_DξX1(pr,pξ,Θ,Ξ)
        dbg(3,"DξX1 constructed")
        DpX2 = construct_DpX2(pr,pξ,Θ,Ξ)
        dbg(3,"DpX2 constructed")
        DξX2 = construct_DξX2(pr,pξ,Ξ)
        dbg(3,"DξX2 constructed")


        G = [DpX1 DξX1;
            DpX2 DξX2]
        dbg(3,"G constructed")
        dbg(2,"G = $(Array(G))")
        dbg(2,"cond(G) = $(cond(Array(G)))")

        X = construct_X(pr,pξ,Θ,Ξ)
        dbg(2,"X = $(X)")
        dbg(3,"X constructed")
        nX = norm(X)
        dbg(3,"about to solve system")
        uv = G\-X
        dbg(3,"solved system")
        dbg(4,"uv = $(uv)")

        u = uv[1:d]
        v = uv[d+1:end]

        # constuct Tangent vectors
        U = sum(u.*Θ)

        if N isa TangentBundle
            VT = sum(v.*Ξ)
            V = TBTVector(getBase(ξ),VT)
        else
            V = sum(v.*Ξ)
        end

        p = exp(M,p,U)
        ξ = ξ + V
        pξ = tuple(p,ξ)
        dbg(4,"p = $(p)")
        dbg(4,"ξ = $(ξ)")

        # checkout new result
        costn = getcost(pr,p)
        change = abs(costn-cost)
        cost = costn
        X = construct_X(pr,pξ,Θ,Ξ)
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

function construct_X(pr::P,pξ::Tuple{<:MPoint,<:TVector},
    Θ,Ξ) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

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
    ξn = proxg(ξ4,τ) # right format for Tangent Bundle Tangent vector

    η1 = parallelTransport(N,n,Λm,ξ) # ∈ TΛ(m)N
    η2 = - σ * DΛadj(m,η1) # officially ∈ T*mM, but embedded in TmM
    η3 = parallelTransport(M,m,p,η2) # ∈ TpM
    p1 = exp(M,p,η3)
    pn = proxf(p1,σ)

    X1 = - log(M,p,pn)
    if N isa TangentBundle
        X2 = getTangent(ξ) - getTangent(ξn)
    else
        X2 = ξ - ξn
    end

    X = Float64[]
    for Θᵢ ∈ Θ
        Xᵢ = dot(M,p,Θᵢ,X1)
        push!(X,Xᵢ)
    end
    for Ξᵢ ∈ Ξ
        if N isa TangentBundle
            Xᵢ = dot(N.manifold,getBase(n),Ξᵢ,X2)
        else
            Xᵢ = dot(N,n,Ξᵢ,X2)
        end
        push!(X,Xᵢ)
    end

    return X

end

function construct_DpX1(pr::P,pξ::Tuple{<:MPoint,<:TVector},
    Θ) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

    dims = manifoldDimension(M)

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
    for (j,Θⱼ) ∈ enumerate(Θ)

        Gⱼ = DyGeo(M,m,p,1/2,Θⱼ)
        Fⱼ = 2*DyLog(M,qb,qₚ,Gⱼ)
        Eⱼ = DξExp(M,qb,q₅,Fⱼ)
        D₂ⱼ = - DyLog(M,p,q₄,Eⱼ)
        D₁ⱼ = - DxLog(M,p,q₄,Θⱼ)
        Dⱼ = D₁ⱼ + D₂ⱼ
        C₂ⱼ = DξExp(M,p,q₃,Dⱼ)

        C₁ⱼ = DxExp(M,p,q₃,Θⱼ)
        Cⱼ = C₁ⱼ + C₂ⱼ
        Bⱼ = Dcproxf(q₂,σ,Cⱼ)
        A₂ⱼ = -DyLog(M,p,q₁,Bⱼ)
        A₁ⱼ = -DxLog(M,p,q₁,Θⱼ)
        Aⱼ =  A₁ⱼ + A₂ⱼ

        DpX1ⱼ = spzeros(dims)
        for (i,Θᵢ) ∈ enumerate(Θ)
            DpX1ᵢⱼ = dot(M,p,Θᵢ,Aⱼ)
            if DpX1ᵢⱼ !=0.
                DpX1ⱼ[i] = DpX1ᵢⱼ
            end
        end

        dropzeros!(DpX1ⱼ)
        DpX1[:,j] = DpX1ⱼ

    end

    return DpX1

end

function construct_DξX1(pr::P,pξ::Tuple{<:MPoint,<:TVector},
    Θ,Ξ) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

    dims = manifoldDimension(M)
    if N isa TangentBundle
        NN = N.manifold
        dualdims = manifoldDimension(NN)
    else
        dualdims = manifoldDimension(N)
    end

    Λ = pr.operator
    DΛadj = pr.operatorAdjD
    Λm = Λ(m)

    σ = pr.σ

    proxf = pr.primalProx
    Dcproxf = pr.primalProxCdiff

    p,ξ = pξ

    η1 = parallelTransport(N,n,Λm,ξ) # ∈ TΛ(m)N
    η2 = - σ * DΛadj(m,η1)
    q₃ = parallelTransport(M,m,p,η2)
    q₂ = exp(M,p,q₃)
    q₁ = proxf(q₂,σ)

    DξX1 = spzeros(dims,dualdims)

    for (j,Ξⱼ) ∈ enumerate(Ξ)

        if N isa TangentBundle
            Ξⱼ = TBTVector(getBase(ξ),Ξⱼ)
        end

        h1 = parallelTransport(N,n,Λm,Ξⱼ) # ∈ TΛ(m)N
        h2 = - σ * DΛadj(m,h1) # officially ∈ T*mM, but embedded in TmM
        Hⱼ = parallelTransport(M,m,p,h2)
        C₂ⱼ = DξExp(M,p,q₃,Hⱼ)
        Bⱼ = Dcproxf(q₂,σ,C₂ⱼ)
        A₂ⱼ = -DyLog(M,p,q₁,Bⱼ)

        DξX1ⱼ = spzeros(dims)
        for (i,Θᵢ) ∈ enumerate(Θ)
            DξX1ᵢⱼ = dot(M,p,Θᵢ,A₂ⱼ)
            if DξX1ᵢⱼ !=0.
                DξX1ⱼ[i] = DξX1ᵢⱼ
            end
        end

        dropzeros!(DξX1ⱼ)
        DξX1[:,j] = DξX1ⱼ

    end

    return DξX1

end

function construct_DpX2(pr::P,pξ::Tuple{<:MPoint,<:TVector},
    Θ,Ξ) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

    dims = manifoldDimension(M)
    if N isa TangentBundle
        NN = N.manifold
        dualdims = manifoldDimension(NN)
    else
        dualdims = manifoldDimension(N)
    end

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
    η₁ = ξ + τ*η₂

    DpX2 = spzeros(dualdims,dims)
    for (j,Θⱼ) ∈ enumerate(Θ)

        Mⱼ = DyLog(M,m,p,Θⱼ)
        # Lⱼ = Mⱼ
        Kⱼ = τ*parallelTransport(N,Λm,n,DΛ(m,Mⱼ))
        Jⱼ = Dcproxg(η₁,τ,Kⱼ)

        DpX2ⱼ = spzeros(dualdims)
        for (i,Ξᵢ) ∈ enumerate(Ξ)
            if N isa TangentBundle
                DpX2ᵢⱼ = -dot(N.manifold,getBase(n),Ξᵢ,Jⱼ)
                if DpX2ᵢⱼ !=0.
                    DpX2ⱼ[i] += DpX2ᵢⱼ
                end
            else
                DpX2ᵢⱼ = -dot(N,n,Ξᵢ,Jⱼ)
                if DpX2ᵢⱼ !=0.
                    DpX2ⱼ[i] += DpX2ᵢⱼ
                end
            end
        end

        dropzeros!(DpX2ⱼ)
        DpX2[:,j] = DpX2ⱼ

    end

    return DpX2

end

function construct_DξX2(pr::P,pξ::Tuple{<:MPoint,<:TVector},
    Ξ) where {P <: tProblem}

    M = pr.M
    m = pr.m
    N = pr.N
    n = pr.n

    if N isa TangentBundle
        NN = N.manifold
        dualdims = manifoldDimension(NN)
    else
        dualdims = manifoldDimension(N)
    end

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
    η₁ = ξ + τ*η₂

    DξX2 = spzeros(dualdims,dualdims)
    for (j,Ξⱼ) ∈ enumerate(Ξ)

        if N isa TangentBundle
            Jⱼ = Dcproxg(η₁,τ,TBTVector(getBase(ξ),Ξⱼ))
        else
            Jⱼ = Dcproxg(η₁,τ,Ξⱼ)
        end


        Iⱼ = Ξⱼ - Jⱼ

        DξX2ⱼ = spzeros(dualdims)
        for (i,Ξᵢ) ∈ enumerate(Ξ)
            if N isa TangentBundle
                DξX2ᵢⱼ = dot(N.manifold,getBase(n),Ξᵢ,Iⱼ)
                if DξX2ᵢⱼ !=0.
                    DξX2ⱼ[i] += DξX2ᵢⱼ
                end
            else
                DξX2ᵢⱼ = dot(N,n,Ξᵢ,Iⱼ)
                if DξX2ᵢⱼ !=0.
                    DξX2ⱼ[i] += DξX2ᵢⱼ
                end
            end
        end

        dropzeros!(DξX2ⱼ)
        DξX2[:,j] = DξX2ⱼ

    end

    return DξX2

end
