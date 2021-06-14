include(joinpath("..","jltools","ExpTools.jl"))
include(joinpath("..","jltools","Dbg.jl"))
include(joinpath("..","functions","ManTVproxes.jl"))
include(joinpath("..","functions","ManTVdualbasepoint.jl"))
include(joinpath("..","plans","problem.jl"))

using Manopt, Plots

pyplot()

function main()

        exp_begin(aux_folders = ["experiments","demos","functions","jltools","plans","solvers"],
        prefix = "demo_lRCPA_S1_")
        dbglevel(2)

        dbg(0,"Results folder is " * exp_prefix())

        # Settings
        n = 1 # d = 2n
        α = 5.

        dataColor = RGBA{Float64}(colorant"#BBBBBB")
        s2dColor = RGBA{Float64}(colorant"#EE7733") # data Color: Tol Vibrant Orange
        s1Color = RGBA{Float64}(colorant"#0077BB") # control point data color: Tol Virbant Blue
        nColor = RGBA{Float64}(colorant"#33BBEE") # tangent vector: Tol Vibrant Teal
        mColor = RGBA{Float64}(colorant"#EE2211") # base point: Tol Vibrant Red
        lColor = RGBA{Float64}(colorant"#BB5511") # slice points: Tol Vibrant Red


        MM = Circle()
        M = Power(MM, (2*n,))

        θ = 2π/3
        f = vcat(fill(S1Point(θ),n),fill(S1Point(π-θ),n))
        data = PowPoint(f)
        t = range(1,2*n,length=2*n)

        θstart = -π/3
        start = PowPoint(fill(S1Point(θstart),2*n))
        if n==1
                δ = min(1/2, α/distance(MM,f[1],f[2]))
                exact = PowPoint([geodesic(MM,f[1],f[2],δ), geodesic(MM,f[2],f[1],δ)])
        end
        m = PowPoint(fill(S1Point(0.),2*n))

        # Optimize
        fR,r = lRCPA_TV(M,m,data,α,start)

        # Results

        # prepare TV contour plot
        if n ==1
                steps = 50
                dp = 2*π/steps
                prange = -π:dp:π
                p1 = [S1Point(i) for i in -π:dp:π]
                p2 = p1

                P1 = repeat(reshape(p1, 1, :), steps, 1)
                P2 = repeat(p2, 1, steps)
                ftv = (p,q) -> costL2TV(M,data,α,PowPoint([p,q]))


                Z = map(ftv,P1,P2)

                cont = contour(prange,prange,Z, levels=50)

                # plot of the optimization landscape along a geodesic
                lsteps = 100
                dl = 2*π/lsteps
                lrange = -π:dl:π

                X = log(M,fR,data)
                g = norm(MM,getValue(fR)[1],getValue(X)[1])
                Xn = X/g
                l1 = [getValue(exp(M,fR,l*Xn))[1] for l in lrange]
                l2 = [getValue(exp(M,fR,l*Xn))[2] for l in lrange]
                slicecost = map(ftv,l1,l2)
                tvslice = plot(lrange,slicecost, color = lColor, lab = "cost")
                xticks!([-π,-π/2,0,π/2,π], ["-π", "- π/2", "0", "π/2", "π"])
                xlabel!("Δp1")

                scatter!(tvslice,[0.],[costL2TV(M,data,α,fR)],
                # scatter!(tvslice,[0.],[ftv(getValue(fR)[1],getValue(fR)[2])],
                markersize=6, markercolor = nColor, markerstrokecolor=nColor,
                lab="reconstruction")
                scatter!(tvslice,[g],[costL2TV(M,data,α,exp(M,fR,X))],
                markersize=6, markercolor = dataColor, markerstrokecolor=dataColor,
                lab="original")
                # plot the slice on the contourplot

                for i in 1:lsteps
                        scatter!(cont,[getValue(l1[i])],
                        [getValue(l2[i])],
                        markersize=1, markercolor = lColor, markerstrokecolor=lColor,
                        lab="")
                end
        end

        # make other plots
        scene = scatter(t,getValue.(f),
        markersize=4, markercolor = dataColor, markerstrokecolor=dataColor,
        lab="original")

        # plot exact solution for comparison
        if n ==1
                scatter!(scene, t,getValue.(getValue(exact)),
                markersize=4, markercolor = s2dColor, markerstrokecolor=s2dColor,
                lab="exact solution")
                scatter!(cont,[getValue.(getValue(data))[1]],[getValue.(getValue(data))[2]],
                markersize=4, markercolor = dataColor, markerstrokecolor=dataColor,
                lab="original")
                scatter!(cont,[getValue.(getValue(m))[1]],[getValue.(getValue(m))[2]],
                markersize=4, markercolor = mColor, markerstrokecolor=mColor,
                lab="m")
                scatter!(cont,[getValue.(getValue(exact))[1]],[getValue.(getValue(exact))[2]],
                markersize=4, markercolor = s2dColor, markerstrokecolor=s2dColor,
                lab="exact solution")
        end
        # add iterates to see progress
        scatter!(scene,t,getValue.(getValue(start)),
        markersize=2, markercolor=s1Color , markerstrokecolor=s1Color,
        lab= "")
        if n==1
                scatter!(cont,[getValue.(getValue(start))[1]],[getValue.(getValue(start))[2]],
                markersize=4, markercolor = s1Color, markerstrokecolor=s1Color,
                lab="")
        end
        for i in 1:min(size(r)[1]-1,50)
                scatter!(scene,t,getValue.(getValue(r[i][4])),
                markersize=2, markercolor=s1Color , markerstrokecolor=s1Color,
                lab= "")
                if n ==1
                        scatter!(cont,[getValue.(getValue(r[i][4]))[1]],
                        [getValue.(getValue(r[i][4]))[2]],
                        markersize=2, markercolor = s1Color, markerstrokecolor=s1Color,
                        lab="")
                end
        end
        scatter!(scene,t,getValue.(getValue(fR)),
        markersize=2, markercolor = nColor, markerstrokecolor=nColor,
        lab="reconstruction")
        yticks!([-π,-π/2,0,π/2,π], ["-π", "- π/2", "0", "π/2", "π"])
        if n ==1
                scatter!(cont,[getValue.(getValue(fR))[1]],[getValue.(getValue(fR))[2]],
                markersize=2, markercolor = nColor, markerstrokecolor=nColor,
                lab="reconstruction")
                yticks!([-π,-π/2,0,π/2,π], ["-π", "- π/2", "0", "π/2", "π"])
                xticks!([-π,-π/2,0,π/2,π], ["-π", "- π/2", "0", "π/2", "π"])
        end

        ϵ = zeros(size(r)[1],2)
        for i in 1:size(r)[1]
                ϵ[i,1] = r[i][1]
                ϵ[i,2] = r[i][2]
        end
        errorfig = plot(ϵ[:,1], ϵ[:,2], yaxis=:log, lab="Cost")



        exp_savefig(scene, "result", true)

        exp_savefig(errorfig, "cost", true)

        if n==1
                exp_savefig(cont, "optimization_landscape", true)

                exp_savefig(tvslice, "optimization_slice", true)
        end
        exp_save("r_Iteration-Cost-Change-Iterate", "r", r)

        exp_end()

        return r

end

struct lrcpaProblem <: tProblem
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
        σ::Float64
        τ::Float64
        θ::Float64
        γ::Float64
        lrcpaProblem(M::Power,m::PowPoint,
        N::Union{Power,TangentBundle},n::Union{PowPoint,TBPoint},
        cF::Function,pP::Function,dP::Function,
        Λ::Function,DΛ::Function,DΛadj::Function)  = new(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,1/2,1/2,1.,0.)
        lrcpaProblem(M::Power,m::PowPoint,
        N::Union{Power,TangentBundle},n::Union{PowPoint,TBPoint},
        cF::Function,pP::Function,dP::Function,
        Λ::Function,DΛ::Function,DΛadj::Function,
        σ::Float64,τ::Float64)  = new(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,σ,τ,1.,0.)
        lrcpaProblem(M::Power,m::PowPoint,
        N::Union{Power,TangentBundle},n::Union{PowPoint,TBPoint},
        cF::Function,pP::Function,dP::Function,
        Λ::Function,DΛ::Function,DΛadj::Function,
        σ::Float64,τ::Float64,θ::Float64,γ::Float64)  =
        new(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,σ,τ,θ,γ)
end

function lRCPA(pr::P,x0::PowPoint;y0::Union{PowTVector,TBTVector,Missing}=missing,
        tol=1e-9,return_dual::Bool=false,max_it::Int=500,
        save_iterates::Bool=false) where {P <: tProblem}

        M = pr.M
        m = pr.m
        N = pr.N
        n = pr.n

        Λ = pr.operator
        DΛ = pr.operatorD
        DΛadj = pr.operatorAdjD
        Λm = Λ(m)

        # Setup algorithm
        p = x0
        if !ismissing(y0)
                ξ = y0
        else
                ξ = zeroTVector(N,n)
        end
        ξb = ξ

        σ = pr.σ
        τ = pr.τ
        θ = pr.θ
        γ = pr.γ

        proxf = pr.primalProx
        proxg = pr.dualProx

        change = 1
        cost = getcost(pr,p)
        if save_iterates
                r = Tuple{Int32,Float64,Float64,PowPoint,TBTVector,Float64}[]
        else
                r = Tuple{Int32,Float64,Float64,Float64}[]
        end
        nX0 = getlRCPAerror(pr,p,ξ)
        relerror = 1.
        dbg(1,"# 0 | Change: - | RelError: 1. | F(p) = $(cost)")
        k = 1
        if dbglevel_atleast(2)
                maxit = 5
        else
                maxit = max_it
        end

        while k<=maxit && relerror >tol
                dbg(4,"k = $(k)")
                dbg(4,"m = $(m)")
                dbg(4,"n = $(n)")
                dbg(4,"p = $(p)")
                dbg(4,"ξb = $(ξb)")
                # TODO parallel transport
                ξ1 = parallelTransport(N,n,Λm,ξb)
                dbg(4,"ξ1 = $(ξ1)")
                # TODO AdjDforwardLogs
                η1 = - σ * DΛadj(m,ξ1)
                dbg(4,"η1 = $(η1)")
                # TODO parallel transport
                η2 = parallelTransport(M,m,p,η1)
                dbg(4,"η2 = $(η2)")
                # TODO exponential map
                p1 = exp(M,p,η2)
                dbg(4,"p1 = $(p1)")
                # TODO prox primal
                p = proxf(p1,σ)
                dbg(4,"p = $(p)")

                # TODO log p new
                η3 = log(M,m,p)
                dbg(4,"η3 = $(η3)")
                # TODO differential forwardLogs
                ξ2 = DΛ(m,η3)
                # ξ2 = TBTVector(forwardLogs(M,m),DforwardLogs(M,m,η3)) # probably add base point here
                dbg(4,"ξ2 = $(ξ2)")
                # TODO parallel transport
                ξ3 = parallelTransport(N,Λm,n,ξ2)
                dbg(4,"ξ3 = $(ξ3)")
                # TODO sum  with former iterate
                ξ4 = ξ + τ*ξ3
                dbg(4,"ξ4 = $(ξ4)")
                # TODO prox dual
                ξ5 = proxg(ξ4,τ)
                dbg(4,"ξ5 = $(ξ5)")

                θ = 1/sqrt(1 + 2*γ*σ)
                σ = σ*θ
                τ = τ/θ

                ξ6 = ξ5 + θ*(ξ5 - ξ)
                dbg(4,"ξ6 = $(ξ6)")

                ξb = ξ6
                dbg(4,"ξb = $(ξb)")

                ξ = ξ5
                dbg(4,"ξ = $(ξ)")

                # checkout new result
                costn = getcost(pr,p)
                change = abs(costn-cost)
                cost = costn
                nX = getlRCPAerror(pr,p,ξ)
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
        if !return_dual
                return p,r
        else
                return p,ξ,r
        end

end

function getlRCPAerror(pr,p,ξ)

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

        nn = getBase(n)
        Λpp = getBase(Λ(p))

        X1 = - log(M,p,pn)
        nX1 = norm(M,p,X1)
        if N isa TangentBundle
                X2 = getTangent(ξ) - getTangent(ξn)
                nX2 = norm(N.manifold,getBase(n),X2)
        else
                X2 = ξ - ξn
                nX2 = norm(N,n,X2)
        end

        nX = sqrt(nX1^2 + nX2^2)

        return nX

end
