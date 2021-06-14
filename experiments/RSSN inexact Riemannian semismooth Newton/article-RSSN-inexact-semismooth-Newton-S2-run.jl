include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))
include(joinpath("..","..","functions","ManTV.jl"))
include(joinpath("..","..","functions","ManTVdualbasepoint.jl"))
include(joinpath("..","..","functions","ManTVRproxes.jl"))
include(joinpath("..","..","functions","ManTVproxCdifferentials.jl"))
include(joinpath("..","..","solvers","lRCPA.jl"))
include(joinpath("..","..","solvers","ilRSSN-ML.jl"))

using Plots, Manopt, Images, Random, ColorSchemes, CPUTime

pyplot()

function run_RSSN_inexact_S2(α,δ;postfix="")

   if !(postfix=="")
      postfix = "-"*postfix
   end

   # Import data
   Orig = artificialS2Lemniscate(SnPoint(1/sqrt(2)*[1.,0.,1.]))
   orig = PowPoint(Orig)

   # Settings
   d = size(Orig)

   MM = Sphere(2)
   M = Power(MM, d)
   m = PowPoint(fill(SnPoint(1/sqrt(2)*[1.,0.,1.]),d))

   N = TangentBundle(Power(MM,d))
   n = TBPoint(ManTVdualbasepoint(m),forwardLogs(M,m))

   Random.seed!(31)

   Data = addNoise.(Ref(MM),Orig,Ref(δ))
   data = PowPoint(Data)

   anisotropic = true

   # lRCPA step

   cF = p -> ManTV(M,data,α,p,anisotropic)
   pP = (p,σ) -> proxDistance(M,σ/α,data,p)
   dP = (ξ,τ) -> ManTVproxDual(N,n,ξ,anisotropic)
   Λ = p -> TBPoint(ManTVdualbasepoint(p),forwardLogs(M,p))
   DΛ = (p,η) -> TBTVector(zeroTVector(N.manifold,getBase(n)),DforwardLogs(M,p,η))
   DΛadj = (p,η) -> AdjDforwardLogs(M,p,getTangent(η))
   σ = τ = 0.35
   θ = 1.
   γ = 0.2
   pr = lrcpaProblem(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,σ,τ,θ,γ)

   # Optimize
   start = data
   CPUtic()
   p,ξ,rprior = lRCPA(pr,start,return_dual=true,tol=1e-1)
   t0 = CPUtoc()
   relerror = rprior[end][4]

   # prepare lRSSN
   cF = p -> ManTV(M,data,α,p)
   pP = (p,σ) -> proxDistance(M,σ/α,data,p)
   dP = (ξ,τ) -> ManTVproxDual(N,n,ξ)
   Λ = p -> TBPoint(ManTVdualbasepoint(p),forwardLogs(M,p))
   DΛ = (p,η) -> TBTVector(zeroTVector(N.manifold,getBase(n)),
   DforwardLogs(M,p,η))
   DΛadj = (p,η) -> AdjDforwardLogs(M,p,getTangent(η))
   pPCd = (j,p,σ,η) -> DxGeo(M.manifold,p,data[j],σ/(α + σ),η)
   dPCd = (j,K,ξ,τ,η) -> ManTVproxDualCdifferential(N.manifold.manifold,
   getBase(n)[j,K],d,tuple(j,K),ξ,η)
   αs = [
   (nX,k) -> 0.,
   (nX,k) -> 1/5*nX*relerror,
   (nX,k) -> 1/(5*k)*nX*relerror
   ]

   fRs = Array{SnPoint{Float64}}[]
   rs = Array{Tuple}[]
   timing = zeros(3)

   start = p
   dualstart = ξ
   for (k,αk) in enumerate(αs)
      pr = ilrssnProblem(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,pPCd,dPCd,σ,τ,αk)

      # Optimize
      CPUtic()
      fR,r = ilRSSN_ML(pr,start,y0=dualstart,tol=1e-7/relerror,max_it=50)
      timing[k] = CPUtoc()

      push!(fRs,getValue(fR))
      push!(rs,r)
   end

   exp_save("res-ilRSSN-S2"*postfix, "res", (fRs,rs),
   "rprior", rprior,
   "orig", Orig,
   "data", Data,
   "relerror", relerror,
   "TVprior", cF(data),
   "TVstart", cF(p),
   "t0", t0,
   "timing", timing
   )

end
