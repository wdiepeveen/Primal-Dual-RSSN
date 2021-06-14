include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))
include(joinpath("..","..","functions","ManTV.jl"))
include(joinpath("..","..","functions","ManTVdualbasepoint.jl"))
include(joinpath("..","..","functions","ManTVproxes.jl"))
include(joinpath("..","..","functions","ManTVproxCdifferentials.jl"))
include(joinpath("..","..","solvers","lRCPA.jl"))
include(joinpath("..","..","solvers","lRSSN-ML.jl"))

using Manopt

function run_S2_RSSN(ℓ,α,σ,τ;postfix="",max_it=20,tol=1e-8,dual_prior=false)

   dbg(1,"solve S2  with RSSN")

   if !(postfix=="")
      postfix = "-"*postfix
   end

   # Settings
   d = (2*ℓ,)

   MM = Sphere(2)
   M = Power(MM, d)
   m = PowPoint(fill(SnPoint([1., 0., 0.]),d))
   N = TangentBundle(Power(MM,d))
   n = TBPoint(ManTVdualbasepoint(m),forwardLogs(M,m))

   f₁ = SnPoint(1/sqrt(2)*[1., 1., 0.])
   f₂ = SnPoint(1/sqrt(2)*[1., -1., 0.])
   f = vcat(fill(f₁,ℓ),fill(f₂,ℓ))
   data = PowPoint(f)

   δ = min(1/2, α/(ℓ*distance(MM,f₁,f₂)))
   exact = PowPoint(vcat(fill(geodesic(MM,f₁,f₂,δ),ℓ), fill(geodesic(MM,f₂,f₁,δ),ℓ)))

   dualstart = zeroTVector(N,n)

   if dual_prior==true

      cF = p -> ManTV(M,data,α,p)
      pP = (p,σ) -> proxDistance(M,σ/α,data,p)
      dP = (ξ,τ) -> ManTVproxDual(N,n,ξ)
      Λ = p -> TBPoint(ManTVdualbasepoint(p),forwardLogs(M,p))
      DΛ = (p,η) -> TBTVector(zeroTVector(N.manifold,getBase(n)),DforwardLogs(M,p,η))
      DΛadj = (p,η) -> AdjDforwardLogs(M,p,getTangent(η))

      pr = lrcpaProblem(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,σ,τ)

      # Optimize
      fR,dualstart,r = lRCPA(pr,data,return_dual=true,max_it=1)
   end

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

   pr = lrssnProblem(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,pPCd,dPCd,σ,τ)

   # Optimize
   start = data
   fR = m
   r = Tuple[]
   CPUtic()
   if dual_prior==true
      fR,r = lRSSN_ML(pr,start,y0=dualstart,max_it=max_it,tol=tol)
   else
      fR,r = lRSSN_ML(pr,start,max_it=max_it,tol=tol)
   end
   t = CPUtoc()

   TVstart = cF(start)
   error = distance(M,fR,exact)

   exp_save("res-RSSN-S2"*postfix, "res", (fR,r),
   "data", data,
   "exact", exact,
   "TVstart", TVstart,
   "t", t,
   "error", error)

end
