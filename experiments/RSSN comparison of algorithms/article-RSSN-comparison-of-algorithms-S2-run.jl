include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))
include(joinpath("..","..","functions","ManTV.jl"))
include(joinpath("..","..","functions","ManTVdualbasepoint.jl"))
include(joinpath("..","..","functions","ManTVRproxes.jl"))
include(joinpath("..","..","functions","ManTVRproxCdifferentials.jl"))
include(joinpath("..","..","solvers","lRCPA.jl"))
include(joinpath("..","..","solvers","lRSSN-ML.jl"))

using Plots, LinearAlgebra, CPUTime, LaTeXStrings

pyplot()

function run_RSSN_comparison_algorithms_S2(λ;postfix="",return_result::Bool=false)

   if !(postfix=="")
      postfix = "-"*postfix
   end

   # Import data
   data = artificialS2RotationsImage(20)

   # Settings
   d = size(data)
   α = 1.5

   MM = Sphere(2)
   M = Power(MM, d)
   m = PowPoint(fill(SnPoint([0.,0.,1.]),d))

   N = TangentBundle(Power(MM,tuple(d...,length(d))))
   n = TBPoint(ManTVdualbasepoint(m),forwardLogs(M,m))

   anisotropic = false

   # lRCPA step

   cF = p -> ManTV(M,data,α,p,anisotropic)
   pP = (p,σ) -> proxDistance(M,σ/α,data,p)
   dP = (ξ,τ) -> ManTVRproxDual(N,n,ξ,λ*τ,anisotropic)
   Λ = p -> TBPoint(ManTVdualbasepoint(p),forwardLogs(M,p))
   DΛ = (p,η) -> TBTVector(zeroTVector(N.manifold,getBase(n)),DforwardLogs(M,p,η))
   DΛadj = (p,η) -> AdjDforwardLogs(M,p,getTangent(η))
   σ = τ = 0.35
   θ = 1.
   γ = 0.2
   lrcpaP = lrcpaProblem(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,σ,τ,θ,γ)

   # Optimize (for starting points)
   CPUtic()
   p,ξ,r = lRCPA(lrcpaP,data,return_dual=true,max_it=20,tol=1/2)
   t0 = CPUtoc()
   relerror = r[end][4]

   # lRSSN solve
   cF = p -> ManTV(M,data,α,p,anisotropic)
   pP = (p,σ) -> proxDistance(M,σ/α,data,p)
   dP = (ξ,τ) -> ManTVRproxDual(N,n,ξ,λ*τ,anisotropic)
   Λ = p -> TBPoint(ManTVdualbasepoint(p),forwardLogs(M,p))
   DΛ = (p,η) -> TBTVector(zeroTVector(N.manifold,getBase(n)),
                   DforwardLogs(M,p,η))
   DΛadj = (p,η) -> AdjDforwardLogs(M,p,getTangent(η))
   pPCd = (j,p,σ,η) -> DxGeo(M.manifold,p,data[j],σ/(α + σ),η)
   dPCd = (j,K,ξ,τ,η) -> ManTVRproxDualCdifferential(N.manifold.manifold,
                           getBase(n)[j,K],d,tuple(j,K),ξ,η,λ*τ,anisotropic)

   lrssnP = lrssnProblem(M,m,N,n,cF,pP,dP,Λ,DΛ,DΛadj,pPCd,dPCd,σ,τ)

   # Optimize
   start = p
   dualstart = ξ
   # tols = [0.4,0.3,0.2]
   tols = [1e-2,1e-4,1e-6]
   timing = zeros(2,length(tols),2)

   max_it = 100
   # max_it = 1

   # clear cache
   fR_lrcpa,r_lrcpa = lRCPA(lrcpaP,start,y0=dualstart,tol=1e-1,max_it=3)
   fR_lrssn,r_lrssn = lRSSN_ML(lrssnP,start,y0=dualstart,tol=1e-1,max_it=3)

   for (i,tol) in enumerate(tols)

      dbg(1,"#$(i) | Solve with lRCPA")
      CPUtic()
      fR_lrcpa,r_lrcpa = lRCPA(lrcpaP,start,y0=dualstart,tol=tol/relerror,max_it=100*max_it)
      lrcpat = CPUtoc()

      timing[1,i,1] = lrcpat #+ t0
      timing[1,i,2] = length(r_lrcpa)

      dbg(1,"#$(i) | Solve with RSSN")
      CPUtic()
      fR_lrssn,r_lrssn = lRSSN_ML(lrssnP,start,y0=dualstart,tol=tol/relerror,max_it=max_it)
      lrssnt = CPUtoc()

      timing[2,i,1] = lrssnt #+ t0
      timing[2,i,2] = length(r_lrssn)

   end

   exp_save("res-S2"*postfix,
      "data", data,
      "TVstart", cF(data),
      "TVp", cF(start),
      "Xlrcpa", cF(fR_lrcpa),
      "Xlrssn", cF(fR_lrssn),
      "tols", tols,
      "relerror", relerror,
      "timing", timing,
      "t0", t0,
      "r-start", r,
      "lrcpa-res", tuple(fR_lrcpa,r_lrcpa),
      "lrssn-res", tuple(fR_lrssn,r_lrssn))

   if return_result
      return timing, (fR_lrcpa,r_lrcpa), (fR_lrssn,r_lrssn)
   end

end
