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

function run_RSSN_scaling_S2(λ,scales;postfix="",return_result::Bool=false)

   if !(postfix=="")
      postfix = "-"*postfix
   end

   timing = zeros(2,length(scales),2)
   t0 = zeros(length(scales),2)

   for (i,scale) in enumerate(scales)
      dbg(1,"Scale = $(scale)")


   # Import data
      data = artificialS2RotationsImage(scale)  # default 20

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

      cF = p -> ManTV(M,data,α,p,anisotropic) # /scale^2
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
      dbg(1,"Doing Pre-steps")
      CPUtic()
      p,ξ,r = lRCPA(lrcpaP,data,return_dual=true,max_it=1000,tol=1/2)
      t0i = CPUtoc()
      t0[i,1] = t0i
      t0[i,2] =  length(r)
      relerror = r[end][4]

      # lRSSN solve
      cF = p -> ManTV(M,data,α,p,anisotropic) #/scale^2
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

      tol = 1e-6  # TODO was 1e-6
      max_it = 25 # TODO was 100
      # max_it = 1

      # clear cache
      dbg(1,"Clear cache")
      fR_lrcpa,r_lrcpa = lRCPA(lrcpaP,start,y0=dualstart,tol=1e-1,max_it=3)
      fR_lrssn,r_lrssn = lRSSN_ML(lrssnP,start,y0=dualstart,tol=1e-1,max_it=3)


      dbg(1,"# $(i) | Solve with lRCPA")
      CPUtic()
      fR_lrcpa,r_lrcpa = lRCPA(lrcpaP,start,y0=dualstart,tol=tol/relerror,max_it=1)  # was 50*max_it
      lrcpat = CPUtoc()

      timing[1,i,1] = lrcpat #+ t0
      timing[1,i,2] = length(r_lrcpa)

      dbg(1,"# $(i) | Solve with RSSN")
      CPUtic()
      fR_lrssn,r_lrssn = lRSSN_ML(lrssnP,start,y0=dualstart,tol=tol/relerror,max_it=max_it)
      lrssnt = CPUtoc()

      timing[2,i,1] = lrssnt #+ t0
      timing[2,i,2] = length(r_lrssn)

      exp_save("res-S2-scale-$(scale)"*postfix,
         "data", data,
         "TVstart", cF(data),
         "TVp", cF(start),
         "Xlrcpa", cF(fR_lrcpa),
         "Xlrssn", cF(fR_lrssn),
         "relerror", relerror,
         # "timing", timing,
         # "t0", t0,
         "r-start", r,
         "lrcpa-res", tuple(fR_lrcpa,r_lrcpa),
         "lrssn-res", tuple(fR_lrssn,r_lrssn))

   end

   exp_save("res-S2"*postfix,
      "timing", timing,
      "t0", t0)

   if return_result
      return timing, (fR_lrcpa,r_lrcpa), (fR_lrssn,r_lrssn)
   end

end
