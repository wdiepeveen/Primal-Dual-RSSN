include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))
include(joinpath("..","..","functions","ManTV.jl"))
include(joinpath("..","..","functions","ManTVdualbasepoint.jl"))
include(joinpath("..","..","functions","ManTVRproxes.jl"))
include(joinpath("..","..","functions","ManTVproxCdifferentials.jl"))
include(joinpath("..","..","solvers","lRCPA.jl"))
include(joinpath("..","..","solvers","ilRSSN-ML.jl"))

using Plots, Manopt, Images, Random, ColorSchemes, LaTeXStrings

pyplot()

function process_RSSN_inexact_S2(;postfix="",view_result::Bool=true,path::Union{String,Missing}=missing)

   if !(postfix=="")
      postfix = "-"*postfix
   end

   if ismissing(path)
      res = load(exp_prefix() *"res-ilRSSN-S2"*postfix*".jld")
   else
      res = load(path *"res-ilRSSN-S2"*postfix*".jld")
      cp(path *"res-ilRSSN-S2"*postfix*".jld",exp_prefix() *"res-ilRSSN-S2"*postfix*".jld")
   end

   fRs,rs = res["res"]
   rprior = res["rprior"]

   orig = res["orig"]
   data = res["data"]

   TVprior = res["TVprior"]
   TVstart = res["TVstart"]
   relerror = res["relerror"]

   t0 = res["t0"]
   timing = res["timing"]

   dataColor = RGBA{Float64}(colorant"#BBBBBB")
   s2dColor = RGBA{Float64}(colorant"#EE7733") # data Color: Tol Vibrant Orange
   s1Color = RGBA{Float64}(colorant"#0077BB") # control point data color: Tol Virbant Blue

   # Results

   ϵₗ = zeros(length(rprior),4)
   for i in 1:length(rprior)
      ϵₗ[i,1] = rprior[i][1]
      ϵₗ[i,2] = rprior[i][2]
      ϵₗ[i,3] = rprior[i][3]+1e-16
      ϵₗ[i,4] = rprior[i][4]+1e-16
   end
   L = length(rprior)

   ϵ₁ = zeros(length(rs[1]),4)
   cRate₁ = zeros(length(rs[1]))
   for i in 1:length(rs[1])
      ϵ₁[i,1] = L+rs[1][i][1]
      ϵ₁[i,2] = rs[1][i][2]
      ϵ₁[i,3] = rs[1][i][3]+1e-16
      ϵ₁[i,4] = (rs[1][i][4]+1e-16)*relerror
      if i >= 3
         cRate₁[i] = log(ϵ₁[i,4]/ϵ₁[i-1,4])/log(ϵ₁[i-1,4]/ϵ₁[i-2,4])
      elseif i==2
         cRate₁[i] = log(ϵ₁[i,4]/ϵ₁[i-1,4])/log(ϵ₁[i-1,4]/ϵₗ[L,4])
      elseif i==1
         cRate₁[i] = log(ϵ₁[i,4]/ϵₗ[L,4])/log(ϵₗ[L,4]/ϵₗ[L-1,4])
      end
   end

   ϵ₂ = zeros(length(rs[2]),4)
   cRate₂ = zeros(length(rs[2]))
   for i in 1:length(rs[2])
      ϵ₂[i,1] = L+rs[2][i][1]
      ϵ₂[i,2] = rs[2][i][2]
      ϵ₂[i,3] = rs[2][i][3]+1e-16
      ϵ₂[i,4] = (rs[2][i][4]+1e-16)*relerror
      if i >= 3
         cRate₂[i] = log(ϵ₂[i,4]/ϵ₂[i-1,4])/log(ϵ₂[i-1,4]/ϵ₂[i-2,4])
      elseif i==2
         cRate₂[i] = log(ϵ₂[i,4]/ϵ₂[i-1,4])/log(ϵ₂[i-1,4]/ϵₗ[L,4])
      elseif i==1
         cRate₂[i] = log(ϵ₂[i,4]/ϵₗ[L,4])/log(ϵₗ[L,4]/ϵₗ[L-1,4])
      end
   end

   ϵ₃ = zeros(length(rs[3]),4)
   cRate₃ = zeros(length(rs[3]))
   for i in 1:length(rs[3])
      ϵ₃[i,1] = L+rs[3][i][1]
      ϵ₃[i,2] = rs[3][i][2]
      ϵ₃[i,3] = rs[3][i][3]+1e-16
      ϵ₃[i,4] = (rs[3][i][4]+1e-16)*relerror
      if i >= 3
         cRate₃[i] = log(ϵ₃[i,4]/ϵ₃[i-1,4])/log(ϵ₃[i-1,4]/ϵ₃[i-2,4])
      elseif i==2
         cRate₃[i] = log(ϵ₃[i,4]/ϵ₃[i-1,4])/log(ϵ₃[i-1,4]/ϵₗ[L,4])
      elseif i==1
         cRate₃[i] = log(ϵ₃[i,4]/ϵₗ[L,4])/log(ϵₗ[L,4]/ϵₗ[L-1,4])
      end
   end


   costfig = plot([0;ϵₗ[:,1]], [TVprior;ϵₗ[:,2]], lab="")

   plot!(costfig, [L;ϵ₁[:,1]], [TVstart;ϵ₁[:,2]], lab="a₁")

   plot!(costfig, [L;ϵ₂[:,1]], [TVstart;ϵ₂[:,2]], lab="a₂")

   plot!(costfig, [L;ϵ₃[:,1]], [TVstart;ϵ₃[:,2]], lab="a₃")
   xlabel!(costfig,"Iterations")
   ylabel!(costfig,"Cost")

   # changefig = plot(ϵ[:,1], ϵ[:,3], yaxis=:log, lab="")

   relerrorfig = plot([0;ϵₗ[:,1]], [1.;ϵₗ[:,4]], yaxis=:log, lab="")

   plot!(relerrorfig, [L;ϵ₁[:,1]], [relerror;ϵ₁[:,4]], yaxis=:log, lab="a₁")

   plot!(relerrorfig, [L;ϵ₂[:,1]], [relerror;ϵ₂[:,4]], yaxis=:log, lab="a₂")

   plot!(relerrorfig, [L;ϵ₃[:,1]], [relerror;ϵ₃[:,4]], yaxis=:log, lab="a₃")
   xlabel!(relerrorfig,"Iterations")
   ylabel!(relerrorfig,L"\epsilon_{rel}")

   convfig = plot([0;ϵ₂[:,1] .- L], ones(length(rs[2])+1),lab="q=1")

   scatter!(convfig,ϵ₁[:,1] .- L,cRate₁, markersize=8, lab="a₁")

   scatter!(convfig,ϵ₂[:,1] .- L, cRate₂, markersize=7, lab="a₂")

   scatter!(convfig,ϵ₃[:,1] .- L, cRate₃, markersize=6, lab="a₃")

   xlims!( (0,maximum(length.(rs))) )
   xlabel!(convfig,"Iterations")
   ylabel!(convfig,"q")
   ylims!((0,2.5))

   asyExportS2Data(exp_prefix()*"original-S2.asy"; data=PowPoint(orig))
   asyExportS2Data(exp_prefix()*"data-S2.asy"; data=PowPoint(data))

   for (k,fR) in enumerate(fRs)

      if view_result
         renderAsymptote(exp_prefix()*"result-S2-alpha-$(k)"*postfix*".asy",
         asyExportS2Signals;
         points=[data,fR],
         curves=[orig],
         colors=Dict(:points => [s2dColor,s1Color],:curves => [dataColor]),
         dotSize = 2.5,lineWidth = 0.75, cameraPosition = (2.5,1.,.5),
         render=4)
      else
         asyExportS2Signals(exp_prefix()*"result-S2-alpha-$(k)"*postfix*".asy",
         points=[data,fR],
         curves=[orig],
         colors=Dict(:points => [s2dColor,s1Color],:curves => [dataColor]),
         dotSize = 2.5,lineWidth = 0.75, cameraPosition = (2.5,1.,.5))
      end
   end

   table = "t₀ = $(t0) | t₁ = $(timing[1]) | t₂ = $(timing[2]) | t₃ = $(timing[3])"

   exp_savefig(costfig, "cost-S2"*postfix, true)

   # exp_savefig(changefig, "change", true)

   exp_savefig(relerrorfig, "relerror-S2"*postfix, true)

   exp_savefig(convfig, "convrate-S2"*postfix, true)

   exp_savetable(table, "timing-S2"*postfix)

end
