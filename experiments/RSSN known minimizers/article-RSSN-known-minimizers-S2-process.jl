include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))
include(joinpath("..","..","functions","ManTV.jl"))
include(joinpath("..","..","functions","ManTVdualbasepoint.jl"))
include(joinpath("..","..","functions","ManTVproxes.jl"))
include(joinpath("..","..","functions","ManTVproxCdifferentials.jl"))
# include(joinpath("..","..","solvers","eRSSN-ML.jl"))

using Manopt, Plots, LaTeXStrings

pyplot()

function process_S2(;postfix="",view_result::Bool=true,path::Union{String,Missing}=missing)

   if !(postfix=="")
      postfix = "-"*postfix
   end

   if ismissing(path)
      resc = load(exp_prefix() *"res-RSSN-S2-cold.jld")
      resw = load(exp_prefix() *"res-RSSN-S2-warm.jld")
   else
      resc = load(path *"res-RSSN-S2-cold.jld")
      resw = load(path *"res-RSSN-S2-warm.jld")
      cp(path *"res-RSSN-S2-cold.jld",exp_prefix() *"res-RSSN-S2-cold.jld")
      cp(path *"res-RSSN-S2-warm.jld",exp_prefix() *"res-RSSN-S2-warm.jld")
   end

   data = resc["data"]
   exact = resc["exact"]

   TVstart = resc["TVstart"]

   tc = resc["t"]
   erc = resc["error"]
   tw = resw["t"]
   erw = resw["error"]

   fRc,rc = resc["res"]
   fRw,rw = resw["res"]

   dbg(1,"process S2 solutions")

   table = "runtime cold start: $(tc) |runtime warm start: $(tw)"


   if view_result
      renderAsymptote(exp_prefix()*"original-S2.asy", asyExportS2Data; data=data, render=4)
      renderAsymptote(exp_prefix()*"exact-S2.asy", asyExportS2Data; data=exact, render=4)
      renderAsymptote(exp_prefix()*"result-RSSN-S2-cold.asy", asyExportS2Data; data=fRc, render=4)
      renderAsymptote(exp_prefix()*"result-RSSN-S2-warm.asy", asyExportS2Data; data=fRw, render=4)
   else
      asyExportS2Data(exp_prefix()*"original-S2.asy";data=data)
      asyExportS2Data(exp_prefix()*"exact-S2.asy";data=exact)
      asyExportS2Data(exp_prefix()*"result-RSSN-S2-cold.asy";data=fRc)
      asyExportS2Data(exp_prefix()*"result-RSSN-S2-warm.asy";data=fRw)
   end

   # Results
   ??c = zeros(size(rc)[1],4)
   for i in 1:size(rc)[1]
      ??c[i,1] = rc[i][1]
      ??c[i,2] = rc[i][2]
      ??c[i,3] = rc[i][3]+1e-16
      ??c[i,4] = rc[i][4]+1e-16
   end

   ??w = zeros(size(rw)[1],4)
   for i in 1:size(rw)[1]
      ??w[i,1] = rw[i][1]
      ??w[i,2] = rw[i][2]
      ??w[i,3] = rw[i][3]+1e-16
      ??w[i,4] = rw[i][4]+1e-16
   end

   costfig = plot([0;??c[:,1]], [TVstart;??c[:,2]], lab="cold start")
   plot!(costfig, [0;??w[:,1]], [TVstart;??w[:,2]], lab="warm start")
   xlabel!("Iterations")
   ylabel!("Cost")
   xlims!((0,10))

   changefig = plot(??c[:,1], ??c[:,3], yaxis=:log, lab="cold start")
   plot!(changefig, ??w[:,1], ??w[:,3], yaxis=:log, lab="warm start")
   xlabel!("Iterations")
   ylabel!("Change")
   xlims!((0,10))

   rerrorfig = plot([0;??c[:,1]], [1;??c[:,4]], yaxis=:log, lab="cold start")
   plot!(rerrorfig, [0;??w[:,1]], [1;??w[:,4]], yaxis=:log, lab="warm start")
   xlabel!("Iterations")
   ylabel!(L"\epsilon_{rel}")
   xlims!((0,10))

   exp_savetable(table,"runtimes-S2"*postfix)

   exp_savefig(costfig, "cost-S2"*postfix, true)

   exp_savefig(changefig, "change-S2"*postfix, true)

   exp_savefig(rerrorfig, "relerror-S2"*postfix, true)

   return [erc, erw]
end
