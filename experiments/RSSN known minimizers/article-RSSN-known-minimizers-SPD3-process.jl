include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))

using Manopt, Plots, LaTeXStrings

pyplot()

function process_SPD3(;postfix="",view_result::Bool=true,path::Union{String,Missing}=missing)

   if !(postfix=="")
      postfix = "-"*postfix
   end

   if ismissing(path)
      resc = load(exp_prefix() *"res-RSSN-SPD3-cold.jld")
      resw = load(exp_prefix() *"res-RSSN-SPD3-warm.jld")
   else
      resc = load(path *"res-RSSN-SPD3-cold.jld")
      resw = load(path *"res-RSSN-SPD3-warm.jld")
      cp(path *"res-RSSN-SPD3-cold.jld",exp_prefix() *"res-eRSSN-SPD3-cold.jld")
      cp(path *"res-RSSN-SPD3-warm.jld",exp_prefix() *"res-lRSSN-SPD3-warm.jld")
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

   dbg(1,"process SPD solutions")

   table = "runtime cold start: $(tc) |runtime warm start: $(tw)"


   if view_result
      renderAsymptote(exp_prefix()*"original-SPD3.asy", asyExportSPDData; data=data, render=4)
      renderAsymptote(exp_prefix()*"exact-SPD3.asy", asyExportSPDData; data=exact, render=4)
      renderAsymptote(exp_prefix()*"result-RSSN-SPD3-cold.asy", asyExportSPDData; data=fRc, render=4)
      renderAsymptote(exp_prefix()*"result-RSSN-SPD3-warm.asy", asyExportSPDData; data=fRw, render=4)
   else
      asyExportSPDData(exp_prefix()*"original-SPD3.asy";data=data)
      asyExportSPDData(exp_prefix()*"exact-SPD3.asy";data=exact)
      asyExportSPDData(exp_prefix()*"result-RSSN-SPD3-cold.asy";data=fRc)
      asyExportSPDData(exp_prefix()*"result-RSSN-SPD3-warm.asy";data=fRw)
   end

   # Results
   ϵc = zeros(size(rc)[1],4)
   for i in 1:size(rc)[1]
      ϵc[i,1] = rc[i][1]
      ϵc[i,2] = rc[i][2]
      ϵc[i,3] = rc[i][3]+1e-16
      ϵc[i,4] = rc[i][4]+1e-16
   end

   ϵw = zeros(size(rw)[1],4)
   for i in 1:size(rw)[1]
      ϵw[i,1] = rw[i][1]
      ϵw[i,2] = rw[i][2]
      ϵw[i,3] = rw[i][3]+1e-16
      ϵw[i,4] = rw[i][4]+1e-16
   end

   costfig = plot([0;ϵc[:,1]], [TVstart;ϵc[:,2]], lab="cold start")
   plot!(costfig, [0;ϵw[:,1]], [TVstart;ϵw[:,2]], lab="warm start")
   xlabel!("Iterations")
   ylabel!("Cost")
   xlims!((0,10))

   changefig = plot(ϵc[:,1], ϵc[:,3], yaxis=:log, lab="cold start")
   plot!(changefig, ϵw[:,1], ϵw[:,3], yaxis=:log, lab="warm start")
   xlabel!("Iterations")
   ylabel!("Change")
   xlims!((0,10))

   rerrorfig = plot([0;ϵc[:,1]], [1;ϵc[:,4]], yaxis=:log, lab="cold start")
   plot!(rerrorfig, [0;ϵw[:,1]], [1;ϵw[:,4]], yaxis=:log, lab="warm start")
   xlabel!("Iterations")
   ylabel!(L"\epsilon_{rel}")
   xlims!((0,10))

   exp_savetable(table,"runtimes-SPD3"*postfix)

   exp_savefig(costfig, "cost-SPD3"*postfix, true)

   exp_savefig(changefig, "change-SPD3"*postfix, true)

   exp_savefig(rerrorfig, "relerror-SPD3"*postfix, true)

   return [erc, erw]

end
