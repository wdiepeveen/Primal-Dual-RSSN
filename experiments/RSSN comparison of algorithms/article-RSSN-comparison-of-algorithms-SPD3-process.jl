include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))

using Plots, CPUTime, LaTeXStrings

pyplot()

function process_RSSN_comparison_algorithms_SPD3(λ;postfix="",view_result::Bool=false,path::Union{String,Missing}=missing)

   if !(postfix=="")
      postfix = "-"*postfix
   end

   max_xlim = 500

   dataColor = RGBA{Float64}(colorant"#BBBBBB")
   s2dColor = RGBA{Float64}(colorant"#EE7733") # data Color: Tol Vibrant Orange
   s1Color = RGBA{Float64}(colorant"#0077BB") # control point data color: Tol Virbant Blue
   pdhgColor = RGBA{Float64}(colorant"#33BBEE") # tangent vector: Tol Vibrant Teal
   ssnColor = RGBA{Float64}(colorant"#EE2211") # tangent vector: Tol Vibrant Teal
   assnColor = RGBA{Float64}(colorant"#33BB22") # tangent vector: Tol Vibrant Teal

   # Results
   if ismissing(path)
      res = load(exp_prefix() *"res-SPD3" * postfix * ".jld")
   else
      res = load(path *"res-SPD3" * postfix * ".jld")
      cp(path *"res-SPD3" * postfix * ".jld",exp_prefix() *"res-SPD3" * postfix * ".jld")
   end
   data = res["data"]

   TVstart = res["TVstart"]
   TVp = res["TVp"]
   TVlrcpa = res["Xlrcpa"]
   TVlrssn = res["Xlrssn"]
   tols = res["tols"]
   relerror = res["relerror"]

   timing = res["timing"]
   t0 = res["t0"]
   r = res["r-start"]
   fR_lrcpa,r_lrcpa = res["lrcpa-res"]
   fR_lrssn,r_lrssn = res["lrssn-res"]

   sc = 5
   scaling = (sc, sc, sc)
   if view_result
      renderAsymptote(exp_prefix()*"original-SPD3.asy", asyExportSPDData; data=data, render=4, scaleAxes = scaling)
      renderAsymptote(exp_prefix()*"result-lRCPA-SPD3"*postfix*".asy", asyExportSPDData; data=fR_lrcpa, render=4, scaleAxes = scaling)
      renderAsymptote(exp_prefix()*"result-RSSN-SPD3"*postfix*".asy", asyExportSPDData; data=fR_lrssn, render=4, scaleAxes = scaling)
   else
      asyExportSPDData(exp_prefix()*"original-SPD3.asy";data=data, scaleAxes = scaling)
      asyExportSPDData(exp_prefix()*"result-lRCPA-SPD3"*postfix*".asy";data=fR_lrcpa, scaleAxes = scaling)
      asyExportSPDData(exp_prefix()*"result-RSSN-SPD3"*postfix*".asy";data=fR_lrssn, scaleAxes = scaling)
   end

   table = raw"\begin{table}[h!]
   \begin{tabular}{l|l|l|c|l|c|l|c}
   \hline \hline
   $\lambda$ = "* "$(λ)" * raw" & $ t_0 $ =  " *"$(t0)" *raw"&\multicolumn{2}{l|}{$\epsilon = 1e-2$} & \multicolumn{2}{l|}{$\epsilon = 1e-4$} & \multicolumn{2}{l}{$\epsilon = 1e-6$} \\ \hline
   Method & $X(p^*)$ & Time  & \# Iterations  & Time   & \# Iterations   & Time   & \# Iterations  \\ \hline
   " *
   "lRCPA  & $(Float32(TVlrcpa)) &  $(timing[1,1,1]) &  $(Int(timing[1,1,2]))  &  $(timing[1,2,1]) &  $(Int(timing[1,2,2])) & $(timing[1,3,1]) & $(Int(timing[1,3,2])) "*
   raw" \\ \hline
   "*
   "lRSSN  & $(Float32(TVlrssn))  &  $(timing[2,1,1]) &  $(Int(timing[2,1,2])) & $(timing[2,2,1]) & $(Int(timing[2,2,2]))  & $(timing[2,3,1]) & $(Int(timing[2,3,2]))  " *
   raw"\\ \hline \hline
   \end{tabular}
   \end{table}"

   L = length(r)

   prior = zeros(L,3)
   prior_cRate = zeros(L-2)
   prior[:,1] = collect(1:L)
   for i in 1:L
      prior[i,2] = r[i][2]
      prior[i,3] = r[i][4]+1e-16
      if i >= 3
         prior_cRate[i-2] = log(prior[i,3]/prior[i-1,3])/log(prior[i-1,3]/prior[i-2,3])
      end
   end


   lrcpa = zeros(size(r_lrcpa)[1],3)
   lrcpa_cRate = zeros(length(r_lrcpa))
   lrcpa[:,1] = collect(L+1:L+size(r_lrcpa)[1])
   for i in 1:size(r_lrcpa)[1]
      lrcpa[i,2] = r_lrcpa[i][2]
      lrcpa[i,3] = (r_lrcpa[i][4]+1e-16)*relerror
      if i >= 3
         lrcpa_cRate[i] = log(lrcpa[i,3]/lrcpa[i-1,3])/log(lrcpa[i-1,3]/lrcpa[i-2,3])
      elseif i==2
         lrcpa_cRate[i] = log(lrcpa[i,3]/lrcpa[i-1,3])/log(lrcpa[i-1,3]/prior[L,3])
      elseif i==1
         lrcpa_cRate[i] = log(lrcpa[i,3]/prior[L,3])/log(prior[L,3]/prior[L-1,3])
      end
   end

   lrssn = zeros(size(r_lrssn)[1],3)
   lrssn_cRate = zeros(length(r_lrssn))
   lrssn[:,1] = collect(L+1:L+size(r_lrssn)[1])
   for i in 1:size(r_lrssn)[1]
      lrssn[i,2] = r_lrssn[i][2]
      lrssn[i,3] = (r_lrssn[i][4]+1e-16)*relerror
      if i >= 3
         lrssn_cRate[i-2] = log(lrssn[i,3]/lrssn[i-1,3])/log(lrssn[i-1,3]/lrssn[i-2,3])
      elseif i==2
         lrssn_cRate[i] = log(lrssn[i,3]/lrssn[i-1,3])/log(lrssn[i-1,3]/prior[L,3])
      elseif i==1
         lrssn_cRate[i] = log(lrssn[i,3]/prior[L,3])/log(prior[L,3]/prior[L-1,3])
      end
   end

   costfig = plot([L;lrcpa[:,1]], [TVp;lrcpa[:,2]], lab="lRCPA")

   plot!(costfig,[L;lrssn[:,1]], [TVp;lrssn[:,2]], lab="PD-RSSN (proposed)")
   plot!(costfig,[0;prior[:,1]], [TVstart;prior[:,2]], lab="")
   xlims!( (0,L+min(max_xlim,max(length(r_lrssn),250))) )
   xlabel!(costfig,"Iterations")
   ylabel!(costfig,"Cost")

   convfig = scatter(lrcpa[:,1],lrcpa_cRate, markersize=2, lab="lRCPA")

   scatter!(convfig,lrssn[:,1], lrssn_cRate, markersize=6, lab="PD-RSSN (proposed)")

   scatter!(convfig,[0;prior[:,1]],[NaN;NaN;NaN;prior_cRate], markersize=2, lab="")
   xlims!( (0,L+min(max_xlim,max(length(r_lrssn),250))) )
   xlabel!(convfig,"Iterations")
   ylabel!(convfig,"q")
   ylims!((0,3))

   rerrorfig = plot([L;lrcpa[:,1]], [relerror;lrcpa[:,3]], yaxis=:log, lab="lRCPA")

   plot!(rerrorfig,[L;lrssn[:,1]], [relerror;lrssn[:,3]], yaxis=:log, lab="PD-RSSN (proposed)")

   plot!(rerrorfig,[0;prior[:,1]], [1;prior[:,3]], yaxis=:log, lab="")
   xlims!( (0,L+min(max_xlim,max(length(r_lrssn),250))) )
   xlabel!(rerrorfig,"Iterations")
   ylabel!(rerrorfig,L"\epsilon_{rel}")

   timingfig = plot(tols, timing[1,:,1], xaxis=:log, linecolor = pdhgColor, lab="lRCPA")

   plot!(timingfig, tols, timing[2,:,1], xaxis=:log, linecolor = ssnColor, lab="PD-RSSN (proposed)")

   exp_savetable(table, "result-SPD3"*postfix)

   exp_savefig(costfig, "cost-SPD3"*postfix, true)

   exp_savefig(convfig, "convrate-SPD3"*postfix, true)

   exp_savefig(rerrorfig, "relerror-SPD3"*postfix, true)

   exp_savefig(timingfig, "timing-SPD3"*postfix, true)

   exp_end()

end
