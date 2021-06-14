include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))

using Plots, CPUTime, LaTeXStrings

pyplot()

function process_RSSN_scaling_SPD3(λ,scales;postfix="",view_result::Bool=false,path::Union{String,Missing}=missing)

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
   # First non-scalin related
   if ismissing(path)
      res = load(exp_prefix() *"res-SPD3" * postfix * ".jld")
   else
      res = load(path *"res-SPD3" * postfix * ".jld")
      cp(path *"res-SPD3" * postfix * ".jld",exp_prefix() *"res-SPD3" * postfix * ".jld")
   end

   timing = res["timing"]
   t0 = res["t0"]

   table = raw"\begin{table}[h!]
   \begin{tabular}{c|c|c|c|c|c|c}
   \hline \hline
   $\lambda =$ "* "$(λ)" * raw" &\multicolumn{2}{c|}{Pre-steps}  &\multicolumn{2}{c|}{lRCPA} & \multicolumn{2}{c}{PD-RSSN} \\ \hline
   d & Time  & \# Iterations  & Time   & \# Iterations   & Time   & \# Iterations  \\ \hline
   "
   for (i,scale) in enumerate(scales)
      table *="$(scale) & $(t0[i,1]) & $(Int(t0[i,2])) & $(timing[1,i,1]) &  $(Int(timing[1,i,2]))  &  $(timing[2,i,1]) &  $(Int(timing[2,i,2]))"*
   raw" \\ \hline
   "
   end
   table *= raw" \hline
   \end{tabular}
   \end{table}"

   # Scaling related
   costfig = plot()
   convfig = scatter()
   rerrorfig = plot()

   for scale in scales
      if ismissing(path)
         res = load(exp_prefix() *"res-SPD3-scale-$(scale)" * postfix * ".jld")
      else
         res = load(path *"res-SPD3-scale-$(scale)" * postfix * ".jld")
         cp(path *"res-SPD3-scale-$(scale)" * postfix * ".jld",exp_prefix() *"res-SPD3-scale-$(scale)" * postfix * ".jld")
      end


      data = res["data"]

      TVstart = res["TVstart"]
      TVp = res["TVp"]
      TVlrcpa = res["Xlrcpa"]
      TVlrssn = res["Xlrssn"]
      relerror = res["relerror"]

      r = res["r-start"]
      fR_lrcpa,r_lrcpa = res["lrcpa-res"]
      fR_lrssn,r_lrssn = res["lrssn-res"]

      if view_result
         renderAsymptote(exp_prefix()*"original-SPD3-scale-$(scale).asy", asyExportSPDData; data=data, render=4)
         renderAsymptote(exp_prefix()*"result-lRCPA-SPD3-scale-$(scale)"*postfix*".asy", asyExportSPDData; data=fR_lrcpa, render=4)
         renderAsymptote(exp_prefix()*"result-RSSN-SPD3-scale-$(scale)"*postfix*".asy", asyExportSPDData; data=fR_lrssn, render=4)
      else
         asyExportSPDData(exp_prefix()*"original-SPD3-scale-$(scale).asy";data=data)
         asyExportSPDData(exp_prefix()*"result-lRCPA-SPD3-scale-$(scale)"*postfix*".asy";data=fR_lrcpa)
         asyExportSPDData(exp_prefix()*"result-RSSN-SPD3-scale-$(scale)"*postfix*".asy";data=fR_lrssn)
      end

      L = length(r)

      prior = zeros(L,3)
      # prior_cRate = zeros(L-2)
      prior[:,1] = collect(1:L)
      for i in 1:L
         prior[i,2] = r[i][2]
         prior[i,3] = r[i][4]+1e-16
      #    if i >= 3
      #       prior_cRate[i-2] = log(prior[i,3]/prior[i-1,3])/log(prior[i-1,3]/prior[i-2,3])
      #    end
      end


      lrcpa = zeros(size(r_lrcpa)[1],3)
      # lrcpa_cRate = zeros(length(r_lrcpa))
      lrcpa[:,1] = collect(L+1:L+size(r_lrcpa)[1])
      for i in 1:size(r_lrcpa)[1]
         lrcpa[i,2] = r_lrcpa[i][2]
         lrcpa[i,3] = (r_lrcpa[i][4]+1e-16)*relerror
      #    if i >= 3
      #       lrcpa_cRate[i] = log(lrcpa[i,3]/lrcpa[i-1,3])/log(lrcpa[i-1,3]/lrcpa[i-2,3])
      #    elseif i==2
      #       lrcpa_cRate[i] = log(lrcpa[i,3]/lrcpa[i-1,3])/log(lrcpa[i-1,3]/prior[L,3])
      #    elseif i==1
      #       lrcpa_cRate[i] = log(lrcpa[i,3]/prior[L,3])/log(prior[L,3]/prior[L-1,3])
      #    end
      end

      lrssn = zeros(size(r_lrssn)[1],3)
      # lrssn_cRate = zeros(length(r_lrssn))
      lrssn[:,1] = collect(L+1:L+size(r_lrssn)[1])
      for i in 1:size(r_lrssn)[1]
         lrssn[i,2] = r_lrssn[i][2]
         lrssn[i,3] = (r_lrssn[i][4]+1e-16)*relerror
      #    if i >= 3
      #       lrssn_cRate[i-2] = log(lrssn[i,3]/lrssn[i-1,3])/log(lrssn[i-1,3]/lrssn[i-2,3])
      #    elseif i==2
      #       lrssn_cRate[i] = log(lrssn[i,3]/lrssn[i-1,3])/log(lrssn[i-1,3]/prior[L,3])
      #    elseif i==1
      #       lrssn_cRate[i] = log(lrssn[i,3]/prior[L,3])/log(prior[L,3]/prior[L-1,3])
      #    end
      end

      plot!(costfig,[L;lrcpa[:,1]], [TVp;lrcpa[:,2]], lab="lRCPA scale=$(scale)")
      plot!(costfig,[L;lrssn[:,1]], [TVp;lrssn[:,2]], lab="PD-RSSN scale=$(scale)")
      plot!(costfig,[0;prior[:,1]], [TVstart;prior[:,2]], lab="")
      xlims!( (0,L+min(max_xlim,max(length(r_lrssn),round(length(r_lrcpa)/10)))) )
      xlabel!(costfig,"Iterations")
      ylabel!(costfig,"Cost")

      # convfig = scatter(lrcpa[:,1],lrcpa_cRate, markersize=2, lab="lRCPA")
      #
      # scatter!(convfig,lrssn[:,1], lrssn_cRate, markersize=6, lab="PD-RSSN (proposed)")
      #
      # scatter!(convfig,[0;prior[:,1]],[NaN;NaN;NaN;prior_cRate], markersize=2, lab="")
      # xlims!( (0,L+min(max_xlim,max(length(r_lrssn),round(length(r_lrcpa)/10)))) )
      # xlabel!(convfig,"Iterations")
      # ylabel!(convfig,"q")
      # ylims!((0,3))

      plot!(rerrorfig,[L;lrcpa[:,1]], [relerror;lrcpa[:,3]], yaxis=:log, lab="lRCPA scale=$(scale)")
      plot!(rerrorfig,[L;lrssn[:,1]], [relerror;lrssn[:,3]], yaxis=:log, lab="PD-RSSN scale=$(scale)")
      plot!(rerrorfig,[0;prior[:,1]], [1;prior[:,3]], yaxis=:log, lab="")
      xlims!( (0,L+min(max_xlim,max(length(r_lrssn),round(length(r_lrcpa)/10)))) )
      xlabel!(rerrorfig,"Iterations")
      ylabel!(rerrorfig,L"\epsilon_{rel}")

   end

   timingfig = plot(scales, timing[2,:,1], lab="PD-RSSN (proposed)")  # was linecolor = ssnColor,
   # timingfig = plot(scales, timing[1,:,1], lab="lRCPA")  # was linecolor = pdhgColor,
   # plot!(timingfig, scales, timing[2,:,1], lab="PD-RSSN (proposed)")  # was linecolor = ssnColor,
   xlabel!(timingfig,"Data scale")
   ylabel!(timingfig,"t (s)")

   iterfig = plot(scales, timing[2,:,2], lab="PD-RSSN (proposed)")  # was linecolor = ssnColor,
   # iterfig = plot(scales, timing[1,:,2], lab="lRCPA")  # was linecolor = pdhgColor,
   # plot!(iterfig, scales, timing[2,:,2], lab="PD-RSSN (proposed)")  # was linecolor = ssnColor,
   xlabel!(iterfig,"Data scale")
   ylabel!(iterfig, "Iterations")


   exp_savetable(table, "result-SPD3"*postfix)

   exp_savefig(costfig, "cost-SPD3"*postfix, true)

   # exp_savefig(convfig, "convrate-SPD3"*postfix, true)

   exp_savefig(rerrorfig, "relerror-SPD3"*postfix, true)

   exp_savefig(timingfig, "timing-SPD3"*postfix, true)

   exp_savefig(iterfig, "iter-SPD3"*postfix, true)

   exp_end()

end
