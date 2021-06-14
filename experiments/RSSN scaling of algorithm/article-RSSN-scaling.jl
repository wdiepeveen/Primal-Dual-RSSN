include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))
include("article-RSSN-scaling-S2-run.jl")
include("article-RSSN-scaling-S2-process.jl")

include("article-RSSN-scaling-SPD3-run.jl")
include("article-RSSN-scaling-SPD3-process.jl")

function main()

   exp_begin(aux_folders = ["demos","experiments","functions","jltools","plans","report","solvers"],
   prefix = "article-RSSN-scaling-of-algorithm-")
   dbglevel(1)

   dbg(0,"Results folder is " * exp_prefix())

   λ = 1e-6

   scales = [5]  # was [1,2,3,4,5]
   S2_scales = [10,15,20,25,30]
   SPD3_scales = [4,7,10,13,16]

   # path = "results/article-RSSN-comparison-of-algorithms-2020-12-23-18-16-52/"
   # path2 = "results/report-RSSN-comparison-of-algorithms-2020-08-25-16-32-09/"
   # path3 = "results/report-RSSN-comparison-of-algorithms-2020-08-26-11-14-08/"
   # path4 = "results/report-RSSN-comparison-of-algorithms-2020-08-12-12-01-07/"

   # S2 Problems

   run_RSSN_scaling_S2(λ,S2_scales)

   process_RSSN_scaling_S2(λ,S2_scales,view_result=false)  #,path=path)


   # P(3)

   run_RSSN_scaling_SPD3(λ,SPD3_scales)

   process_RSSN_scaling_SPD3(λ,SPD3_scales,view_result=false)  #,path=path)


   exp_end()

end
