include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))
include("article-RSSN-comparison-of-algorithms-S2-run.jl")
include("article-RSSN-comparison-of-algorithms-S2-process.jl")

include("article-RSSN-comparison-of-algorithms-SPD3-run.jl")
include("article-RSSN-comparison-of-algorithms-SPD3-process.jl")

function main()

   exp_begin(aux_folders = ["demos","experiments","functions","jltools","plans","report","solvers"],
   prefix = "article-RSSN-comparison-of-algorithms-")
   dbglevel(1)

   dbg(0,"Results folder is " * exp_prefix())

   λ = 1e-6

   # path = "results/article-RSSN-comparison-of-algorithms-2021-06-14-11-03-44/"
   # path2 = "results/report-RSSN-comparison-of-algorithms-2020-08-25-16-32-09/"
   # path3 = "results/report-RSSN-comparison-of-algorithms-2020-08-26-11-14-08/"
   # path4 = "results/report-RSSN-comparison-of-algorithms-2020-08-12-12-01-07/"

   # S2 Problems

   run_RSSN_comparison_algorithms_S2(λ)

   process_RSSN_comparison_algorithms_S2(λ,view_result=false)  #,path=path)

   # P(3)
   run_RSSN_comparison_algorithms_SPD3(λ)

   process_RSSN_comparison_algorithms_SPD3(λ,view_result=false)  #,path=path)


   exp_end()

end
