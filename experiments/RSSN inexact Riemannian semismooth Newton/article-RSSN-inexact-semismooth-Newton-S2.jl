include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))

include("article-RSSN-inexact-semismooth-Newton-S2-run.jl")
include("article-RSSN-inexact-semismooth-Newton-S2-process.jl")


function main()

   exp_begin(aux_folders = ["demos","experiments","functions","jltools","plans","report","solvers"],
   prefix = "article-RSSN-inexact-semismooth-Newton-S2-")
   dbglevel(1)

   dbg(0,"Results folder is " * exp_prefix())


   α = 0.5
   δ = 0.1
   ViewResults=false

   # path = "results/report-RSSN-inexact-semismooth-Newton-S2-2020-08-19-13-12-23/"

   run_RSSN_inexact_S2(α,δ)

   process_RSSN_inexact_S2(view_result=ViewResults)#,path=path)

   exp_end()

end
