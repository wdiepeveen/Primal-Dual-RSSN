include(joinpath("..","..","jltools","ExpTools.jl"))
include(joinpath("..","..","jltools","Dbg.jl"))

include("article-RSSN-known-minimizers-S2-run.jl")
include("article-RSSN-known-minimizers-S2-process.jl")

include("article-RSSN-known-minimizers-SPD3-run.jl")
include("article-RSSN-known-minimizers-SPD3-process.jl")

function main()

    exp_begin(aux_folders = ["experiments","demos","functions","jltools","plans","report","solvers"],
    prefix = "article-RSSN-known-minimizers-")
    dbglevel(1)

    dbg(0,"Results folder is " * exp_prefix())

    ℓ = 10
    α = 5. #3.

    # S2 parameters
    σ = τ = 0.5
    max_it = 50
    tol = 1e-10
    ViewResult = true

    relerrors = zeros(2,2)

# S2 problem
        # without dual start

        run_S2_RSSN(ℓ,α,σ,τ,max_it=max_it,tol=tol,dual_prior=false,postfix="cold")

        # with dual start

        run_S2_RSSN(ℓ,α,σ,τ,max_it=max_it,tol=tol,dual_prior=true,postfix="warm")

        relerrors[1,:] = process_S2(view_result=false)#ViewResult)

# P3 problem

# without dual start

        run_SPD3_RSSN(ℓ,α,σ,τ,max_it=max_it,tol=tol,postfix="cold")

# with dual start

        run_SPD3_RSSN(ℓ,α,σ,τ,max_it=max_it,tol=tol,dual_prior=true,postfix="warm")

        relerrors[2,:] = process_SPD3(view_result=false)#ViewResult)

table = raw"\begin{table}[h!]
\centering
\begin{tabular}{c|c|c}
\hline\hline
 &  Cold start  & Warm start  \\ \hline
$\mathcal{M}$ &   $d_{\mathcal{M}}(p^*,\tilde{p})$ & $d_{\mathcal{M}}(p^*,\tilde{p})$ \\ \hline
$S^2$ &"*"$(float32(relerrors[1,1]))&  $(float32(relerrors[1,2])) "* raw"\\ \hline
$\mathcal{P}(3)$ &"*"$(float32(relerrors[2,1]))&  $(float32(relerrors[2,2]))"*
raw"\\ \hline\hline
\end{tabular}
\end{table}"

exp_savetable(table, "result")


end
