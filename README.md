Primal-Dual Riemannian Semi-smooth Newton
=========================================

Solve manifold-valued image processing problems with a duality-based higher
order method for non-smooth optimisation.

        [1] W. Diepeveen, J. Lellmann.  
        An Inexact Semi-smooth Newton Method on Riemannian Manifolds with Application to Duality-based Total Variation Denoising.
        arXiv preprint arXiv:2102.10309. 2021 Feb 20.

Setup
-----

The recommended (and tested) setup is based on MacOS 11.4 running Julia 1.5.3
from Atom 1.57.0. In that case, the repo can be set up using the following lines
in the REPL window:

    julia> ]
    (v1.5) pkg> activate .
    (env)  pkg> instantiate


Reproducing the experiments in [1]
----------------------------------

The following scripts have been used to produce the results in the
article *Duality-based Higher-order Non-smooth Optimization on Manifolds* (citation [1] above).
After compiling, each script can be executed using the following lines in the REPL window:

        julia> main()

The plots are directly generated after running the script. The manifold-valued images
can be rendered with Asymptote.

* 6.1. Signal with Known Minimizers (Fig. 2 and 3):

        experiments/RSSN known minimizers/article-RSSN-known-minimizers.jl

* 6.2. PD-RSSN for Solving Regularized l^2-TV (Fig. 4 and 5):

        experiments/RSSN comparison of algorithms/article-RSSN-comparison-of-algorithms.jl

* 6.4. Primal-Dual Inexact Riemannian Semi-smooth Newton (PD-IRSSN) (Fig. 6):

        experiments/RSSN inexact Riemannian semismooth Newton/article-RSSN-inexact-semismooth-Newton-S2.jl
