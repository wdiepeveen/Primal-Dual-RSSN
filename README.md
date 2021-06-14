Primal-Dual Riemannian Semi-smooth Newton
=========================================

Solve manifold-valued image processing problems with a duality-based higher
order method for non-smooth optimisation.

        [1] W. Diepeveen, J. Lellmann.  
        Duality-based Higher-order Non-smooth Optimization on Manifolds.
        arXiv preprint arXiv:2102.10309. 2021 Feb 20.

Setup
-----

The recommended (and tested) setup is based on MacOS 11.4 running Julia 1.5.3
from Atom 1.57.0. In that case, the repo can be setup using the following lines
in the REPL window:

    julia> ]
    (v1.5) pkg> activate .
    (env)  pkg> instantiate


Reproducing the experiments in [1]
----------------------------------

The following options have been used to produce the results in the
article *Duality-based Higher-order Non-smooth Optimization on Manifolds* (citation [1] above).

* 6.1. Signal with Known Minimizers (Fig. 2 and 3):

        (env) $ python -m mflift rcom tv

* 6.2. PD-RSSN for Solving Regularized l^2-TV (Fig. 3, left to right):

        (env) $ python -m mflift flat_1d tv --data-params "labels=(10,2)"
        (env) $ python -m mflift flat_1d tv --data-params "labels=(10,4)"
        (env) $ python -m mflift flat_1d tv --data-params "labels=(2,20)"
        (env) $ python -m mflift flat_1d rof

* Total variation denoising of a curve on the sphere (Fig. 4, left to right):

        (env) $ python -m mflift sphere_1d tv --data-params "dimsubls=2,dimres=10"
        (env) $ python -m mflift sphere_1d tv --data-params "dimsubls=5,dimres=40"
        (env) $ python -m mflift sphere_1d tv --data-params "dimsubls=25,dimres=75"

* Tikhonov denoising of a curve on the Klein bottle (Fig. 10):

        (env) $ python -m mflift klein_1d quadratic --model-params "lbd=50.0"

* Tikhonov inpainting of a 2-d signal of rotations SO(3) (Fig. 11):

        (env) $ python -m mflift cam_2d_inpaint quadratic --data-params "mode='complete'"

* Denoising of surface normals in DEM (Fig. 12, top to bottom):

        (env) $ python -m mflift bull tv --model-params "lbd=0.4"
        (env) $ python -m mflift bull huber --model-params "lbd=0.75,alph=0.1"
        (env) $ python -m mflift bull quadratic --model-params "lbd=3.0"

* Denoising of high resolution cyclic InSAR measurements (Fig. 13):

        (env) $ python -m mflift insar tv --model-params "lbd=0.6"
        (env) $ python -m mflift insar huber --model-params "lbd=0.75,alph=0.1"
        (env) $ python -m mflift insar quadratic --model-params "lbd=1.0"
