# IMEX-MRI: Implicit-Explicit Multirate Infinitesimal methods #

#### Rujeko Chinomona and Daniel R. Reynolds ####
Department of Mathematics, Southern Methodist University

This repository contains MATLAB files for test problems highlighted in the paper, [Rujeko Chinomona & Daniel R. Reynolds, "Implicit-explicit multirate infinitesimal methods," arXiv:2007.09776, 2020](https://arxiv.org/abs/2007.09776).

IMEX-MRI methods of orders three (`IMEX-MRI3a`, `IMEX-MRI3b`) and four (`IMEX-MRI4`) are implemented on the additively split multirate problem:

  ```y' = f^E(t,y) + f^I(t,y) + f^F(t,y)```

where ```f^E``` denotes the slow-explicit part, ```f^I``` the slow-implicit part  and ```f^F``` denotes the fast part of the right hand side. In addition, we also run comparison tests with operator splitting methods that are normally used for these types of equations i.e. first order "Lie-Trotter operator splitting" and second order "Strang-Marchuk splitting".

We include two test problems:
 * A modified Kvaerno-Prothero-Robinson (KPR) problem (`driver_prothero_robinson.m`), adapted from [A.Sandu, SIAM J. Numer. Anal., 2019](https://doi.org/10.1137/18M1205492)
 * A stiff brusselator test problem with advection, diffusion and reaction (`driver_brusselatorPDE.m`)

For both tests we output maximum errors, root-mean-square errors and their respective convergence rates. We also output the total number of fast and slow steps taken, the number of slow nonlinear solves, and the number of fast nonlinear solves.  Coefficients for each methods are stored in the files `IMEXMRI3a_high_sym.mat`,  `IMEXMRI3b_high_sym.mat` and  `IMEXMRI4_high_sym.mat`, or can instead be extracted from the text file `imexmricoeffs.txt`. The main implementation routine for the IMEX-MRI methods is `solve_IMEX_MRI.m`. This calls other routines to evolve the inner fast ODE (`solve_IRK.m` or `solve_ERK.m`, depending on the type of fast method) based on explicit or diagonally-implicit Runge--Kutta coefficients stored in `butcher.m`.  Here, we use explicit Runge-Kutta methods for the fast modified ODE on the KPR problem and implicit Runge-Kutta methods on the brusselator test problem, the order of accuracy of each fast RK method is chosen to match the IMEX-MRI method order.
