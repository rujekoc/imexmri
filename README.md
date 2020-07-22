# IMEX-MRI: Implicit-Explicit Multirate Infinitesimal methods # 

#### Rujeko Chinomona and Daniel R. Reynolds ####
Department of Mathematics, Southern Methodist University

This repository contains MATLAB files for test problems highlighted in the paper, [Rujeko Chinomona & Daniel R. Reynolds, "Implicit-explicit multirate infinitesimal methods," arXiv:2007.09776, 2020](https://arxiv.org/abs/2007.09776).
   
IMEX-MRI methods of orders three (`IMEX-MRI3a`, `IMEX-MRI3b`) and four (`IMEX-MRI4`) are implemented on the additively split multirate problem: 

  ```y' = f^E(t,y) + f^I(t,y) + f^F(t,y)``` 

where ```f^E``` denotes the slow-explicit part, ```f^I``` the slow-implicit part  and ```f^F``` denotes the fast part of the right hand side. In addition, we also run comparison tests with operator splitting methods that are normally used for these types of equations i.e. first order `Lie-Trotter operator splitting` and second order `Strang splitting`.

We consider two test problems:
 * Kvaerno-Prothero-Robinson (KPR) problem (`driver_prothero_robinson.m`) from [A.Sandu, SIAM J. Numer. Anal., 2019](https://doi.org/10.1137/18M1205492)
 * Brusselator test problem with advection and diffusion (`driver_brusselatorPDE.m`)
 
 For both tests we output maximum errors, root-mean-square errors and their respective convergence rates. We also output number of fast and slow steps taken, number of slow nonlinear solves, and fast nonlinear solves. 
Coefficients for the method are stored in `.mat` files and can also be extracted from (`imexmricoeffs.txt`). The main implementation routine for IMEX-MRI methods is (`solve_IMEX_MRI.m`). This calls other routines to evolve the inner fast ODE (`solve_IRK.m`,`solve_ERK.m`) with coefficients stored in (`butcher.m`). The current implementation uses implicit Runge-Kutta methods of order equal to IMEX-MRI methods as the solver for the fast modified ODE for the brusselator test problem and explicit Runge-Kutta methods for the KPR problem.

