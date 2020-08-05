function [tvals,Y,ns,nf,cnits,clits,snits,fnits,ierr] = solve_IMEX_MRI(fe,fi,ff,Jfi,Jff,tvals,Y0,co,G,W,Bi,hs,hf)
  % usage: [tvals,Y,ns,nf,cnits,clits,snits,fnits,ierr] = solve_IMEX_MRI(fe,fi,ff,Jfi,Jff,tvals,Y0,co,G,W,Bi,hs,hf)
  %
  % Fixed time step IMEX Multirate Infinitesimal Step (IMEX-MRI) solver for the vector-valued ODE problem
  %     y' = fe(t,Y) + fi(t,Y) + ff(t,Y), t >= t0, y in R^n,
  %     Y0 = [y1(t0), y2(t0), ..., yn(t0)]'.
  % This solver is designed to handle a slow method that is any one of:
  % (a) explicit,
  % (b) solve-decoupled implicit (i.e., alternates DIRK at the slow
  %     time scale with subcycling at the fast time scale), and
  % (c) solve-decoupled IMEX (i.e., alternates DIRK at the slow
  %     time scale with subcycling at the fast time scale)
  % The individual time steps are performed using the step_IMEXMRI function;
  % that in turn iterates over the stages for the slow-time-scale method,
  % and calls structure-specific routines for each stage.
  %
  % Inputs:
  %     fe     = function handle for (slow-explicit) ODE RHS
  %     fi     = function handle for (slow-implicit) ODE RHS
  %     ff     = function handle for (fast) ODE RHS
  %     Jfi     = function handle for Jacobian of fi
  %     Jff     = function handle for Jacobian of ff (only required if an
  %              implicit inner method is supplied)
  %     tvals  = array of desired output times, [t0, t1, t2, ..., tN]
  %     Y0     = solution vector at start of step (column vector of length n)
  %     co     = Outer 'slow' method abcissae.  The MRI method requires that
  %              these be sorted, i.e.
  %                    0 <= co(1) <= co(2) <= ... <= co(end) <= 1
  %              and (potentially) padded such that co(end)=1
  %     G      = cell array of MRI "Gamma" matrices, {G0, G1, ...}
  %     W      = cell array of MRI "Omega" matrices, {W0, W1, ...}
  %     Bi     = Butcher table for a single step of an 'inner' (fast) method
  %              (can be ERK, DIRK or IRK)
  %                 Bi = [ci Ai;
  %                       qi bi;
  %                       pi di ]
  %              All components have the same role as with Bi; we
  %              assume that Bi encodes a method with si stages
  %     hs     = step size to use for slow time scale
  %     hf     = desired step size to use for fast time scale,
  %                  hf <= hs*min_{i}(co(i+1)-co(i))
  %              Note: this is only a target step size; in fact we
  %              will determine each substepping interval and find
  %              hinner <= hi such that we take an integer number of
  %              steps to subcycle up to each outer stage time.
  %
  % Outputs:
  %     tvals  = the same as the input array tvals
  %     Y      = [y1(t0+hs), y2(t0+hs), ..., yn(t0+hs)].
  %     ns     = number of 'slow' time steps taken by method
  %     nf     = number of 'fast' time steps taken by method
  %     cnits  = actual number of coupled nonlinear iterations required.
  %     clits  = actual number of coupled linear iterations required.
  %     snits  = actual number of slow-only nonlinear iterations required.
  %     fnits  = actual number of fast-only nonlinear iterations required.
  %     ierr   = flag denoting success (0) or failure (1)
  %
  % Rujeko Chinomona & Daniel R. Reynolds
  % Department of Mathematics
  % Southern Methodist University
  % July 2020
  % All Rights Reserved

  % initialize output arrays
  N = length(tvals)-1;
  n = length(Y0);
  Y = zeros(n,N+1);
  Y(:,1) = Y0;

  % initialize diagnostics and error flag
  ns = 0;
  nf = 0;
  cnits = 0;
  clits = 0;
  snits = 0;
  fnits = 0;
  ierr = 0;

  % set the solver parameters
  ONEMSM = 1-sqrt(eps);  % coefficient to account for floating-point roundoff

  % initialize temporary variables
  t = tvals(1);
  Ynew = Y0;

  % iterate over output time steps
  for tstep = 1:N

    % loop over internal time steps to get to desired output time
    while (t < tvals(tstep+1)*ONEMSM)

      % bound internal time step
      h = min([hs, tvals(tstep+1)-t]);   % stop at output times

      % call IMEX-MRI stepper to do the work, increment counters
      [Ynew,m,cni,cli,sni,fni,ierr] = step_IMEXMRI(fe,fi,ff,Jfi,Jff,t,Ynew,co,G,W,Bi,h,hf);
      ns    = ns + 1;
      nf    = nf + m;
      cnits = cnits + cni;
      clits = clits + cli;
      snits = snits + sni;
      fnits = fnits + fni;

      % handle error flag
      if (ierr ~= 0)
        fprintf('Error: solve_IMEX_MRI encountered an error at t = %g, aborting.\n',t);
        return;
      end

      % update current time
      t = t + h;

    end

    % store updated solution in output array
    Y(:,tstep+1) = Ynew;

  end  % time output loop

  % end function
end


function [Y,m,cnits,clits,snits,fnits,ierr] = step_IMEXMRI(fe,fi,ff,Jfi,Jff,t0,Y0,co,G,W,Bi,hs,hf)
  % usage: [Y,m,cnits,clits,snits,fnits,ierr] = step_IMEXMRI(fe,fi,ff,Jfi,Jff,t0,Y0,co,G,W,Bi,hs,hf)
  %
  % This routine performs a single step of an imex-at-slow multirate
  % infinitesimal step (IMEX-MRI) method for the vector-valued ODE problem
  %     y' = fs(t,Y) + ff(t,Y), t >= t0, y in R^n,
  %     Y0 = [y1(t0), y2(t0), ..., yn(t0)]'.
  % This driver assumes that the IMEX-MRI method has a "solve-decoupled" structure
  % coupling the fast/slow time scales.  Thus at each slow stage one may either solve
  % for a diagonally-implicit stage using a nonlinear root-finding problem, or perform
  % evolution of a modified IVP at the fast time scale.  Logic is included in the stage
  % loop to determine the stage type, and then calls appropriate functions to handle
  % each case.
  %
  % Inputs:
  %     fe     = function handle for (slow-explicit) ODE RHS
  %     fi     = function handle for (slow-implicit) ODE RHS
  %     ff     = function handle for (fast) ODE RHS
  %     Jfi     = function handle for Jacobian of fi
  %     Jff     = function handle for Jacobian of ff (only required if an
  %              implicit inner method is supplied)
  %     t0     = value of independent variable at start of step
  %     Y0     = solution vector at start of step (column vector of length n)
  %     co     = Outer 'slow' method abcissae.  The MRI method requires that
  %              these be sorted, i.e.
  %                    0 <= co(1) <= co(2) <= ... <= co(end) <= 1
  %              and (potentially) padded such that co(end)=1
  %     G      = cell array of MRI "Gamma" matrices, {G0, G1, ...}
  %     W      = cell array of MRI "Omega" matrices, {W0, W1, ...}
  %     Bi     = Butcher table for a single step of an 'inner' (fast) method
  %              (can be ERK, DIRK or IRK)
  %                 Bi = [ci Ai;
  %                       qi bi;
  %                       pi di ]
  %              All components have the same role as with Bi; we
  %              assume that Bi encodes a method with si stages
  %     hs     = step size to use for slow time scale
  %     hf     = desired step size to use for fast time scale,
  %                  hf <= hs*min_{i}(co(i+1)-co(i))
  %              Note: this is only a target step size; in fact we
  %              will determine each substepping interval and find
  %              hinner <= hi such that we take an integer number of
  %              steps to subcycle up to each outer stage time.
  %
  % Outputs:
  %     Y     = [y1(t0+hs), y2(t0+hs), ..., yn(t0+hs)].
  %     m     = actual number of fast substeps used.
  %     cnits = actual number of coupled nonlinear iterations required.
  %     clits = actual number of coupled linear iterations required.
  %     snits = actual number of slow-only nonlinear iterations required.
  %     fnits = actual number of fast-only nonlinear iterations required.
  %     ierr  = flag denoting success (0) or failure (1)

  % check that abcissae satisfy MRI assumptions
  so = length(co);          % number of stages
  Delta_co = zeros(so,1);
  Delta_co(2:end) = co(2:end)-co(1:end-1);
  if ( (abs(co(1)) > 100*eps) || (min(Delta_co) < -100*eps) || (abs(co(so)-1) > 100*eps) )
    error('Error: co does not satisfy MRI assumptions')
  end

  % initialize outputs
  m = 0;
  cnits = 0;
  clits = 0;
  snits = 0;
  fnits = 0;
  n = length(Y0);
  Y = reshape(Y0,n,1);

  % examine G to determine implicitness of each stage
  Gabs = zeros(so,so);
  for k=1:length(G)
    Gabs = Gabs + abs(G{k});
  end
  Gdiag = diag(Gabs);

  % initialize temporary variables
  FE = zeros(n,so);        % storage for all slow-explicit RHS values
  FI = zeros(n,so);        % storage for all slow-implicit RHS values

  % first outer stage is always explicit
  FE(:,1) = fe(t0,Y);
  FI(:,1) = fi(t0,Y);

  % iterate over remaining outer stages
  for stage = 2:length(co)

    % determine which subroutine to call based on stage structure
    %    slow stage type: ERK or DIRK?
    if (Gdiag(stage) > 100*eps)    % DIRK

      % fast time scale evolution?
      if (Delta_co(stage) > 0)  % fast evolution
        fprintf('solve-coupled structure is not supported')
      else                        % no fast evolution
        [Y,FE,FI,sni,ierr] = dirk_nofast(fe,fi,Jfi,Y,t0,hs,FE,FI,co,G,W,stage);
        snits = snits + sni;
      end

    else                           % ERK

      % fast time scale evolution?
      if (Delta_co(stage) > 0)  % fast evolution

        [Y,FE,FI,mi,fni,ierr] = erk_fast(fe,fi,ff,Jff,Y,t0,hs,hf,FE,FI,co,G,W,Bi,stage);
        m = m + mi;
        fnits = fnits + fni;

      else                        % no fast evolution

        [Y,FE,FI,ierr] = erk_nofast(fe,fi,Y,t0,hs,FE,FI,co,G,W,stage);

      end
    end

    % if the stage update encountered an error, return immediately
    if (ierr ~= 0)
      return
    end

  end

  % all stages complete; final stage is already stored in 'Y', and
  % counters are already up-to-date, so just return
end


%------- Utility routines -------%
function [z,FE,FI,nits,ierr] = dirk_nofast(fe,fi,Jfi,z,t,hs,FE,FI,co,G,W,stage)
  % usage: [z,FE,FI,nits,ierr] = dirk_nofast(fe,fi,Jfi,z,t,hs,FE,FI,co,G,W,stage)
  %
  % This routine performs a single "standard" DIRK stage solve.
  %

  % ensure that z is a column vector
  z = reshape(z,length(z),1);

  % determine effective DIRK coefficients
  nGammas = length(G);
  Ai = zeros(1,stage);
  Ci = zeros(1,stage);
  for j=1:stage
    for k=1:nGammas
      Ai(j) = Ai(j) + G{k}(stage,j)/k;
      Ci(j) = Ci(j) + W{k}(stage,j)/k;
    end
  end

  % set Newton solver RHS and Jacobian routines
  tcur = t + co(stage)*hs;
  rhs = z;
  for j = 1:stage-1
    rhs = rhs + hs*Ai(j)*FI(:,j) + hs*Ci(j)*FE(:,j);
  end
  F = @(Z,Fdata) Z - rhs - hs*Ai(stage)*fi(tcur, Z);
  J = @(Z,Fdata) eye(length(z)) - hs*Ai(stage)*Jfi(tcur, Z);

  % call Newton solver to compute new stage solution, update statistics
  [z,nits,ierr] = newton(F, J, z, 1, 1e-12, 1e-12, 20);

  % if Newton method failed, set error flag and return
  if (ierr ~= 0)
    return;
  end

  % construct new stage RHS and return
  FE(:,stage) = fe(tcur,z);
  FI(:,stage) = fi(tcur,z);

  % if solution or RHS include NaNs, set error flag and return
  if ((max(isnan(z))>0) || (max(isnan(FE(:,stage)))>0) || (max(isnan(FI(:,stage)))>0))
    ierr = 1;
    return;
  end

end


function [z,FE,FI,ierr] = erk_nofast(fe,fi,z,t,hs,FE,FI,co,G,W,stage)
  % usage: [z,FE,FI,ierr] = erk_nofast(fe,fi,z,t,hs,FE,FI,co,G,W,stage)
  %
  % This routine performs a single "standard" ERK stage.
  %

  % initialize return flag
  ierr = 0;

  % determine effective ERK coefficients
  nGammas = length(G);
  Ai = zeros(1,stage-1);
  Ci = zeros(1,stage-1);
  for j=1:stage-1
    for k=1:nGammas
      Ai(j) = Ai(j) + G{k}(stage,j)/k;
      Ci(j) = Ci(j) + W{k}(stage,j)/k;
    end
  end

  % ensure that z is a column vector
  z = reshape(z,length(z),1);

  % compute updated stage
  tcur = t + co(stage)*hs;
  for j = 1:stage-1
    z = z + hs*Ai(j)*FI(:,j) + hs*Ci(j)*FE(:,j);
  end

  % construct new stage RHS and return
  FE(:,stage) = fe(tcur,z);
  FI(:,stage) = fi(tcur,z);

  % if solution or RHS include NaNs, set error flag and return
  if ((max(isnan(z))>0) || (max(isnan(FE(:,stage)))>0) || (max(isnan(FI(:,stage)))>0))
    ierr = 1;
    return;
  end

end


function [z,FE,FI,mi,fnits,ierr] = erk_fast(fe,fi,ff,Jff,z,t,hs,hf,FE,FI,co,G,W,Bi,stage)
  % usage: [z,FE,FI,mi,fnits,ierr] = erk_fast(fe,fi,ff,Jff,z,t,hs,hf,FE,FI,co,G,W,Bi,stage)
  %
  % This routine performs a single stage of an MRI method wherein the
  % slow time scale is imex, and the fast time scale requires
  % evolution.
  %

  % initialize return flag
  ierr = 0;

  % dummy fast stability function, tolerances to enforce fixed timestepping
  estab = @(t,y) hs;
  rtol  = 1e20;
  atol  = 1e20;

  % extract stage time intervals
  so = length(co);          % number of stages
  Delta_co = zeros(so,1);
  Delta_co(2:end) = co(2:end)-co(1:end-1);

  % extract RK method information from Bi
  si = size(Bi,2)-1;        % number of stages
  Ai = Bi(1:si,2:si+1);     % RK coefficients
  innerRK = 0;
  if (max(max(abs(triu(Ai)))) > 0)        % implicit components exist
    innerRK = 1;
    if (max(max(abs(triu(Ai,1)))) > 0)   % method is IRK
      innerRK = 2;
    end
  end

  % initialize outputs
  fnits = 0;
  z = reshape(z,length(z),1);   % ensure that z is a column vector


  %--- perform fast evolution ---%

  % construct 'inner' ODE for this stage
  %   RHS function
  fni = RHSFast(stage,t,hs,ff,G,W,Delta_co,FE,FI);
  %   time interval
  tcur = t + co(stage-1)*hs;  % 'current' time (prior to start of this routine)
  tspan = [tcur, tcur + Delta_co(stage)*hs];

  %   num internal time steps
  ni = ceil(Delta_co(stage)*hs/hf);

  %   step size
  hi = Delta_co(stage)*hs/ni;

  % call inner RK method solver to perform substepping
  if (innerRK == 2)         % IRK inner method
    [~, V, mi, nits, ~] = solve_IRK(fni, Jff, tspan, z, Bi, rtol, atol, hi, hi, hi);
    fnits = fnits + nits;
  elseif (innerRK == 1)     % DIRK inner method
    [~, V, mi, nits, ~] = solve_DIRK(fni, Jff, tspan, z, Bi, rtol, atol, hi, hi, hi);
    fnits = fnits + nits;
  else                      % ERK inner method
    [~, V, mi, ~] = solve_ERK(fni, estab, tspan, z, Bi, rtol, atol, hi, hi, hi);
    % [~,V,mi] = solve_ERKfast(fni,tspan,z,Bi,hi);
  end
  z = V(:,end);

  % construct new slow stage RHS and return
  tcur = t + co(stage)*hs;
  FE(:,stage) = fe(tcur,z);
  FI(:,stage) = fi(tcur,z);


  % if solution or RHS include NaNs, set error flag and return
  if ((max(isnan(z))>0) || (max(isnan(FE(:,stage)))>0) || (max(isnan(FI(:,stage)))>0))
    ierr = 1;
    return;
  end

end


function [fi] = RHSFast(stage,tn,H,ff,G,W,DeltaC,FE,FI)
  % usage: [fi] = RHSFast(stage,tn,H,ff,G,W,DeltaC,FE,FI)
  %
  % Inputs:
  %    stage - current slow stage
  %    tn - time at start of current slow step
  %    H - slow step size
  %    ff - function handle for fast ODE RHS
  %    G - cell array of MRI "Gamma" matrices, {G0, G1, ...}
  %    W - cell array of MRI "Gamma" matrices, {W0, W1, ...}
  %    DeltaC - slow abcissae increments
  %    FE - previously-evaluated slow-explicit RHS vectors
  %    FI - previously-evaluated slow-implicit RHS vectors
  %
  % This routine evaluates the modified fast ODE RHS for an MRI
  % method.

  % determine number of slow stages
  so = length(DeltaC);

  % determine number of MRI "Gamma" matrices
  nGammas = length(G);

  % starting 'time' for current stage solve
  tprev = tn + sum(DeltaC(1:stage-1))*H;

  % utility routine to convert tau (MIS-like) to theta (MRI-like)
  theta = @(tau) (tau-tprev)/DeltaC(stage);

  % utility routine to construct column vector with powers of (theta/H)
  thetapow = @(theta) ((theta/H).^(0:nGammas-1))';

  % Gamma matrix for this fast RHS
  Gamma = zeros(so,nGammas);
  Omega = zeros(so,nGammas);
  for col=1:nGammas
    for row=1:stage-1
      Gamma(row,col) = G{col}(stage,row)/DeltaC(stage);
      Omega(row,col) = W{col}(stage,row)/DeltaC(stage);
    end
  end

  % construct fast function handle
  fi = @(tau,v) ff(tau,v) + FI*Gamma*thetapow(theta(tau)) + FE*Omega*thetapow(theta(tau));

end
