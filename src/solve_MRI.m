function [tvals,Y,ns,nf,cnits,clits,snits,fnits,ierr] = solve_MRI(fs,ff,Js,Jf,tvals,Y0,co,G,Bi,hs,hf)
% usage: [tvals,Y,ns,nf,cnits,clits,snits,fnits,ierr] = solve_MRI(fs,ff,Js,Jf,tvals,Y0,co,G,Bi,hs,hf)
%
% Fixed time step Multirate Infinitesimal Step (MRI) solver for the vector-valued ODE problem
%     y' = fs(t,Y) + ff(t,Y), t >= t0, y in R^n,
%     Y0 = [y1(t0), y2(t0), ..., yn(t0)]'.
% This solver is designed to handle a slow method that is any one of:
% (a) explicit,
% (b) solve-decoupled implicit (i.e., alternates DIRK at the slow
%     time scale with subcycling at the fast time scale), and
% (c) solve-coupled implicit (i.e., solves a fully-coupled system
%     for _both_ DIRK at the slow time scale _and_ subcycling at
%     the fast time scale.
% The individual time steps are performed using the step_MRI function;
% that in turn iterates over the stages for the slow-time-scale method,
% and calls structure-specific routines for each stage.
%
% Inputs:
%     fs     = function handle for (slow) ODE RHS
%     ff     = function handle for (fast) ODE RHS
%     Js     = function handle for Jacobian of fs
%     Jf     = function handle for Jacobian of ff (only required if an
%              implicit inner method is supplied)
%     tvals  = array of desired output times, [t0, t1, t2, ..., tN]
%     Y0     = solution vector at start of step (column vector of length n)
%     co     = Outer 'slow' method abcissae.  The MRI method requires that
%              these be sorted, i.e.
%                    0 <= co(1) <= co(2) <= ... <= co(end) <= 1
%              and (potentially) padded such that co(end)=1
%     G      = cell array of MRI "Gamma" matrices, {G0, G1, ...}
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
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% October 2019
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

      % call MIS stepper to do the work, increment counters
      [Ynew,m,cni,cli,sni,fni,ierr] = step_MRI(fs,ff,Js,Jf,t,Ynew,co,G,Bi,h,hf);
      ns    = ns + 1;
      nf    = nf + m;
      cnits = cnits + cni;
      clits = clits + cli;
      snits = snits + sni;
      fnits = fnits + fni;

      % handle error flag
      if (ierr ~= 0)
         fprintf('Error: solve_MRI encountered an error at t = %g, aborting.\n',t);
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




function [Y,m,cnits,clits,snits,fnits,ierr] = step_MRI(fs,ff,Js,Jf,t0,Y0,co,G,Bi,hs,hf)
% usage: [Y,m,cnits,clits,snits,fnits,ierr] = step_MRI(fs,ff,Js,Jf,t0,Y0,co,G,Bi,hs,hf)
%
% This routine performs a single step of an implicit-at-slow multirate
% infinitesimal step (MIS) method for the vector-valued ODE problem
%     y' = fs(t,Y) + ff(t,Y), t >= t0, y in R^n,
%     Y0 = [y1(t0), y2(t0), ..., yn(t0)]'.
% This driver assumes that the 'outer' Butcher table is
% diagonally-implicit, and thus each slow MIS stage must be
% computed through solution of a nonlinear root-finding problem.
%
% This driver can handle both "solve-decoupled" and fully-coupled
% approaches for the fast/slow time scales, calling appropriate functions
% to handle each case.
%
% Inputs:
%     fs     = function handle for (slow) ODE RHS
%     ff     = function handle for (fast) ODE RHS
%     Js     = function handle for Jacobian of fs
%     Jf     = function handle for Jacobian of ff (only required if an
%              implicit inner method is supplied)
%     t0     = value of independent variable at start of step
%     Y0     = solution vector at start of step (column vector of length n)
%     co     = Outer 'slow' method abcissae.  The MRI method requires that
%              these be sorted, i.e.
%                    0 <= co(1) <= co(2) <= ... <= co(end) <= 1
%              and (potentially) padded such that co(end)=1
%     G      = cell array of MRI "Gamma" matrices, {G0, G1, ...}
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
Delta_co = co(2:end)-co(1:end-1);
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
   Gabs(2:end,:) = Gabs(2:end,:) + abs(G{k});
end
Gdiag = diag(Gabs);

% initialize temporary variables
Fs = zeros(n,length(co));        % storage for all slow RHS values

% first outer stage is always explicit
Fs(:,1) = fs(t0,Y);

% iterate over remaining outer stages
for stage = 2:length(co)

   % determine which subroutine to call based on stage structure
   %    slow stage type: ERK or DIRK?
   if (Gdiag(stage) > 100*eps)    % DIRK

      % fast time scale evolution?
      if (Delta_co(stage-1) > 0)  % fast evolution

         [Y,Fs,mi,cni,cli,fni,ierr] = dirk_fast(fs,ff,Jf,Y,t0,hs,hf,Fs,co,G,Bi,stage);
         m = m + mi;
         cnits = cnits + cni;
         clits = clits + cli;
         fnits = fnits + fni;

      else                        % no fast evolution

         [Y,Fs,sni,ierr] = dirk_nofast(fs,Js,Y,t0,hs,Fs,co,G,stage);
         snits = snits + sni;

      end

   else                           % ERK

      % fast time scale evolution?
      if (Delta_co(stage-1) > 0)  % fast evolution

         [Y,Fs,mi,fni,ierr] = erk_fast(fs,ff,Jf,Y,t0,hs,hf,Fs,co,G,Bi,stage);
         m = m + mi;
         fnits = fnits + fni;

      else                        % no fast evolution

         [Y,Fs,ierr] = erk_nofast(fs,Y,t0,hs,Fs,co,G,stage);

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



function [z,Fs,nits,ierr] = dirk_nofast(fs,Js,z,t,hs,Fs,co,G,stage)
% usage: [z,Fs,nits,ierr] = dirk_nofast(fs,Js,z,t,hs,Fs,co,G,stage)
%
% This routine performs a single "standard" DIRK stage solve.
%
% Input/Output argument names match step_MRI.

% ensure that z is a column vector
z = reshape(z,length(z),1);

% determine effective DIRK coefficients
nGammas = length(G);
Ai = zeros(1,stage);
for j=1:stage
   for k=1:nGammas
      Ai(j) = Ai(j) + G{k}(stage-1,j)/k;
   end
end

% set Newton solver RHS and Jacobian routines
tcur = t + co(stage)*hs;
rhs = z;
for j = 1:stage-1
  rhs = rhs + hs*Ai(j)*Fs(:,j);
end
F = @(z,Fdata) z - rhs - hs*Ai(stage)*fs(tcur, z);
J = @(z,Fdata) eye(length(z)) - hs*Ai(stage)*Js(tcur, z);

% call Newton solver to compute new stage solution, update statistics
[z,nits,ierr] = newton(F, J, z, 1, 1e-10, 1e-10, 20);

% if Newton method failed, set error flag and return
if (ierr ~= 0)
  return;
end

% construct new stage RHS and return
Fs(:,stage) = fs(tcur,z);

% if solution or RHS include NaNs, set error flag and return
if ((max(isnan(z))>0) || (max(isnan(Fs(:,stage)))>0))
  ierr = 1;
  return;
end

end




function [z,Fs,ierr] = erk_nofast(fs,z,t,hs,Fs,co,G,stage)
% usage: [z,Fs,ierr] = erk_nofast(fs,z,t,hs,Fs,co,G,stage)
%
% This routine performs a single "standard" ERK stage.
%
% Input/Output argument names match step_MRI.

% initialize return flag
ierr = 0;

% determine effective ERK coefficients
nGammas = length(G);
Ai = zeros(1,stage-1);
for j=1:stage-1
   for k=1:nGammas
      Ai(j) = Ai(j) + G{k}(stage-1,j)/k;
   end
end

% ensure that z is a column vector
z = reshape(z,length(z),1);

% compute updated stage
tcur = t + co(stage)*hs;
for j = 1:stage-1
  z = z + hs*Ai(j)*Fs(:,j);
end

% construct new stage RHS and return
Fs(:,stage) = fs(tcur,z);

% if solution or RHS include NaNs, set error flag and return
if ((max(isnan(z))>0) || (max(isnan(Fs(:,stage)))>0))
  ierr = 1;
  return;
end

end




function [z,Fs,mi,fnits,ierr] = erk_fast(fs,ff,Jf,z,t,hs,hf,Fs,co,G,Bi,stage)
% usage: [z,Fs,mi,fnits,ierr] = erk_fast(fs,ff,Jf,z,t,hs,hf,Fs,co,G,Bi,stage)
%
% This routine performs a single stage of an MIS method wherein the
% slow time scale is explicit, and the fast time scale requires
% evolution.
%
% Input/Output argument names match step_MRI.

% initialize return flag
ierr = 0;

% dummy fast stability function, tolerances to enforce fixed timestepping
estab = @(t,y) hs;
rtol  = 1e20;
atol  = 1e20;

% extract stage time intervals
Delta_co = co(2:end)-co(1:end-1);

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
fi = RHSFast(stage,t,hs,ff,G,Delta_co,Fs);
%   time interval
tcur = t + co(stage-1)*hs;  % 'current' time (prior to start of this routine)
tspan = [tcur, tcur + Delta_co(stage-1)*hs];
%   num internal time steps
ni = ceil(Delta_co(stage-1)*hs/hf);
%   step size
hi = Delta_co(stage-1)*hs/ni;

% call inner RK method solver to perform substepping
if (innerRK == 2)         % IRK inner method
  [~, V, mi, nits, ~] = solve_IRK(fi, Jf, tspan, z, Bi, rtol, atol, hi, hi, hi);
  fnits = fnits + nits;
elseif (innerRK == 1)     % DIRK inner method
  [~, V, mi, nits, ~] = solve_DIRK(fi, Jf, tspan, z, Bi, rtol, atol, hi, hi, hi);
  fnits = fnits + nits;
else                      % ERK inner method
  [~, V, mi, ~] = solve_ERK(fi, estab, tspan, z, Bi, rtol, atol, hi, hi, hi);
end
z = V(:,end);


% construct new slow stage RHS and return
tcur = t + co(stage)*hs;
Fs(:,stage) = fs(tcur,z);


% if solution or RHS include NaNs, set error flag and return
if ((max(isnan(z))>0) || (max(isnan(Fs(:,stage)))>0))
  ierr = 1;
  return;
end

end




function [z,Fs,mi,cnits,clits,fnits,ierr] = dirk_fast(fs,ff,Jf,z,t,hs,hf,Fs,co,G,Bi,stage)
% usage: [z,Fs,mi,cnits,clits,fnits,ierr] = dirk_fast(fs,ff,Jf,z,t,hs,hf,Fs,co,G,Bi,stage)
%
% This routine performs a single stage of a "coupled" MIS method, i.e.,
% a stage that is DIRK on the slow time scale _and_ involves
% evolution of the fast time scale, in a fully-coupled fashion.
%
% Here, that problem is solved using a 'brute-force', unpreconditioned inexact
% Newton-Krylov solver, that uses finite-difference approximations for each
% Jacobian-vector product.
%
% Input/Output argument names match step_MRI.

% initialize outputs
z = reshape(z,length(z),1);   % ensure that z is a column vector

% extract stage time intervals
Delta_co = co(2:end)-co(1:end-1);

% extract RK method information from Bi
[~, Bcols] = size(Bi);
si = Bcols - 1;           % number of stages
Ai = Bi(1:si,2:si+1);     % RK coefficients
innerRK = 0;
if (max(max(abs(triu(Ai)))) > 0)        % implicit components exist
   innerRK = 1;
   if (max(max(abs(triu(Ai,1)))) > 0)   % method is IRK
      innerRK = 2;
   end
end

% set up reusable data for Newton solver residual routine
Fdat.zold = z;
Fdat.Fs = Fs;
Fdat.stage = stage;
Fdat.G = G;
Fdat.Jf = Jf;
Fdat.Bi = Bi;
Fdat.tn = t;
Fdat.hs = hs;
Fdat.Delta_co = Delta_co;
Fdat.tspan = [t + co(stage-1)*hs, t + co(stage)*hs];
Fdat.hi = Delta_co(stage-1)*hs / ceil(Delta_co(stage-1)*hs/hf);
Fdat.tcur = t + co(stage)*hs;
Fdat.fs = fs;
Fdat.ff = ff;
Fdat.innerRK = innerRK;
Fdat.fnits = 0;

% solve for updated state solution
global mi;
mi = 0;
Nmaxit = 20;
Ntol = 1e-10;
eta = 1e-3;
Kmaxit = min(100,length(z));
P = @(x,Pdat) x;    % no preconditioning
[z,~,cnits,clits,ierr] = newton_gmres(@Fcoupled,Fdat,z,Nmaxit,Ntol,eta,Kmaxit,P,0,0);

% extract 'fast' statistics and updated Fs storage, and return
fnits = Fdat.fnits;
Fs(:,stage) = fs(Fdat.tcur, z);

% if solution or RHS include NaNs, set error flag and return
if ((max(isnan(z))>0) || (max(isnan(Fs(:,stage)))>0))
  ierr = 1;
  return;
end

end




function [F] = Fcoupled(z,Fdat)
% usage: [F] = Fcoupled(z,Fdat)
%
% This routine evalutes the nonlinear residual for the coupled ImpMIS stage solve.

% access global fast time step counter
global mi;

% fill 'current' column of Fs with slow RHS evaluated at guess
Fdat.Fs(:,Fdat.stage) = Fdat.fs(Fdat.tcur, z);

% update fast RHS with current Fs values
fi = RHSFast(Fdat.stage,Fdat.tn,Fdat.hs,Fdat.ff,Fdat.G,Fdat.Delta_co,Fdat.Fs);

% dummy fast stability function, tolerances to enforce fixed timestepping
estab = @(t,y) Fdat.hs;
rtol  = 1e20;
atol  = 1e20;

% call inner RK method solver, and update work counters
if (Fdat.innerRK == 2)         % IRK inner method
  [~, V, miloc, nits, ~] = solve_IRK(fi, Fdat.Jf, Fdat.tspan, Fdat.zold, Fdat.Bi, rtol, atol, Fdat.hi, Fdat.hi, Fdat.hi);
  Fdat.fnits = Fdat.fnits + nits;
elseif (Fdat.innerRK == 1)     % DIRK inner method
  [~, V, miloc, nits, ~] = solve_DIRK(fi, Fdat.Jf, Fdat.tspan, Fdat.zold, Fdat.Bi, rtol, atol, Fdat.hi, Fdat.hi, Fdat.hi);
  Fdat.fnits = Fdat.fnits + nits;
else                      % ERK inner method
  [~, V, miloc, ~] = solve_ERK(fi, estab, Fdat.tspan, Fdat.zold, Fdat.Bi, rtol, atol, Fdat.hi, Fdat.hi, Fdat.hi);
end
mi = mi + miloc;

% compute nonlinear residual
F = z - V(:,end);

end




function [fi] = RHSFast(stage,tn,H,ff,G,DeltaC,Fs)
% usage: [fi] = RHSFast(stage,tn,H,ff,G,DeltaC,Fs)
%
% Inputs:
%    stage - current slow stage
%    tn - time at start of current slow step
%    H - slow step size
%    ff - function handle for fast ODE RHS
%    G - cell array of MRI "Gamma" matrices, {G0, G1, ...}
%    DeltaC - slow abcissae increments
%    Fs - previously-evaluated slow RHS vectors
%
% This routine evaluates the modified fast ODE RHS for an MRI
% method.

% determine number of slow stages
so = length(DeltaC)+1;

% determine number of MRI "Gamma" matrices
nGammas = length(G);

% starting 'time' for current stage solve
tprev = tn + sum(DeltaC(1:stage-2))*H;

% utility routine to convert tau (MIS-like) to theta (MRI-like)
theta = @(tau) (tau-tprev)/DeltaC(stage-1);

% utility routine to construct column vector with powers of (theta/H)
thetapow = @(theta) ((theta/H).^(0:nGammas-1))';

% Gamma matrix for this fast RHS
Gamma = zeros(so,nGammas);
for col=1:nGammas
   for row=1:stage
      Gamma(row,col) = G{col}(stage-1,row)/DeltaC(stage-1);
   end
end

% construct fast function handle
fi = @(tau,v) ff(tau,v) + Fs*Gamma*thetapow(theta(tau));

end
