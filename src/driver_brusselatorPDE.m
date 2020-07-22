% driver for 1D brusselator test problem (advection-diffusion-reaction system):
%      u' = d1*u_xx + a1*u_x + a - (w+1)*u + u^2*v,
%      v' = d2*v_xx + a2*v_x + w*u - u^2*v,
%      w' = d3*w_xx + a3*w_x + (b-w)/ep - w*u
% with parameters a=0.6, b=2, ep = 1e-2
% diffusion coefficients d1=d2=d3=0.01 and,
% advection coefficients a1=a2=a3=0.001,
% over the spatial interval [0,1] and the time interval [0,10].
% let y = [u;v;w].
% RHS is split into 3 i.e. y' = fe(t,y) + fi(t,y) + ff(t,y)
% fe(t,y) = [a1*u_x; a2*v_x; a3*w_x]
% fi(t,y) = [d1*u_xx; d2*v_xx; d3*w_xx]
% ff(t,y) = [a - (w+1)*u + u^2*v; w*u - u^2*v; (b-w)/ep + w*u]
% Initial conditions are
% u(x,0) = a + 0.1*sin(pi*x),
% v(x,0) = b/a + 0.1*sin(pi*x),
% w(x,0) = b + 0.1*sin(pi*x)
% and stationary boundary conditions. We use a mesh with 100 spatial zones.
%
% A reference solution is obtained using ode15s.
% Test problem is solved with 1st order operator splitting,
% 2nd order Strang splitting, 3rd order IMEXMRI3a and IMEXMRI3b, and 4th order
% IMEXMRI4.
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% July 2020
% All Rights Reserved

% Set up test problem
global Pdata
[Y0,t0,tf] = brusselatorPDE_setup();
fn = @brusselatorPDE;
fe = @brusselatorPDE_explicit;
fi = @brusselatorPDE_implicit;
ff = @brusselatorPDE_fast;
fs = @(t,y) brusselatorPDE_explicit(t,y) + brusselatorPDE_implicit(t,y);
Jn = @brusselatorPDE_Jacobian;
Js = @brusselatorPDE_JacobianS;
Ji = @brusselatorPDE_JacobianI;
Je = @brusselatorPDE_JacobianE;
Jf = @brusselatorPDE_JacobianF;

% time parameters
n     = tf/2+1;                          % number of output times
tout  = linspace(t0,tf,n);               % output times
hs    = (tout(2)-tout(1))*0.5.^(0:10);   % slow time steps
m     = 10;                              % time scale separation factor

% set reference solution
hmin = 1e-7;
hmax = 1.0;
rtol = 2.5e-14;
atol = 1e-14*ones(length(Y0),1);
opts = odeset('RelTol',rtol, 'AbsTol',atol,'InitialStep',hmin, 'MaxStep',hmax);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
Ytrue = Ytrue';

% Visualization of reference solution
nx = Pdata.nx;
xspan = linspace(0,1,nx)';              % spatial discretization

%------------------------------------------------------------------------------%

filename = 'IMEXMRI3a_high_sym.mat';
mrisolver = 'IMEXMRI3a';
Q = matfile(filename);
G = {double(Q.G0)};
W = {double(Q.W0)};
c = double(Q.c);
innersolver = 'ESDIRK-3-3';
doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)

%------------------------------------------------------------------------------%

filename = 'IMEXMRI3b_high_sym.mat';
mrisolver = 'IMEXMRI3b';
Q = matfile(filename);
G = {double(Q.G0)};
W = {double(Q.W0)};
c = double(Q.c);
innersolver = 'ESDIRK-3-3';
doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)

%------------------------------------------------------------------------------%
filename = 'IMEXMRI4_high_sym.mat';
mrisolver = 'IMEXMRI4';
Q = matfile(filename);
G = {double(Q.G0),double(Q.G1)};
W = {double(Q.W0),double(Q.W1)};
c = double(Q.c);
innersolver = 'Cash(5,3,4)-SDIRK';
hs = (tout(2)-tout(1))/32*0.5.^(0:7);
doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)

%------------------------------------------------------------------------------%
% Set up tables for Lie-Trotter splitting
erksolver = 'ERK-1-1';
irksolver = 'IRK-1-1';
hs = (tout(2)-tout(1))*0.5.^(0:15);
dofirstordersplit(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,erksolver,irksolver)

%--------------------------------------------------------------------------------%
% Set up tables for Strang splitting
erksolver = 'Heun-Euler-ERK';
irksolver = 'LobattoIIIA-2-2-IRK';
hs = (tout(2)-tout(1))*0.5.^(0:15);
doStrang(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,erksolver,irksolver)

%------------------------------------------------------------------------------%
function doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)

  % Set innersolver
  Bi = butcher(innersolver);

  % Initialize storage
  numhs = length(hs);
  err_max = zeros(1,numhs);
  err_rms = zeros(1,numhs);
  max_rate = zeros(1,numhs);
  rms_rate = zeros(1,numhs);
  nslow = zeros(1,numhs);
  nfast = zeros(1,numhs);
  snits = zeros(1,numhs);
  fnits = zeros(1,numhs);

  fprintf('Running with %s and innersolver %s\n',mrisolver,innersolver);
  fprintf('|    hs       |   m   |    hf       | max err     | max rate    | rms err     | rms rate    |\n');
  fprintf('|-------------------------------------------------------------------------------------------|\n');

  for j = 1:numhs
    % Set fast time step
    hf = hs(j)/m;

    % Call to solver
    [~,Y,ns,nf,~,~,sni,fni,ierr] = solve_IMEX_MRI(fe,fi,ff,Ji,Jf,tout,Y0,c,G,W,Bi,hs(j),hf);
    if (ierr ~= 0)   % on a solver failure, just skip to the next h
      err_max(j) = 1e4;
      err_rms(j) = 1e4;
      fprintf('    Solver failure, skipping to next h value\n')
      continue
    end

    % Error calculation
    err_max(j) = max(max(abs(Y-Ytrue)));
    err_rms(j) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));

    % Rate of convergence calculation
    if j > 1
      max_rate(j) = log(err_max(j-1)/err_max(j))/log(hs(j-1)/hs(j));
      rms_rate(j) = log(err_rms(j-1)/err_rms(j))/log(hs(j-1)/hs(j));
    end

    % Output error information
    fprintf('| %.5e | %5d | %.5e | %.5e | %.5e | %.5e | %.5e |\n',hs(j),m,hf,err_max(j), max_rate(j),err_rms(j),rms_rate(j))

    % Print out other solver stats
    fprintf('    num slow steps    = %i\n',ns);
    fprintf('    num fast steps    = %i\n',nf);

    if (sni > 0)
      fprintf('    slow Newton iters = %i\n', sni);
      snits(j) = sni;
    end
    if (fni > 0)
      fprintf('    fast Newton iters = %i\n', fni);
      fnits(j) = fni;
    end
  end
  fprintf('|---------------------------------------------------------------------------------------------------------|\n');

  % Best-fit convergence rate
  p1 = polyfit(log(hs),log(err_max),1);
  fprintf('best-fit max rate = %g\n', p1(1));
  p2 = polyfit(log(hs),log(err_rms),1);
  fprintf('best-fit rms rate = %g\n', p2(1));
end

%------------------------------------------------------------------------------%
function dofirstordersplit(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,erksolver,irksolver)

  % Set solvers
  B = butcher(erksolver);
  D = butcher(irksolver);

  % Initialize storage
  numhs = length(hs);
  Y          = zeros(length(Y0),length(tout));
  Y(:,1)     = Y0;
  err_max    = zeros(1,numhs);
  err_rms    = zeros(1,numhs);
  max_rate   = zeros(1,numhs);
  rms_rate   = zeros(1,numhs);
  snits      = zeros(1,numhs);
  fnits      = zeros(1,numhs);
  nfe        = zeros(1,numhs);
  nfi        = zeros(1,numhs);
  nff        = zeros(1,numhs);

  ONEMSM = 1-sqrt(eps);  % coefficient to account for floating-point roundoff

  fprintf('|    hs       |   m   |    hf       | max err     | max rate    | rms err     | rms rate    |\n');
  fprintf('|-------------------------------------------------------------------------------------------|\n');
  for j = 1:numhs
    % Loop over tout values
    for i = 2:length(tout)
      tn = tout(i-1);
      Y0step = Y(:,i-1);
      while (tn < tout(i)*ONEMSM)

        h = min([hs(j), tout(i)-tn]);   % stop at output times
        hf = h/m;                       % fast scale times step
        tspan = [tn,tn+h];              % time span

        % dummy fast stability function, tolerances to enforce fixed timestepping
        estab = @(t,y) h;
        rtol  = 1e20;
        atol  = 1e20;

        % solve slow-explicit
        [~,Ye,nse,~] = solve_ERK(fe,estab,tspan,Y0step,B,rtol,atol,h,h,h);
        nfe(j) = nfe(j) + nse;

        % solve slow-implicit
        [~,Yi,nsi,sni,~] = solve_IRK(fi,Ji,tspan,Ye(:,end),D,rtol,atol,h,h,h);
        nfi(j) = nfi(j) + nsi;
        snits(j) = snits(j) + sni;

        % solve fast
        [~,Yf,nsf,fni,~] = solve_IRK(ff,Jf,tspan,Yi(:,end),D,rtol,atol,hf,hf,hf);
        nff(j) = nff(j) + nsf;
        fnits(j) = fnits(j) + fni;

        tn = tn + h;
        Y0step = Yf(:,end);
      end
      Y(:,i) = Yf(:,end);
    end

    % Error calculation
    err_max(j) = max(max(abs(Y-Ytrue)));
    err_rms(j) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));


    % Rate of convergence calculation
    if j > 1
      max_rate(j) = log(err_max(j-1)/err_max(j))/log(hs(j-1)/hs(j));
      rms_rate(j) = log(err_rms(j-1)/err_rms(j))/log(hs(j-1)/hs(j));
    end

    % Output error information
    fprintf('| %.5e | %5d | %.5e | %.5e | %.5e | %.5e | %.5e |\n',hs(j),m,hf,err_max(j), max_rate(j),err_rms(j),rms_rate(j))

    % Print out other solver stats
    fprintf('    num slow-explicit steps   = %i\n', nfe(j));
    fprintf('    num slow-implicit steps   = %i\n', nfi(j));
    fprintf('    num fast steps            = %i\n', nff(j));
    fprintf('    slow Newton iters         = %i\n', snits(j));
    fprintf('    fast Newton iters         = %i\n', fnits(j));
  end
  % Best-fit convergence rate
  p1 = polyfit(log(hs),log(err_max),1);
  fprintf('best-fit max rate = %g\n', p1(1));
  p2 = polyfit(log(hs),log(err_rms),1);
  fprintf('best-fit rms rate = %g\n', p2(1));
end

%------------------------------------------------------------------------------%

function doStrang(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,erksolver,irksolver)
  % Set solvers
  B = butcher(erksolver);
  D = butcher(irksolver);

  % Initialize storage
  numhs = length(hs);
  Y          = zeros(length(Y0),length(tout));
  Y(:,1)     = Y0;
  err_max    = zeros(1,numhs);
  err_rms    = zeros(1,numhs);
  max_rate   = zeros(1,numhs);
  rms_rate   = zeros(1,numhs);
  snits      = zeros(1,numhs);
  fnits      = zeros(1,numhs);
  nfe        = zeros(1,numhs);
  nfi        = zeros(1,numhs);
  nff        = zeros(1,numhs);

  ONEMSM = 1-sqrt(eps);  % coefficient to account for floating-point roundoff

  fprintf('|    hs       |   m   |    hf       | max err     | max rate    | rms err     | rms rate    |\n');
  fprintf('|-------------------------------------------------------------------------------------------|\n');
  for j = 1:numhs
    % Loop over tout values
    for i = 2:length(tout)
      tn = tout(i-1);
      Y0step = Y(:,i-1);
      while (tn < tout(i)*ONEMSM)

        h = min([hs(j), tout(i)-tn]);   % stop at output times
        hf = h/m;                       % fast scale times step
        tspan = [tn,tn+h];              % time spans
        th1span = [tn,tn + h/2];
        th2span = [tn + h/2,tn+h];

        % dummy fast stability function, tolerances to enforce fixed timestepping
        estab = @(t,y) h;
        rtol  = 1e20;
        atol  = 1e20;

        % solve slow-explicit
        [~,Ye,nse,~] = solve_ERK(fe,estab,th1span,Y0step,B,rtol,atol,h/2,h/2,h/2);
        nfe(j) = nfe(j) + nse;

        % solve slow-implicit
        [~,Yi,nsi,sni,~] = solve_DIRK(fi,Ji,th1span,Ye(:,end),D,rtol,atol,h/2,h/2,h/2);
        nfi(j) = nfi(j) + nsi;
        snits(j) = snits(j) + sni;

        % solve fast implicitly
        [~,Yf,nsf,fni,~] = solve_DIRK(ff,Jf,tspan,Yi(:,end),D,rtol,atol,hf,hf,hf);
        nff(j) = nff(j) + nsf;
        fnits(j) = fnits(j) + fni;

        % solve slow-implicit
        [~,Yi,nsi,fni,~] = solve_DIRK(fi,Ji,th2span,Yf(:,end),D,rtol,atol,h/2,h/2,h/2);
        nfi(j) = nfi(j) + nsi;
        snits(j) = snits(j) + sni;

        % solve slow-explicit
        [~,Ye,nse,~] = solve_ERK(fe,estab,th2span,Yi(:,end),B,rtol,atol,h/2,h/2,h/2);
        nfe(j) = nfe(j) + nse;

        % increment time, reset Y0
        tn = tn + h;
        Y0step = Ye(:,end);
      end
      Y(:,i) = Ye(:,end);
    end

    % Error calculation
    err_max(j) = max(max(abs(Y-Ytrue)));
    err_rms(j) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));

    % Rate of convergence calculation
    if j > 1
      max_rate(j) = log(err_max(j-1)/err_max(j))/log(hs(j-1)/hs(j));
      rms_rate(j) = log(err_rms(j-1)/err_rms(j))/log(hs(j-1)/hs(j));
    end

    % Output error information
    fprintf('| %.5e | %5d | %.5e | %.5e | %.5e | %.5e | %.5e |\n',hs(j),m,hf,err_max(j), max_rate(j),err_rms(j),rms_rate(j))

    % Print out other solver stats
    fprintf('    num slow-explicit steps   = %i\n', nfe(j));
    fprintf('    num slow-implicit steps   = %i\n', nfi(j));
    fprintf('    num fast steps            = %i\n', nff(j));
    fprintf('    slow Newton iters         = %i\n', snits(j));
    fprintf('    fast Newton iters         = %i\n', fnits(j));
  end
  % Best-fit convergence rate
  p1 = polyfit(log(hs),log(err_max),1);
  fprintf('best-fit max rate = %g\n', p1(1));
  p2 = polyfit(log(hs),log(err_rms),1);
  fprintf('best-fit rms rate = %g\n', p2(1));
end
