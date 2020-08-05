% driver for Kvaerno Prothero Robinson problem from
% Sandu 2019: A Class of Multirate Infinitesimal GARK methods
% Test problem is nonlinear and given by
%
%  y^{f}' = \lambda^{f} *A(t,y^{f})
%             + (1-\epsilon)/\alpha * (lambda^{f} - \lambda^{s})* B(t,y^{s})
%             - \omega*sin(\omega*t)/(2y^{f})
%  y^{s}' = -\alpha*\epsilon(lambda^{f} - \lambda^{s})*A(t,y^{f}) + \lambda^{s}*B(t,y^{s})
%           - sin(t)/(2*y^{s})
%
% where  A(t,y^{f}) = (-3 + y^{f}*y^{f} - cos(\omega*t))/(2y^{f})
%        B(t,y^{s}) = (-2 + y^{s}*y^{s} - cos(t))/(2y^{s})
%
% alpha = 1, epsilon = 0.1, omega = 20
% \lambda^{f} = -10, \lambda^{s} = -1
%
% Analytical solution
% y^{f}(t) = \sqrt(2+cos(\omega*t)) ,   y^{s} = \sqrt(2 + cos(t))
%
% Time interval [0,5\pi/2]
%
% Test problem is solved with 1st order operator splitting,
% 2nd order Strang splitting, 3rd order IMEXMRI3a and IMEXMRI3b, and 4th order
% IMEXMRI4.
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% July 2020
% All Rights Reserved

% Set up problem
[Y0,t0,tf] = prothero_robinson_setup();
fn = @prothero_robinson;
fe = @prothero_robinson_explicit;
fi = @prothero_robinson_implicit;
ff = @prothero_robinson_fast;
fs = @(t,y) prothero_robinson_explicit(t,y) + prothero_robinson_implicit(t,y);
Jn = @prothero_robinson_Jacobian;
Ji = @prothero_robinson_JacobianI;
Je = @prothero_robinson_JacobianE;
Jf = @prothero_robinson_JacobianF;
Js = @prothero_robinson_JacobianS;

% time parameters
n       = 21;                                    % number of output times
tout    = linspace(t0,tf,n);                     % output times
hmax    = (tf-t0)/(n-1);                         % slow time steps
hs      = hmax*0.5.^(0:7);
m       = 20;                                    % times scale separation factor

% Compute analytic solution
Ytrue = zeros(2,n);
for i = 1:n
  Ytrue(:,i) = prothero_robinson_analy(tout(i));
end

% %------------------------------------------------------------------------------%

filename = 'IMEXMRI3a_high_sym.mat';
mrisolver = 'IMEXMRI3a';
Q = matfile(filename);
G = {double(Q.G0)};
W = {double(Q.W0)};
c = double(Q.c);
innersolver = 'ERK-3-3';
doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)

%------------------------------------------------------------------------------%

filename = 'IMEXMRI3b_high_sym.mat';
mrisolver = 'IMEXMRI3b';
Q = matfile(filename);
G = {double(Q.G0)};
W = {double(Q.W0)};
c = double(Q.c);
innersolver = 'ERK-3-3';
doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)

%------------------------------------------------------------------------------%
filename = 'IMEXMRI4_high_sym.mat';
mrisolver = 'IMEXMRI4';
Q = matfile(filename);
G = {double(Q.G0),double(Q.G1)};
W = {double(Q.W0),double(Q.W1)};
c = double(Q.c);
innersolver = 'ERK-4-4';
doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)

%------------------------------------------------------------------------------%
% Set up tables for Lie-Trotter splitting
erksolver = 'ERK-1-1';
irksolver = 'IRK-1-1';
hs = hmax*0.5.^(0:10);
dofirstordersplit(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,erksolver,irksolver)

%--------------------------------------------------------------------------------%
% Set up tables for Strang splitting
erksolver = 'Heun-Euler-ERK';
irksolver = 'LobattoIIIA-2-2-IRK';
hs = hmax*0.5.^(0:10);
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
    fprintf('| %.5e | %5d | %.5e | %.5e | %.5e | %.5e | %.5e |\n',...
            hs(j),m,hf,err_max(j), max_rate(j),err_rms(j),rms_rate(j))

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

        % solve fast explicitly
        [~,Yf,nsf,~] = solve_ERK(ff,estab,tspan,Yi(:,end),B,rtol,atol,hf,hf,hf);
        nff(j) = nff(j) + nsf;

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
    fprintf('| %.5e | %5d | %.5e | %.5e | %.5e | %.5e | %.5e |\n',...
            hs(j),m,hf,err_max(j), max_rate(j),err_rms(j),rms_rate(j))

    % Print out other solver stats
    fprintf('    num slow-explicit steps   = %i\n', nfe(j));
    fprintf('    num slow-implicit steps   = %i\n', nfi(j));
    fprintf('    num fast steps            = %i\n', nff(j));
    fprintf('    slow Newton iters         = %i\n', snits(j));
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
        [~,Yf,nsf,~] = solve_ERK(ff,estab,tspan,Yi(:,end),B,rtol,atol,hf,hf,hf);
        nff(j) = nff(j) + nsf;

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
    fprintf('| %.5e | %5d | %.5e | %.5e | %.5e | %.5e | %.5e |\n',...
            hs(j),m,hf,err_max(j), max_rate(j),err_rms(j),rms_rate(j))

    % Print out other solver stats
    fprintf('    num slow-explicit steps   = %i\n', nfe(j));
    fprintf('    num slow-implicit steps   = %i\n', nfi(j));
    fprintf('    num fast steps            = %i\n', nff(j));
    fprintf('    slow Newton iters         = %i\n', snits(j));
  end
  % Best-fit convergence rate
  p1 = polyfit(log(hs),log(err_max),1);
  fprintf('best-fit max rate = %g\n', p1(1));
  p2 = polyfit(log(hs),log(err_rms),1);
  fprintf('best-fit rms rate = %g\n', p2(1));
end
