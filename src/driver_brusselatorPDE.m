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
n     = 11;                              % number of output times
tout  = linspace(t0,tf,n);               % output times
hs    = 0.1*0.5.^(0:10);                 % slow time steps
m     = 5;                              % time scale separation factor

% set reference solution
hmin = 1e-7;
rtol = 2.5e-14;
atol = 1e-14*ones(length(Y0),1);
opts = odeset('RelTol',rtol, 'AbsTol',atol,'InitialStep',hmin);
[t,Ytrue] = ode15s(fn, tout, Y0, opts);
Ytrue = Ytrue';

%------------------------------------------------------------------------------%

filename = 'IMEXMRI3a_high_sym.mat';
mrisolver = 'IMEXMRI3a';
Q = matfile(filename);
G = {double(Q.G0)};
W = {double(Q.W0)};
c = double(Q.c);
innersolver = 'ESDIRK-3-3';
doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)
clear c G W Q

%------------------------------------------------------------------------------%

filename = 'IMEXMRI3b_high_sym.mat';
mrisolver = 'IMEXMRI3b';
Q = matfile(filename);
G = {double(Q.G0)};
W = {double(Q.W0)};
c = double(Q.c);
innersolver = 'ESDIRK-3-3';
doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)
clear c G W Q

%------------------------------------------------------------------------------%
filename = 'IMEXMRI4_high_sym.mat';
mrisolver = 'IMEXMRI4';
Q = matfile(filename);
G = {double(Q.G0),double(Q.G1)};
W = {double(Q.W0),double(Q.W1)};
c = double(Q.c);
innersolver = 'Cash(5,3,4)-SDIRK';
doimexmritest(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,G,W,c,mrisolver,innersolver)
clear c G W Q

------------------------------------------------------------------------------%
% Set up tables for Lie-Trotter splitting
erksolver = 'ERK-1-1';
irksolver = 'SDIRK-2-1';
dofirstordersplit(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,erksolver,irksolver)

%------------------------------------------------------------------------------%
% Set up tables for Strang splitting
erksolver = 'Heun-Euler-ERK';
irksolver = 'SDIRK-2-1';
doStrang(fe,fi,ff,Ji,Jf,tout,Y0,Ytrue,hs,m,erksolver,irksolver)

%------------------------------------------------------------------------------%
mrisolver = 'MRI-GARK-ESDIRK34a';
lam = 0.435866521508458999416019;
c = [0; 1/3; 1/3; 2/3; 2/3; 1; 1];
G{1} = [1/3 0 0 0 0 0 0;
   -lam 0 lam 0 0 0 0
   (3-10*lam)/(24*lam-6), 0, (5-18*lam)/(6-24*lam), 0, 0, 0, 0;
   (-24*lam^2+6*lam+1)/(6-24*lam), 0, (-48*lam^2+12*lam+1)/(24*lam-6), 0, lam, 0 0;
   (3-16*lam)/(12-48*lam), 0, (48*lam^2-21*lam+2)/(12*lam-3), 0, (3-16*lam)/4, 0, 0;
   -lam, 0, 0, 0, 0, 0, lam];
innersolver = 'ESDIRK-3-3';
domrigarktest(fs,ff,Js,Jf,tout,Y0,Ytrue,hs,m,G,c,mrisolver,innersolver)
clear c G

%------------------------------------------------------------------------------%
mrisolver = 'MRI-GARK-ESDIRK46a';
hsmod = 0.1*0.5.^(2:10);
c = [0; 1/5; 1/5; 2/5; 2/5; 3/5; 3/5; 4/5; 4/5; 1; 1];
G0 = zeros(10,11);
G1 = zeros(10,11);
G0(1,1) = 1/5;
G0(2,1) = -1/(4.0);
G0(2,3) = 1/(4.0);
G0(3,1) = (1771023115159.0)/(1929363690800.0);
G0(3,3) = -(1385150376999.0)/(1929363690800.0);
G0(4,1) = (914009.0)/(345800.0);
G0(4,3) = -(1000459.0)/(345800.0);
G0(4,5) = 1/(4.0);
G0(5,1) = (18386293581909.0)/(36657910125200.0);
G0(5,3) = (5506531089.0)/(80566835440.0);
G0(5,5) = -(178423463189.0)/(482340922700.0);
G0(6,1) = (36036097.0)/(8299200.0);
G0(6,3) = (4621.0)/(118560.0);
G0(6,5) = -(38434367.0)/(8299200.0);
G0(6,7) = 1/(4.0);
G0(7,1) = -(247809665162987.0)/(146631640500800.0);
G0(7,3) = (10604946373579.0)/(14663164050080.0);
G0(7,5) = (10838126175385.0)/(5865265620032.0);
G0(7,7) = -(24966656214317.0)/(36657910125200.0);
G0(8,1) = (38519701.0)/(11618880.0);
G0(8,3) = (10517363.0)/(9682400.0);
G0(8,5) = -(23284701.0)/(19364800.0);
G0(8,7) = -(10018609.0)/(2904720.0);
G0(8,9) = 1/(4.0);
G0(9,1) = -(52907807977903.0)/(33838070884800.0);
G0(9,3) = (74846944529257.0)/(73315820250400.0);
G0(9,5) = (365022522318171.0)/(146631640500800.0);
G0(9,7) = -(20513210406809.0)/(109973730375600.0);
G0(9,9) = -(2918009798.0)/(1870301537.0);
G0(10,1) = (19.0)/(100);
G0(10,3) = -(73.0)/(300.0);
G0(10,5) = (127.0)/(300.0);
G0(10,7) = (127.0)/(300.0);
G0(10,9) = -(313.0)/(300.0);
G0(10,11) = 1/(4.0);

G1(3,1) = -(1674554930619.0)/(964681845400.0);
G1(3,3) = (1674554930619.0)/(964681845400.0);
G1(4,1) = -(1007739.0)/(172900.0);
G1(4,3) = (1007739.0)/(172900.0);
G1(5,1) = -(8450070574289.0)/(18328955062600.0);
G1(5,3) = -(39429409169.0)/(40283417720.0);
G1(5,5) = (173621393067.0)/(120585230675.0);
G1(6,1) = -(122894383.0)/(16598400.0);
G1(6,3) = (14501.0)/(237120.0);
G1(6,5) = (121879313.0)/(16598400.0);
G1(7,1) = (32410002731287.0)/(15434909526400.0);
G1(7,3) = -(46499276605921.0)/(29326328100160.0);
G1(7,5) = -(34914135774643.0)/(11730531240064.0);
G1(7,7) = (45128506783177.0)/(18328955062600.0);
G1(8,1) = -(128357303.0)/(23237760.0);
G1(8,3) = -(35433927.0)/(19364800.0);
G1(8,5) = (71038479.0)/(38729600.0);
G1(8,7) = (8015933.0)/(1452360.0);
G1(9,1) = (136721604296777.0)/(67676141769600.0);
G1(9,3) = -(349632444539303.0)/(146631640500800.0);
G1(9,5) = -(1292744859249609.0)/(293263281001600.0);
G1(9,7) = (8356250416309.0)/(54986865187800.0);
G1(9,9) = (17282943803.0)/(3740603074.0);
G1(10,1) = (3.0)/(25.0);
G1(10,3) = -(29.0)/(300.0);
G1(10,5) = (71.0)/(300.0);
G1(10,7) = (71.0)/(300.0);
G1(10,9) = -(149)/(300.0);
G = {G0,G1};
innersolver = 'Cash(5,3,4)-SDIRK';
domrigarktest(fs,ff,Js,Jf,tout,Y0,Ytrue,hsmod,m,G,c,mrisolver,innersolver)
clear c G
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
        [~,Yi,nsi,sni,~] = solve_DIRK(fi,Ji,tspan,Ye(:,end),D,rtol,atol,h,h,h);
        nfi(j) = nfi(j) + nsi;
        snits(j) = snits(j) + sni;

        % solve fast
        [~,Yf,nsf,fni,~] = solve_DIRK(ff,Jf,tspan,Yi(:,end),D,rtol,atol,hf,hf,hf);
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
    fprintf('| %.5e | %5d | %.5e | %.5e | %.5e | %.5e | %.5e |\n',...
            hs(j),m,hf,err_max(j), max_rate(j),err_rms(j),rms_rate(j))

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
    fprintf('| %.5e | %5d | %.5e | %.5e | %.5e | %.5e | %.5e |\n',...
            hs(j),m,hf,err_max(j), max_rate(j),err_rms(j),rms_rate(j))

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
