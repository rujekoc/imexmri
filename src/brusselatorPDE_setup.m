function [Y0,t0,tf] = brusselatorPDE_setup()
% usage: [Y0,t0,tf] = brusselatorPDE_setup()
%
% Parameter and initial condition setup for stiff brusselator PDE test problem.
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% July 2020
% All Rights Reserved

  global Pdata;
  Pdata.a = 0.6;
  Pdata.b = 2;
  Pdata.d1 = 0.01;
  Pdata.d2 = 0.01;
  Pdata.d3 = 0.01;
  Pdata.a1 = 1e-3;
  Pdata.a2 = 1e-3;
  Pdata.a3 = 1e-3;
  Pdata.nx = 100;
  Pdata.dx = 1/(Pdata.nx-1);
  Pdata.ep = 1/100;
  xspan = linspace(0,1,Pdata.nx)';
  t0 = 0;
  tf = 10;

  % Set initial conditions
  u0 = Pdata.a + 0.1*sin(pi*xspan);
  v0 = Pdata.b/Pdata.a + 0.1*sin(pi*xspan);
  w0 = Pdata.b + 0.1*sin(pi*xspan);
  Y0 = [u0; v0; w0];

end
