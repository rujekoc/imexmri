function [U0,t0,tf] = prothero_robinson_setup()
  % Usage: [U0,t0,tf] = prothero_robinson_setup()
  %
  % Output:
  % U0 - initial conditions
  % [t0,tf] - time interval


  % define reusable global variables
  global alpha
  global omega
  global epsilon
  global lambda_f
  global lambda_s
  global C


  % set problem parameters
  alpha = 1;
  omega = 20;
  epsilon = 0.1;
  lambda_f = -10;
  lambda_s = -1;

  % set time interval
  t0 = 0;
  tf = 5*pi/2;

  % set initial condition
  U0 = prothero_robinson_analy(0);

end
