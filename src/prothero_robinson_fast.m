function [Udot] = prothero_robinson_fast(t,U)
% Usage: [Udot] = prothero_robinson_fast(t,U)
%
% Fast component of the right-hand side,
% Kvaerno-Prothero-Robinson test problem
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% July 2020
% All Rights Reserved

  global alpha
  global omega
  global epsilon
  global lambda_f
  global lambda_s

  A = (-3 + U(1)*U(1) - cos(omega*t))/2/U(1);
  B = (-2 + U(2)*U(2) -cos(t))/2/U(2);
  Udot = zeros(size(U));
  Udot(1) = lambda_f*A + (1-epsilon)/alpha * (lambda_f - lambda_s)* B - ...
  omega*sin(omega*t)/2/U(1);

end
