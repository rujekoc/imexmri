function [J] = prothero_robinson_JacobianI(t,U)
% Usage: [J] = prothero_robinson_JacobianI(t,U)
%
% Jacobian of the slow-implicit portion of the right hand side,
% Kvaerno-Prothero-Robinson test problem.
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

  dA_dU1 = 1 - (-3 + U(1)*U(1) - cos(omega*t))/2/U(1)/U(1);
  dB_dU2 = 1 - (-2 + U(2)*U(2) - cos(t))/2/U(2)/U(2);
  J = zeros(2,2);
  J(2,1) = -alpha*epsilon*(lambda_f - lambda_s)*dA_dU1;
  J(2,2) = lambda_s*dB_dU2;

end
