function [J] = prothero_robinson_JacobianF(t,U)
% Usage: [J] = prothero_robinson_JacobianF(t,U)
%
% Jacobian of the fast component of the right hand side,
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

  dA_dU1 = 1 - (-3 + U(1)*U(1) - cos(omega*t))/2/U(1)/U(1);
  dB_dU2 = 1 - (-2 + U(2)*U(2) - cos(t))/2/U(2)/U(2);
  J = zeros(2,2);
  J(1,1) = lambda_f*dA_dU1 + omega*sin(omega*t)/2/U(1)/U(1);
  J(1,2) = (1-epsilon)*(lambda_f - lambda_s)/alpha*dB_dU2 ;

end
