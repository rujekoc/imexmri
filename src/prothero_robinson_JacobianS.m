function [J] = prothero_robinson_JacobianS(t,U)
  % [J] = prothero_robinson_JacobianS(t,U)
  %
  % Jacobian of the slow portion of the right hand side


  global alpha
  global omega
  global epsilon
  global lambda_f
  global lambda_s

  dA_dU1 = 1 - (-3 + U(1)*U(1) - cos(omega*t))/2/U(1)/U(1);
  dB_dU2 = 1 - (-2 + U(2)*U(2) - cos(t))/2/U(2)/U(2);
  J = zeros(2,2);
  J(2,1) = -alpha*epsilon*(lambda_f - lambda_s)*dA_dU1;
  J(2,2) = lambda_s*dB_dU2 + sin(t)/2/U(2)/U(2);

end
