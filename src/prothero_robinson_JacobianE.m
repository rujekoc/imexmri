function [J] = prothero_robinson_JacobianE(t,U)
  % [J] = prothero_robinson_JacobianE(t,U)
  %
  % Jacobian of the slow-explicit component of the right hand side

  global alpha
  global omega
  global epsilon
  global lambda_f
  global lambda_s
  global C

  J = zeros(2,2);
  J(2,2) = sin(t)/2/U(2)/U(2);
end
