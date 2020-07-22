function [Udot] = prothero_robinson_implicit(t,U)
  % [Udot] = prothero_robinson_implicit(t,U)
  %
  % Slow-implicit portion of right hand side

  global alpha
  global omega
  global epsilon
  global lambda_f
  global lambda_s
  global C

  A = (-3 + U(1)*U(1) - cos(omega*t))/2/U(1);
  B = (-2 + U(2)*U(2) -cos(t))/2/U(2);
  Udot = zeros(size(U));
  Udot(2) = -alpha*epsilon*(lambda_f - lambda_s)*A + lambda_s*B;
end
