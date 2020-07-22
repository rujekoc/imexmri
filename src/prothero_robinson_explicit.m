function [Udot] = prothero_robinson_explicit(t,U)
  % [Udot] = prothero_robinson_explicit(t,U)
  %
  % Slow-explicit portion of the right hand side


  global alpha
  global omega
  global epsilon
  global lambda_f
  global lambda_s
  global C

  A = (-3 + U(1)*U(1) - cos(omega*t))/2/U(1);
  B = (-2 + U(2)*U(2) -cos(t))/2/U(2);
  Udot = zeros(size(U));
  Udot(2) = -sin(t)/2/U(2);

end