function [Udot] = prothero_robinson(t,U)
  % [Udot] = prothero_robinson(t,U)
  %
  % Full right handside of prothero Robinson
  % y^{f}' = \lambda^{f} *A(t,y^{f})
  %             + (1-\epsilon)/\alpha * (lambda^{f} - \lambda^{s})* B(t,y^{s})
  %             - \omega*sin(\omega*t)/(2y^{f})
  %  y^{s}' = -\alpha*\epsilon(lambda^{f} - \lambda^{s})*A(t,y^{f}) + \lambda^{s}*B(t,y^{s})
  %           - sin(t)/(2*y^{s})
  %
  % where  A(t,y^{f}) = (-3+ y^{f}*y^{f} - cos(\omega*t))/(2y^{f})
  %        B(t,y^{s}) = (-2+y^{s} * y^{s} - cos(t))/(2y^{s})
  %
  % Rujeko Chinomona
  % July 2020

  global alpha
  global omega
  global epsilon
  global lambda_f
  global lambda_s
  global C
  A = (-3 + U(1)*U(1) - cos(omega*t))/2/U(1);
  B = (-2 + U(2)*U(2) -cos(t))/2/U(2);

  % Udot(1) = lambda_f*A + (1-epsilon)/alpha * (lambda_f - lambda_s)* B - ...
  % omega*sin(omega*t)/2/U(1);
  %
  % Udot(2) = -alpha*epsilon(lambda_f - lambda_s)*A + lambda_s*B - ...
  % sin(t)/2/U(2);

  Udot = C*[A;B] - [omega*sin(omega*t)/2/U(1);sin(t)/2/U(2)];

end
