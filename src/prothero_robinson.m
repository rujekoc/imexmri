function [Udot] = prothero_robinson(t,U)
% Usage: [Udot] = prothero_robinson(t,U)
%
% Full right hand side for Kvaerno-Prothero-Robinson test problem
%  y^{f}' = \lambda^{f} *A(t,y^{f})
%             + (1-\epsilon)/\alpha * (lambda^{f} - \lambda^{s})* B(t,y^{s})
%             - \omega*sin(\omega*t)/(2y^{f})
%  y^{s}' = -\alpha*\epsilon(lambda^{f} - \lambda^{s})*A(t,y^{f}) + \lambda^{s}*B(t,y^{s})
%           - sin(t)/(2*y^{s})
%
% where  A(t,y^{f}) = (-3+ y^{f}*y^{f} - cos(\omega*t))/(2y^{f})
%        B(t,y^{s}) = (-2+y^{s} * y^{s} - cos(t))/(2y^{s})
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% July 2020

  global omega
  global epsilon
  global C

  A = (-3 + U(1)*U(1) - cos(omega*t))/2/U(1);
  B = (-2 + U(2)*U(2) -cos(t))/2/U(2);

  Udot = C*[A;B] - [omega*sin(omega*t)/2/U(1);sin(t)/2/U(2)];

end
