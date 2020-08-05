function U = prothero_robinson_analy(t)
% Usage: U = prothero_robinson_analy(t)
%
% Analytical solution, Kvaerno-Prothero-Robinson test problem.
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% July 2020
% All Rights Reserved

  global omega
  U = [sqrt(3 + cos(omega*t)); sqrt(2 + cos(t))];

end
