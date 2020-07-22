function U = prothero_robinson_analy(t)
  % [U] = prothero_robinson_analy(t,U)
  %
  % Analytical solution


  global omega
  U = [sqrt(3 + cos(omega*t)); sqrt(2 + cos(t))];

end
