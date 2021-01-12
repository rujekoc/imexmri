function domrigarktest(fs,ff,Js,Jf,tout,Y0,Ytrue,hs,m,G,c,mrisolver,innersolver)
  % USAGE: domrigarktest(fs,ff,Js,Jf,tout,Y0,Ytrue,hs,m,G,c,mrisolver,innersolver)
  %
  % INPUTS:
  % fs    - slow rhs
  % Js    - Jacobian of slow
  % Jf    - Jacobian of fast
  % tout  - output times
  % Y0    - initial conditions
  % Ytrue - reference solution
  % hs    - slow time step
  % m     - time scale separation factor hf = hs/m
  % G     - coupling coefficients for slow
  % c     - abscissae for method
  % mrisolver - name of IMEX-MRI-GARK method
  % innersolver - name of fast integrator
  %
  % Prints out different stats: max error, rms error, number of slow steps
  % number of fast steps, number of slow Newton iterations, number of fast Newton
  % iterations
  %
  % Rujeko Chinomona
  % Southern Methodist University
  % January 2021

  % Set innersolver
  Bi = butcher(innersolver);

  % Initialize storage
  numhs = length(hs);
  err_max = zeros(1,numhs);
  err_rms = zeros(1,numhs);
  max_rate = zeros(1,numhs);
  rms_rate = zeros(1,numhs);
  nslow = zeros(1,numhs);
  nfast = zeros(1,numhs);
  snits = zeros(1,numhs);
  fnits = zeros(1,numhs);

  fprintf('Running with %s and innersolver %s\n',mrisolver,innersolver);
  fprintf('|    hs       |   m   |    hf       | max err     | max rate    | rms err     | rms rate    | time        |\n');
  fprintf('|---------------------------------------------------------------------------------------------------------|\n');
  for j = 1:numhs
    % Set fast time step
    hf = hs(j)/m;

    % Call to solver
    [tvals,Y,ns,nf,cni,cli,sni,fni,ierr] = solve_MRI(fs,ff,Js,Jf,tout,Y0,c,G,Bi,hs(j),hf);
    if (ierr ~= 0)   % on a solver failure, just skip to the next h
      err_max(k,j) = 1e4;
      err_rms(k,j) = 1e4;
      fprintf('    Solver failure, skipping to next h value\n')
      continue
    end

    % Error calculation
    err_max(j) = max(max(abs(Y-Ytrue)));
    err_rms(j) = sqrt(sum(sum((Y-Ytrue).^2))/numel(Y));

    % Rate of convergence calculation
    if j > 1
      max_rate(j) = log(err_max(j-1)/err_max(j))/log(hs(j-1)/hs(j));
      rms_rate(j) = log(err_rms(j-1)/err_rms(j))/log(hs(j-1)/hs(j));
    end

    % Output error information
    fprintf('| %.5e | %5d | %.5e | %.5e | %.5e | %.5e | %.5e |\n',...
            hs(j),m,hf,err_max(j), max_rate(j),err_rms(j),rms_rate(j))

    % Print out other solver stats
    fprintf('    num slow steps    = %i\n',ns);
    fprintf('    num fast steps    = %i\n',nf);

    if (sni > 0)
      fprintf('    slow Newton iters = %i\n', sni);
      snits(j) = sni;
    end
    if (fni > 0)
      fprintf('    fast Newton iters = %i\n', fni);
      fnits(j) = fni;
    end

  end
  fprintf('|---------------------------------------------------------------------------------------------------------|\n');
  % Best-fit convergence rate
  p1 = polyfit(log(hs),log(err_max),1);
  fprintf('best-fit max rate = %g\n', p1(1));
  p2 = polyfit(log(hs),log(err_rms),1);
  fprintf('best-fit rms rate = %g\n', p2(1));
end
