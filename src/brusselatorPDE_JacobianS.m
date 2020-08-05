function J = brusselatorPDE_JacobianS(t, y)
% usage: J = brusselatorPDE_JacobianS(t, y)
%
% Jacobian of the slow portion of the right hand side,
% stiff brusselator PDE test problem.
%
% Rujeko Chinomona
% Department of Mathematics
% Southern Methodist University
% July 2020
% All Rights Reserved

 J = brusselatorPDE_JacobianI(t,y) + brusselatorPDE_JacobianE(t,y);

end
