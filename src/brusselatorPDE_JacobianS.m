function J = brusselatorPDE_JacobianS(t, y)
% usage: J = brusselatorPDE_JacobianS(t, y)
%
% Jacobian of the slow portion of the right hand side 


 J = brusselatorPDE_JacobianI(t,y) + brusselatorPDE_JacobianE(t,y);

end
