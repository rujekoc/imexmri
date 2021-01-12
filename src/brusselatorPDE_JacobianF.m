function J = brusselatorPDE_JacobianF(t, y)
% usage: J = brusselatorPDE_JacobianF(t, y)
%
% Jacobian of the fast portion of the right hand side,
% stiff brusselator PDE test problem.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% Modified October 2019 (Rujeko Chinomona)
% All Rights Reserved

% extract problem data
global Pdata;
nx  = Pdata.nx;
ep = Pdata.ep;

% extract solution components
u = y(1:nx);
v = y(nx+1:2*nx);
w = y(2*nx+1:3*nx);

% initialize Jacobian blocks terms
Juu = sparse([],[],[],nx,nx,3*nx);
Juv = sparse([],[],[],nx,nx,nx);
Juw = sparse([],[],[],nx,nx,nx);
Jvu = sparse([],[],[],nx,nx,nx);
Jvv = sparse([],[],[],nx,nx,3*nx);
Jvw = sparse([],[],[],nx,nx,nx);
Jwu = sparse([],[],[],nx,nx,nx);
Jwv = sparse([],[],[],nx,nx,nx);
Jww = sparse([],[],[],nx,nx,3*nx);

% reaction components
for j=2:nx-1
   Juu(j,j) = Juu(j,j) - (w(j)+1) + 2*u(j)*v(j);
end
for j=2:nx-1
   Juv(j,j) = Juv(j,j) + u(j)^2;
end
for j=2:nx-1
   Juw(j,j) = Juw(j,j) - u(j);
end

for j=2:nx-1
   Jvv(j,j) = Jvv(j,j) - u(j)^2;
end
for j=2:nx-1
   Jvu(j,j) = Jvu(j,j) + w(j) - 2*u(j)*v(j);
end
for j=2:nx-1
   Jvw(j,j) = Jvw(j,j) + u(j);
end

for j=2:nx-1
   Jww(j,j) = Jww(j,j) - 1/ep - u(j);
end
for j=2:nx-1
   Jwu(j,j) = Jwu(j,j) - w(j);
end
for j=2:nx-1
   Jwv(j,j) = Jwv(j,j) - 0;
end

% combine together into output
J = [Juu, Juv, Juw; Jvu, Jvv, Jvw; Jwu, Jwv, Jww];

% end function
