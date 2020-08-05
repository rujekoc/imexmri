function J = brusselatorPDE_JacobianE(t, y)
% usage: J = brusselatorPDE_JacobianE(t, y)
%
% Jacobian of the slow-explicit component of the right hand side,
% stiff brusselator PDE test problem.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% Modified October 2019 (Rujeko Chinomona)
% All Rights Reserved

% extract problem data
global Pdata;
a  = Pdata.a;
a1 = Pdata.a1;
a2 = Pdata.a2;
a3 = Pdata.a3;
b  = Pdata.b;
d1 = Pdata.d1;
d2 = Pdata.d2;
d3 = Pdata.d3;
nx  = Pdata.nx;
dx = Pdata.dx;
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

% advection components
for j=2:nx-1
   Juu(j,j-1) = Juu(j,j-1) + a1/dx/2.0;
   Juu(j,j+1) = Juu(j,j+1) + a1/dx/2.0;
end

for j=2:nx-1
   Jvv(j,j-1) = Jvv(j,j-1) + a2/dx/2.0;
   Jvv(j,j+1) = Jvv(j,j+1) + a2/dx/2.0;
end

for j=2:nx-1
   Jww(j,j-1) = Jww(j,j-1) + a3/dx/2.0;
   Jww(j,j+1) = Jww(j,j+1) + a3/dx/2.0;
end

% combine together into output
J = [Juu, Juv, Juw; Jvu, Jvv, Jvw; Jwu, Jwv, Jww];

% end function
