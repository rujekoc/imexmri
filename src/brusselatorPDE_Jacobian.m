function J = brusselatorPDE_Jacobian(t, y)
% usage: J = brusselatorPDE_Jacobian(t, y)
%
% Full Jacobian routine for the stiff brusselator PDE test problem.
%
% Daniel R. Reynolds
% Department of Mathematics
% Southern Methodist University
% Modified October 2019 (Rujeko Chinomona)
% All Rights Reserved

% extract problem data
global Pdata;
a1 = Pdata.a1;
a2 = Pdata.a2;
a3 = Pdata.a3;
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

% diffusion components
for j=2:nx-1
   Juu(j,j-1) = Juu(j,j-1) + d1/dx/dx;
   Juu(j,j)   = Juu(j,j) - 2*d1/dx/dx;
   Juu(j,j+1) = Juu(j,j+1) + d1/dx/dx;
end
for j=2:nx-1
   Jvv(j,j-1) = Jvv(j,j-1) + d2/dx/dx;
   Jvv(j,j)   = Jvv(j,j) - 2*d2/dx/dx;
   Jvv(j,j+1) = Jvv(j,j+1) + d2/dx/dx;
end
for j=2:nx-1
   Jww(j,j-1) = Jww(j,j-1) + d3/dx/dx;
   Jww(j,j)   = Jww(j,j) - 2*d3/dx/dx;
   Jww(j,j+1) = Jww(j,j+1) + d3/dx/dx;
end

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

% advection components
for j=2:nx-1
   Juu(j,j-1) = Juu(j,j-1) - a1/dx/2.0;
   Juu(j,j+1) = Juu(j,j+1) + a1/dx/2.0;
end

for j=2:nx-1
   Jvv(j,j-1) = Jvv(j,j-1) - a2/dx/2.0;
   Jvv(j,j+1) = Jvv(j,j+1) + a2/dx/2.0;
end

for j=2:nx-1
   Jww(j,j-1) = Jww(j,j-1) - a3/dx/2.0;
   Jww(j,j+1) = Jww(j,j+1) + a3/dx/2.0;
end

% combine together into output
J = [Juu, Juv, Juw; Jvu, Jvv, Jvw; Jwu, Jwv, Jww];

% end function
