function J = brusselatorPDE_JacobianI(t, y)
% usage: J = brusselatorPDE_JacobianI(t, y)
%
% Jacobian of the slow-implicit portion of the right hand side 
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


% combine together into output
J = [Juu, Juv, Juw; Jvu, Jvv, Jvw; Jwu, Jwv, Jww];

% end function
