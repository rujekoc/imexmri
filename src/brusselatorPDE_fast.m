function dy = brusselatorPDE_fast(t, y)
% usage: dy = brusselatorPDE_fast(t, y)
%
% Fast portion of right hand side 
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

% initialize RHS terms
du = zeros(nx,1);
dv = zeros(nx,1);
dw = zeros(nx,1);

% enforce stationary boundary conditions
du(1) = 0;  du(nx) = 0;  dv(1) = 0;  dv(nx) = 0; dw(1) = 0; dw(nx) = 0;


% reaction components
du(2:nx-1) = du(2:nx-1) + a - (w(2:nx-1)+1).*u(2:nx-1) + u(2:nx-1).*u(2:nx-1).*v(2:nx-1);
dv(2:nx-1) = dv(2:nx-1) + w(2:nx-1).*u(2:nx-1) - u(2:nx-1).*u(2:nx-1).*v(2:nx-1);
dw(2:nx-1) = dw(2:nx-1) + (b - w(2:nx-1))/ep - w(2:nx-1).*u(2:nx-1);


% combine together into output
dy = [du; dv ; dw];

% end function
