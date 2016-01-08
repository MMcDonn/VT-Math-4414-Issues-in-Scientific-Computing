function [t, u] = ode_CN_Broyden(odefun, tspan, y0, Nh, tol, nmax, varargin)
bm1 = 1/2; b0 = 1/2;
%employs the crank-nicolson method with broyden to solve a PDE that has
%been discretized on time and space
t = linspace(tspan(1),tspan(2),Nh+1); tau = t(2) - t(1);
u = zeros(Nh + 1, length(y0)); u(1,:) = y0;
B0 = eye(length(y0), length(y0));
for n = 1:Nh
vn = u(n,:)' + tau*b0*odefun(t(n), u(n,:)', varargin{:});
G = @(z) z - tau*bm1*odefun(t(n+1), z, varargin{:}) - vn;
z0 = u(n,:)' + tau*odefun(t(n), u(n,:)', varargin{:}); % better initial guess?
[z, res, niter, B0] = Broyden(G, B0, z0, tol, nmax);
u(n+1, :) = z';
end
t = t'; % to make t a column vector
return;
