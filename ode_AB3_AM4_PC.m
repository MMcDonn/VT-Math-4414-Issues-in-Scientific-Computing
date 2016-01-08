function [t, u] = ode_AB3_AM4_PC(odefun, tspan, y0m, Nh, varargin)
% prepare the parameters for the methods:
bb0 = 23/12; bb1 = -16/12; bb2 = 5/12;
bm1 = 9/24; b0 = 19/24; b1 = -5/24; b2 = 1/24;

% prepare the nodes:
t = linspace(tspan(1),tspan(2),Nh+1); h = t(2) - t(1);
u = zeros(Nh + 1, size(y0m, 2)); u(1:3,:) = y0m;
% tol = 10^(-12); nmax = 100;

for n = 3:Nh
    fun_n = odefun(t(n), u(n,:)', varargin{:});
    fun_nm1 = odefun(t(n-1), u(n-1,:)', varargin{:});
    fun_nm2 = odefun(t(n-2), u(n-2,:)', varargin{:});
    vn = u(n,:)' + h*(b0*fun_n + b1*fun_nm1 + b2*fun_nm2);
    fun_tmp = @(x) h*bm1*odefun(t(n+1), x, varargin{:}) + vn;
    
    % predictor
    x_old = u(n,:)' + h*(bb0*fun_n + bb1*fun_nm1 + bb2*fun_nm2);
    
    % corrector
    x_new = fun_tmp(x_old);
    u(n+1, :) = x_new';
end
t = t'; % to make t a column vector
return;