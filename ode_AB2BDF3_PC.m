function [t,u] = ode_AB2BDF3_PC(odefun, tspan, y0m, Nh, varargin)
% prepare the parameters for the methods:
%Adams-Bashforth 2:
ab_a0 = 1; ab_b0 = 3/2; ab_b1 = -1/2;
%Backward Difference Formula 3:
bd_bm1 = 6/11; bd_a0 = 18/11; bd_a1 = -9/11; bd_a2 = 2/11;

% prepare the nodes:
t = linspace(tspan(1),tspan(2),Nh+1); h = t(2) - t(1);
u = zeros(Nh + 1, size(y0m, 2)); u(1:3,:) = y0m;

for n = 3:Nh
    fun_n = odefun(t(n), u(n,:)', varargin{:});
    fun_nm1 = odefun(t(n-1), u(n-1,:)', varargin{:});
    BDF3_vn =  bd_a0*u(n,:)' + bd_a1*u(n-1,:)' + bd_a2*u(n-2,:)';
    AB2_vn = ab_a0*u(n,:)';
    fun_tmp = @(x) BDF3_vn + h*bd_bm1*odefun(t(n+1), x, varargin{:});
    
    %Predictor - AB2
    x_old = AB2_vn + h*(ab_b0*fun_n + ab_b1*fun_nm1);
    
    %Corrector - BDF3
    x_new = fun_tmp(x_old);
    u(n+1, :) = x_new';
end
t = t';
end