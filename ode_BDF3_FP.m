function [t,u] = ode_BDF3_FP(odefun, tspan, y0m, Nh, tol, nmax, varargin)
%Prepare the parameters for the method
bm1 = 6/11; a0 = 18/11; a1 = -9/11; a2 = 2/11;

%initializing the solution vectors and includes initial condition
t = linspace(tspan(1),tspan(2),Nh+1); h = t(2) - t(1);
u = zeros(Nh + 1, size(y0m, 2)); u(1:3,:) = y0m;

%iterates the fixed point method for the three point backward difference
%formula
for n = 3:Nh
    vn = a0*u(n,:)' + a1*u(n-1,:)' + a2*u(n-2,:)';
    fun_tmp = @(x) vn + h*bm1*odefun(t(n+1), x, varargin{:});
    x_old = u(n,:)';
    for i = 1:nmax
        x_new = fun_tmp(x_old);
        if norm(x_new - x_old) <= tol
            disp(['number of iterations used = ', int2str(i)]);
            break;
        end
        x_old = x_new;
    end
u(n+1, :) = x_new';
end
t = t';
return;
end