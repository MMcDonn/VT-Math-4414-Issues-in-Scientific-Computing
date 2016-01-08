function [zero,res,niter,B] = Broyden_quasiNewton(fun,B,x0,tol,nmax, varargin)
n = length(x0);
niter = 0; x = x0;
err = tol+1;
F = fun(x, varargin{:});
while (err >= tol && niter < nmax)
    delta = -B\F;
    x = x + delta;
    F = fun(x, varargin{:});
    B = B + F*delta'/(delta'*delta);
    err = norm(delta);
    niter = niter + 1;
end
res = norm(F);
zero = x;
end

