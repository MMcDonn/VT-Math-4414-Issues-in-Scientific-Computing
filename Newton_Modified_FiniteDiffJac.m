function [zero,res,niter,err,J] = Newton_Modified_FiniteDiffJac(fun,x0,x1,tol,nmax,varargin)
%This is a modified Newton's method which employs a finite difference
%Jacobian. fun is the function which is to be evaluated. x0 and x1 are the
%initial guesses for a zero. tol and nmax are the user specified error
%tolerance and iteration limit.
n = length(x1); niter = 0; J = zeros(n,n);
err = norm(x1-x0);
while niter < nmax & err >= tol;
    %h = x1-x0;
    h = 10^(-3);
    F = fun(x1);
    for j = 1:n
        e = zeros(n,1); e(j) = 1;
        J(:,j) = (fun(x1 + h(j)*e) - F)/h(j);
    end
    delta = -J\F;
    x = x1 + delta;
    err = norm(delta);
    zero = x;
    x0 = x1;
    x1 = x;
    niter = niter +1;
    res = norm(fun(zero));
end
end

