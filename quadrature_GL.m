function int_val = quadrature_GL(f_name, a, b, n, varargin)
x = gaussian_quadr_nodes(a, b, n);
w = gaussian_quadr_weights(a, b, n);
int_val = sum(w.*f_name(x, varargin{:}));
return;

for n = 1:10
    approx(n,1) = quadrature_NCC(f,a,b,n)
end