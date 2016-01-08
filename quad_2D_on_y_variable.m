function F = quad_2D_on_y_variable(x, integrand, fun_c, fun_d, ...
M, quad_simple, n_s)
F = zeros(size(x));
for i = 1:length(x);
    f_tmp = @(y) integrand(x(i), y);
    c = fun_c(x(i)); d = fun_d(x(i));
    F(i) = quadrature_comp(f_tmp, c, d, M, quad_simple, n_s);
end
return;