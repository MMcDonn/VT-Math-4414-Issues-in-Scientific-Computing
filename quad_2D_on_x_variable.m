function F = quad_2D_on_x_variable(y, integrand, fun_a, fun_b, ...
M, quad_simple, n_s)
F = zeros(size(y));
for i = 1:length(y);
    f_tmp = @(x) integrand(x, y(i));
    a = fun_a(y(i)); b = fun_b(y(i));
    F(i) = quadrature_comp(f_tmp, a, b, M, quad_simple, n_s);
end
return;