function int_val = quadrature_comp(fun_name, a, b, M, quad_simple_name, n_s, varargin)
x = linspace(a, b, M+1);
int_val = 0;

for i = 1:M
    int_val = int_val + ...
        quad_simple_name(fun_name, x(i), x(i+1), n_s, varargin{:});
end
end