function int_val = quadrature_NCC(f_name,a,b,n,varargin)

%NOTE:
% n is the number 'rectangles' to approximate entire [a,b]
% Trapezoidal formula:     n = 1 (two nodes -> linear approximation)
% Simpson's formula:    n = 2 (three nodes -> quadratic approximation)

nodes = newton_cotes_closed_nodes(a,b,n);
weights = newton_cotes_closed_weights(a,b,n);
int_val = sum(f_name(nodes,varargin{:}).*weights);

end