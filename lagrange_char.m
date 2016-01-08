function phi = lagrange_char(z,x_nodes,k)
x_no_k = [x_nodes(1:k-1),x_nodes(k+1:end)];
coef = poly(x_no_k);
phi = polyval(coef,z)/polyval(coef, x_nodes(k));
end