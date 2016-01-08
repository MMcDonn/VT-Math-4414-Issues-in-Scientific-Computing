function x = gaussian_quadr_nodes(a, b, n)
bx = gaussian_quadr_nodes_ref(n);
x = ((b-a)/2)*bx + (b+a)/2;
return;
