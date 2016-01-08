function p = lagrange_interp(z, x, y, d_order)
%lagrange_char_der.m is used below and is a script created by Tao Lin
%at Virginia Tech
n = size(x,2);
x_nodes = x;
p = zeros(1,size(z,2));
phi = zeros(n,length(z));

if d_order == 0
    for k = 1: n
        phi(k,:) = lagrange_char(z, x_nodes, k);
    end
    for i = 1:n
        p = p + y(i)*phi(i,:);
    end
else
    for k = 1:n
        phi(k,:) = lagrange_char_der(z,x_nodes,k,d_order);
    end
    for i = 1:n
        p = p + y(i)*phi(i,:);
    end
end