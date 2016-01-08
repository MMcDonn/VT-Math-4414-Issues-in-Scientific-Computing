function J = FJacobian(x)
n = max(size(x));
J = zeros(n,n);
J(1,1) = 3; J(1,2) = x(3)*sin(x(2)*x(3)); J(1,3) = x(2)*sin(x(2)*x(3));
J(2,1) = 8*x(1); J(2,2) = 2 - 1250*x(2);
J(3,1) = -exp(-x(1)*x(2))*x(2); J(3,2) = -exp(-x(1)*x(2))*x(1); J(3,3) = 20;
end