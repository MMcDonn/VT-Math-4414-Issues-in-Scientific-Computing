%%
%The following scripts were used in this code:
%int_eval_horner.m and int_ntcoef.m copyrighted by Tao Lin at Virginia Tech

%86
%The following calculates backward difference for the problem
[t,u] = ode_BDF3_FP(f,tspan,y0m,Nh,tol,100);

%The following solves the same problem using adams-moulton 4
[tam4,uam4] = ode_AM4_FixedPoint(f,tspan,y0m,Nh,tol,100);

%The following calculates the points against the debugging key in the HW
%Since t(1) = 0, t(2) = h, t(3) = 2h...
%Given deubbing solutions:

sol = [1; 1.001841470984808; 1.003683435618332; 1.005526024735161];
for i = 1:4
    BDF3_sol(i,1) = u(i);
    AM4_sol(i,1) = uam4(i);
    err_against_key(i,1) = norm(BDF3_sol(i,1) - sol(i));
    err_agaisnt_key(i,2) = norm(AM4_sol(i,1) - sol(i));    
end

%Find the point where t is approximately 0.75
find(abs(t-0.75) < eps)
% t(751) = 0.75
% u(751) = 2.316178684735357
t = tam4;
w = uam4;


%%
%88

%find the t value where t is approximately 5 seconds
find(abs(t_88-5) < eps)
plot(u_88(:, 1), u_88(:, 2), 2, 2, 'r*', 'Linewidth', 1.5)

%interpolate the value for where t is 1.5*pi seconds
tt = 1.5*pi;
%searches for points around desired t to be used for interpolatin nodes
for i = 1:length(t) - 1
if (tt >= t(i)) && (tt <= t(i+1))
I = i; break;
end
end
%creates nodes
t_int = t(I-1:I+2); x_int = w(I-1:I+2, 1); y_int = w(I-1:I+2, 2);
%creates coefficient vectors of polynomials for x and y
coef_x = int_ntcoef(t_int, x_int);
coef_y = int_ntcoef(t_int, y_int);
%evaluates the polynomials at desired t point
x = int_poly_eval_horner(tt, t_int, coef_x)
y = int_poly_eval_horner(tt, t_int, coef_y)
  

%intepolate the value where t is approximately 0.5*pi seconds
tt = 0.5*pi;
for i = 1:length(t) - 1
    if (tt >= t(i)) && (tt <= t(i+1))
    I = i; break;
    end
end
t_int = t(I-1:I+2); x_int = w(I-1:I+2, 1); y_int = w(I-1:I+2, 2);
  
t_int = t(I:I+1); x_int = w(I:I+1, 1); y_int = w(I:I+1, 2);
x_int_der = HW9_prob88_fun(t_int, x_int, c1,c2,b1,b2,d1,d2);
y_int_der = HW9_prob88_fun(t_int, y_int, c1,c2,b1,b2,d1,d2);
type = 0;
x = cubicspline(t_int, x_int, tt, type, x_int_der)
y = cubicspline(t_int, y_int, tt, type, y_int_der)

% x =
%    0.300779691569091
% y =
%    1.549536396145820
   
   
%%
%Calculates the arc length of the trajectory using a 
%composite trapezoidal quadrature

tspan = [0, 29.4602]; w0 = [1.15, 0, 0, 0.0086882909];
options=odeset('RelTol',1.e-6, 'AbsTol', 1.e-6);
[t, w]= ode45(@earth_moon ,tspan ,w0, options);
%Explicitly states the time derivatives of the solution 
x = w(:, 1); xdt = w(:, 2); y = w(:, 3); ydt = w(:, 4);
A = sqrt(xdt.^2 + ydt.^2);
traj_length = quadrature_trap_comp_f_vec_nonunif(A, t);

for i = 1:4
    plot(x(1),y(1),'o')
end
  
