%%
%78
gamma = 0.5; %rate of heat transfer
m = 1; %kilograms
C = 100; %Joules per kilogram per Kelvin
T_e = 200; %Kelvin
S = 24; %Surface area of 1x1x1 cube
eps = 5.670*10^(-8);
tspan = [0,200];
Nh1 = 200/0.5; Nh2 = 200/0.1;
T0 = 180;
f = @(t,T) -(eps*gamma*S*(T.^4 - (T_e).^4))./(m*C);

i = 1;
while u_Nh2(i) < T_specific
    i = i+1;
end
i
u_Nh2(i)

i = 1;
while u_24Surf(i) < T_specific
    i = i+1;
end
i
u_24Surf(i)


%%
%79
%various step sizes to try
h = [1/100; 1/110; 1/120; 1/130; 1/140; 1/150; 1/160; 1/170; ...
    1/180; 1/190; 1/200]
%function initialization
y_ex = @(t) (1/2)*(exp(t) - sin(t) - cos(t));
err = max(abs(u - y_ex(t)));
%calculation for number of steps
Nh = 1./h

%employs third order Runge-Kutta solver
j = 1;
y_ex = @(t) (1/2)*(exp(t) - sin(t) - cos(t));
for i = 1:length(Nh)
[t,u] = ode_rk33(@prob79fun,tspan,y0,Nh(i));
err(i) = max(abs(u - y_ex(t)));
end

%%
%81

%Finding all the time points in (2,5) incremented by 0.5 and assinging it
%to its own solution vector called 'app'
k = 1;
for j = 2:0.5:5
i = find(abs(t-j) < eps); app(k,:) = u(i,:);
k = k+1;
end