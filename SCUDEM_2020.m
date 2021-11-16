clear; clc; close all;
k = 0.025; %relationship between N and A 
b = 0.05; %relationship between A and M
a = -0.01; %relationship between N and M

% model 1 f = @(t,x) [k*x(1)*x(2)+a*x(1)*x(3);-k*x(1)*x(2)+b*x(2)*x(3);-a*x(1)*x(3)-b*x(2)*x(3)];
% model 2 

[t,xa] = ode45(f,[0 1.5],[500 400 300]);

figure(1)
plot(t,xa(:,1),'DisplayName','N(t)')
title('N(t), A(t) and M(t)')
xlabel('t'), ylabel('Number of People')
hold on

plot(t,xa(:,2),'DisplayName','A(t)')
hold on

plot(t,xa(:,3),'DisplayName','M(t)')
legend

figure(2)
plot3(xa(:,1),xa(:,2),xa(:,3))
grid on
title('Solution curve')
