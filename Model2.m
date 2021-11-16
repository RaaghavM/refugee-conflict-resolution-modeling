clear; clc; close all;
k_NA = 0.0004;
k_NM = 0.0002;
k_AN = 0.0002;
k_AM = 0.0002;
k_MN = 0.0002;
k_MA = 0.0002;

%dictionary to define state names as numbers
N = 1;
A = 2;
M = 3;

initial_N = 122.157;
initial_A = 73.88;
initial_M = 103.963; %(total_population-initial_N)/2;
total_population = initial_N+initial_A+initial_M;

P = @(x) x(1)+x(2)+x(3);
move = @(c, x, from, to) c*(x(to)/P(x))*x(from)*x(to);
f = @(t,x) [move(k_AN, x, A, N) + move(k_MN, x, M, N) - move(k_NA, x, N, A) - move(k_NM, x, N, M),
            move(k_NA, x, N, A) + move(k_MA, x, M, A) - move(k_AN, x, A, N) - move(k_AM, x, A, M),
            move(k_NM, x, N, M) + move(k_AM, x, A, M) - move(k_MN, x, M, N) - move(k_MA, x, M, A)];

[t,xa] = ode45(f,[0 300],[initial_N initial_A initial_M]);

figure(1)
plot(t,xa(:,1), '-x','LineWidth', 1.5,'DisplayName','N(t)')
title('N(t), A(t) and M(t)')
xlabel('t'), ylabel('Number of People')
hold on

plot(t,xa(:,2),'-o','LineWidth', 1.5,'DisplayName','A(t)')
hold on

plot(t,xa(:,3),'-+','LineWidth', 1.5,'DisplayName','M(t)')
legend

figure(2)
plot3(xa(:,1),xa(:,2),xa(:,3))
xlabel('N(t)'), ylabel('A(t)'), zlabel('M(t)')
xlim([0 total_population]), ylim([0 total_population]), zlim([0 total_population])
grid on
title('Solution curve')

%move2 = @(c, x, from, to) c*(x(to)/P(x))*x(from)*x(to);
%syms X_N(t) X_A(t) X_M(t)
%eqns = [diff(X_N,t) == k_NA*, diff(X_A,t) == X_A, diff(X_M,t) == X_M];
%[x_N, x_A, x_M] = dsolve(eqns);

figure(3)
[x_N,x_A]=meshgrid(0:1:total_population,0:1:total_population);
x_M = ones(size(x_N));
x_M = x_M .* (total_population);
x_M = x_M - x_N - x_A;
dy = k_NA*(x_A./total_population).*x_A.*x_N + k_MA*(x_A./total_population).*x_A.*x_M - k_AN*(x_N./total_population).*x_N.*x_A - k_AM*(x_M./total_population).*x_A.*x_M;
dx = k_AN*(x_N./total_population).*x_A.*x_N + k_MN*(x_N./total_population).*x_N.*x_M - k_NA*(x_A./total_population).*x_N.*x_A - k_NM*(x_M./total_population).*x_N.*x_M;
dxu=dx./sqrt(dx.^2+dy.^2);
dyu=dy./sqrt(dx.^2+dy.^2);
quiver(x_N,x_A,dxu,dyu)
xlabel('N(t)'),ylabel('A(t)')

k_NA = 0.0004;
k_NM = 0.0002;
k_AN = 0.0002;
k_AM = 0.0002;
k_MN = 0.0002;
k_MA = 0.0002;
k_NA_all = 0:0.0001:0.005;
N_all = zeros(length(k_NA_all));
A_all = zeros(length(k_NA_all));
M_all = zeros(length(k_NA_all));
% for i = 1:length(k_NA_all)
%     syms N A M
%     k_AN = k_NA_all(i);
%     %eqns = [k_AN*N*A + k_MN*N*M == k_NA*A*A + k_NM*M*M, k_NA*N*A + k_MA*M*A == k_AN*N*N + k_AM*M*M, k_NM*M*N + k_AM*M*A == k_MN*N*N + k_MA*A*A];
%     eqns = [k_AN*N*A + k_MN*N*M == k_NA*A*A + k_NM*M*M, k_NA*N*A + k_MA*M*A == k_AN*N*N + k_AM*M*M, k_NM*M*N + k_AM*M*A == k_MN*N*N + k_MA*A*A];
%     %assume(k_AN > 0 & k_MN > 0 & k_NA > 0 & k_NM > 0 & k_MA > 0 & k_AM > 0);
%     assume(N >= 0);
%     assumeAlso(A >= 0);
%     assumeAlso(M >= 0);
%     assumeAlso(N + A + M == total_population);
%     S = solve(eqns, [N,A,M], 'MaxDegree', 3);
%     N_all(i) = S.N;
%     A_all(i) = S.A;
%     M_all(i) = S.M;
% end

figure(4);
plot(k_NA_all, N_all(:,1), 'b','LineWidth', 1.3,'DisplayName', 'N');
title('Unstable Equilibrium vs. k_A_N (P = 30)')
xlabel('k_A_N'), ylabel('Number of People')
hold on
plot(k_NA_all, A_all(:,1), 'r','LineWidth', 1.3,'DisplayName', 'A');
hold on
plot(k_NA_all, M_all(:,1), 'LineWidth', 1.3,'DisplayName', 'M');
legend('N','A','M');

% figure(5);
% plot(k_NA_all, A_all);
% title('A value of unstable equilibrium vs. k_N_A (P = 30)')
% xlabel('k_N_A'), ylabel('A')
