% Reproduce case study A in the SQLR paper

clear all; close all; clc

%% Initialization
global A B P Q R S
A = cell(3, 1); B = cell(3, 1);
A{1} = [-1, 1; -1, 1]; B{1} = [1, 1; 2, 1];
A{2} = [-1, 2; -1, 2]; B{2} = [2, 0; 1, 2];
A{3} = [-1, 3; -1, 3]; B{3} = [1, 0; 2, 2];

Q = 2*eye(2); R = 2*eye(2); S = eye(2);

x0 = [2; 1];

%% Calculate solution to algebraic Riccati equation
P = cell(3, 1);
for i = 1: 3
    [P{i}, ~, ~] = icare(A{i}, B{i}, Q, R, S', eye(2), zeros(2));
end

%% Formulate simulation for continuous-time system
dt = 0.02; T = [0: dt: 2];
x = zeros(2, length(T)); u = zeros(2, length(T)); mode_i = zeros(1, length(T));
x(:, 1) = x0;
% for ii = 1: length(T)-1
%     cost_i = zeros(1, 3);
%     for i = 1: 3
%         lambda = P{i}*x(:, ii);
%         cost_i(i) = lambda'*((A{i} - B{i}*inv(R)*S)*x(:, ii) - 0.5*B{i}*inv(R)*B{i}'*lambda);
%     end
%     [~, mode_i(ii)] = min(cost_i);
%     [T_sol, X_sol] = ode45(@(t, x) swlin_ode(t, x), [T(ii), T(ii+1)], x(:, ii));
%     x(:, ii+1) = X_sol(end, :)';
% end
[T_sol, X_sol] = ode45(@(t, x) swlin_ode(t, x), [0, 2], x(:, 1));

%% Plot results
% figure(1); hold on
% plot(T, mode_i, 'k-', 'LineWidth', 1);
figure(2); hold on
plot(T_sol, X_sol(:, 1), 'b-', 'LineWidth', 1);
plot(T_sol, X_sol(:, 2), 'g-', 'LineWidth', 1);
legend('x_1', 'x_2');
xlabel('Time [s]'); ylabel('State');
% figure(3); hold on
% plot(T, u(1, :), 'b-', 'LineWidth', 1);
% plot(T, u(2, :), 'g-', 'LineWidth', 1);
% legend('u_1', 'u_2');
% xlabel('Time [s]'); ylabel('Control');