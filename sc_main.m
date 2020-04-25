% Reproduce case study A in the SQLR paper

clear all; close all; clc

%% Initialization
global A B P Q R S x0 nx
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
nx = length(x0);
p0_guess = [0; 0];
p0 = fsolve(@single_shoot, p0_guess);

[T_sol, xp_sol] = ode45(@(t, xp) swlin_xp(t, xp, nx), T, [x0; p0]);
x = xp_sol(:, 1: nx);
p = xp_sol(:, nx+1: end);

mode = zeros(length(T_sol), 1);
u = zeros(length(T_sol), 2);
for ii = 1: length(T_sol)
    cost_i = zeros(1, 3);
    for i = 1: 3
        cost_i(i) = p(ii, :)*((A{i} - B{i}*inv(R)*S)*x(ii, :)' - 0.5*B{i}*inv(R)*B{i}'*p(ii, :)');
    end
    [~, mode(ii)] = min(cost_i);
    u_ii = -inv(R)*(S*x(ii, :)' + B{mode(ii)}'*p(ii, :)');
    u(ii, :) = u_ii';
end


%% Plot results
figure(1); hold on
plot(T_sol, mode, 'k-', 'LineWidth', 1);
figure(2); hold on
plot(T_sol, x(:, 1), 'b-', 'LineWidth', 1);
plot(T_sol, x(:, 2), 'g-', 'LineWidth', 1);
legend('x_1', 'x_2');
xlabel('Time [s]'); ylabel('State');
figure(3); hold on
plot(T_sol, u(:, 1), 'b-', 'LineWidth', 1);
plot(T_sol, u(:, 2), 'g-', 'LineWidth', 1);
legend('u_1', 'u_2');
xlabel('Time [s]'); ylabel('Control');
figure(4); hold on
plot(T_sol, p(:, 1), 'b-', 'LineWidth', 1);
plot(T_sol, p(:, 2), 'g-', 'LineWidth', 1);
legend('p_1', 'p_2');
xlabel('Time [s]'); ylabel('Co-State');