% Reproduce case study A in the SQLR paper

clear all; close all; clc

%% Initialization
A = cell(3, 1); B = cell(3, 1);
A{1} = [-1, 1; -1, 1]; B{1} = [1, 1; 2, 1];
A{2} = [-1, 2; -1, 2]; B{2} = [2, 0; 1, 2];
A{3} = [-1, 3; -1, 3]; B{3} = [1, 0; 2, 2];

Q = 2*eye(2); R = 2*eye(2); S = eye(2);

x0 = [2; 1];

%% Formulate simulation for discrete-time system
Ad = cell(3, 1); Bd = cell(3, 1);
dt = 0.005;
for i = 1: 3
    sysc = ss(A{i}, B{i}, [], []);
    sysd = c2d(sysc, dt);
    Ad{i} = sysd.A;
    Bd{i} = sysd.B;
end

%% Calculate solution to algebraic Riccati equation
P = cell(3, 1);
for i = 1: 3
    [P{i}, ~, ~] = idare(Ad{i}, Bd{i}, Q, R, S', eye(2));
end

T = [0: dt: 2];
x = zeros(2, length(T)); u = zeros(2, length(T)); mode_i = zeros(1, length(T));
x(:, 1) = x0;
for ii = 1: length(T)-1
    lambda = cell(3, 1); cost_i = zeros(1, 3);
    for i = 1: 3
        lambda{i} = P{i}*x(:, ii);
        cost_i(i) = lambda{i}'*((Ad{i} - Bd{i}*inv(R)*S)*x(:, ii) - 0.5*Bd{i}*inv(R)*Bd{i}'*lambda{i});
    end
    [~, mode_i(ii)] = min(cost_i);
    u(:, ii) = -inv(R)*(S*x(:, ii) + Bd{mode_i(ii)}'*lambda{mode_i(ii)});
    x(:, ii+1) = Ad{mode_i(ii)}*x(:, ii) + Bd{mode_i(ii)}*u(:, ii);
end

%% Plot results
figure(1); hold on
plot(T, mode_i, 'k-', 'LineWidth', 1);
figure(2); hold on
plot(T, x(1, :), 'b-', 'LineWidth', 1);
plot(T, x(2, :), 'g-', 'LineWidth', 1);
legend('x_1', 'x_2');
xlabel('Time [s]'); ylabel('State');
figure(3); hold on
plot(T, u(1, :), 'b-', 'LineWidth', 1);
plot(T, u(2, :), 'g-', 'LineWidth', 1);
legend('u_1', 'u_2');
xlabel('Time [s]'); ylabel('Control');