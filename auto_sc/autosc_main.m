% Reproduce case study A in the SQLR paper
clear all; close all; clc

%% Initialization
global A Ad x0 nx Qd T
A = cell(2, 1);
A{1} = [-10, 5; 3, -2];
A{2} = [3, 1; 7, -20];

Q = eye(2);

x0 = [1; 6];

%% Formulate simulation for discrete-time system
dt = 0.01; T = [0: dt: 2.5];
Ad = cell(2, 1);
for i = 1: 2
    sys_c = ss(A{i}, [], [], []);
    sys_d = c2d(sys_c, dt);
    Ad{i} = sys_d.A;
end
sysp_c = ss(A{1}, Q, [], []);
sysp_d = c2d(sysp_c, dt);
Qd = sysp_d.B;

nx = length(x0);
p0_guess = [1; 0.4];
p0 = fsolve(@single_shoot, p0_guess);

x = zeros(2, length(T)); p = zeros(2, length(T)); mode = zeros(1, length(T));
x(:, 1) = x0; p(:, 1) = p0;
for i = 1: length(T)-1
    cost_i = zeros(2, 1);
    for j = 1: 2
        cost_i(j) = p(:, i)'*Ad{j}*x(:, i);
    end
    [~, mode(i)] = min(cost_i);
    p(:, i+1) = -Qd*x(:, i) + Ad{mode(i)}'*p(:, i);
    x(:, i+1) = Ad{mode(i)}*x(:, i);
end

%% Plot results
figure(1); hold on
plot(T(1: end-1), mode(1: end-1), 'k-', 'LineWidth', 1);
xlabel('Time [sec]');
ylabel('Mode');
axis([0, 2.5, 0.8, 2.2]);
figure(2); hold on
subplot(4, 1, 1);
plot(T, x(1, :), 'b-', 'LineWidth', 1);
axis([0, 2.5, 0, 2]);
xlabel('Time [s]'); ylabel('x_1');
subplot(4, 1, 2);
plot(T, x(2, :), 'b-', 'LineWidth', 1);
axis([0, 2.5, 0, 10]);
xlabel('Time [s]'); ylabel('x_2');
subplot(4, 1, 3);
plot(T, p(1, :), 'b-', 'LineWidth', 1);
axis([0, 2.5, 0, 10]);
xlabel('Time [s]'); ylabel('p_1');
subplot(4, 1, 4);
plot(T, p(2, :), 'b-', 'LineWidth', 1);
axis([0, 2.5, 0, 10]);
xlabel('Time [s]'); ylabel('p_2');
figure(3); hold on
plot(x(1, :), x(2, :), 'k-', 'LineWidth', 1);
xlabel('x_1'); ylabel('x_2');
axis([-3, 5, 0, 6]);