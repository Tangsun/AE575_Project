% Reproduce case study A in the SQLR paper

clear all; close all; clc
global Ad Bd B_pd Q R S x0 T

%% Initialization
A = cell(3, 1); B = cell(3, 1); eig_A = zeros(2, 3);
controllable = zeros(3, 1); stablizable = ones(3, 1);
Wa = cell(3, 1); Wa_det = zeros(3, 1);
A{1} = [-1, 1; -1, 1]; B{1} = [1, 1; 2, 1];
A{2} = [-1, 2; -1, 2]; B{2} = [2, 0; 1, 2];
A{3} = [-1, 3; -1, 3]; B{3} = [1, 0; 2, 2];
for i = 1: 3
    eig_A(:, i) = eig(A{i});
    controllable(i) = (rank([A{i}, A{i}*B{i}]) == 2);
%     [Wa{i}, ~] = CtrGram(A{i}, B{i});
%     Wa_det(i) = det(Wa{i});
    for j = 1: 2
        rj = rank([A{i} - eig(j, i)*eye(2), B{i}]);
        if rj ~= 2 && eig_A(j, i) >= 0
            stablizable(i) = 0;
        end
    end
end

Q = 2*eye(2); R = 2*eye(2); S = eye(2);

x0 = [2; 1];

%% Formulate simulation for discrete-time system
Ad = cell(3, 1); Bd = cell(3, 1);
dt = 0.02;
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

%% Calculate solution using TPBVP
T = [0: dt: 2];
B_pd = cell(3, 1);
for i = 1: 3
    sys_p = ss(A{i}, [Q, S'], [], []);
    sys_pd = c2d(sys_p, dt);
    B_pd{i} = sys_pd.B;
end

x_s = zeros(2, length(T)); u_s = zeros(2, length(T)); 
p_s = zeros(2, length(T)); mode_s_i = zeros(1, length(T));
x_s(:, 1) = x0; p_0_guess = [3.2; -1.4];
p_s_0 = fsolve(@single_shoot, p_0_guess);
p_s(:, 1) = p_s_0;

for i = 1: length(T)-1
    cost_i = zeros(3, 1);
    for j = 1: 3
        cost_i(j) = p_s(:, i)'*((Ad{j} - Bd{j}*inv(R)*S)*x_s(:, i) - 0.5*Bd{j}*inv(R)*Bd{j}'*p_s(:, i));
    end
    [~, mode] = min(cost_i);
    u_s(:, i) = -inv(R)*(S*x_s(:, i) + Bd{mode}'*p_s(:, i));
    p_s(:, i+1) = -B_pd{mode}*[x(:, i); u(:, i)] - Ad{mode}'*p_s(:, i);
    x_s(:, i+1) = Ad{mode}*x_s(:, i) + Bd{mode}*u_s(:, i);
end
    

%% Plot results
figure(1); hold on
plot(T, mode_i, 'k-', 'LineWidth', 1);
figure(2); hold on
plot(T, x(1, :), 'b-', 'LineWidth', 1);
plot(T, x(2, :), 'b--', 'LineWidth', 1);
plot(T, x_s(1, :), 'r-', 'LineWidth', 1);
plot(T, x_s(2, :), 'r--', 'LineWidth', 1);
legend('x_1 Riccati', 'x_2 Riccati', 'x_1 TPBVP', 'x_2 TPBVP');
xlabel('Time [s]'); ylabel('State');
figure(3); hold on
plot(T, u(1, :), 'b-', 'LineWidth', 1);
plot(T, u(2, :), 'g-', 'LineWidth', 1);
legend('u_1', 'u_2');
xlabel('Time [s]'); ylabel('Control');