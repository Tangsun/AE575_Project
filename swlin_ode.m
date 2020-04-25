function dxdt = swlin_ode(t, x)
    global A B P R S
    R = 2*eye(2); S = eye(2);
    cost_i = zeros(1, 3);
    for i = 1: 3
        lambda = P{i}*x;
        cost_i(i) = lambda'*((A{i} - B{i}*inv(R)*S)*x - 0.5*B{i}*inv(R)*B{i}'*lambda);
    end
    [~, mode_i] = min(cost_i);
    lambda = P{mode_i}*x;
    u = -inv(R)*(S*x + B{mode_i}'*lambda);
    dxdt = A{mode_i}*x + B{mode_i}*u;
end