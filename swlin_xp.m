function dxpdt = swlin_xp(t, xp, nx)
    global A B Q R S
    x = xp(1: nx);
    p = xp(nx+1: end);
    cost_i = zeros(1, 3);
    for i = 1: 3
        cost_i(i) = p'*((A{i} - B{i}*inv(R)*S)*x - 0.5*B{i}*inv(R)*B{i}'*p);
    end
    [~, mode_i] = min(cost_i);
    u = -inv(R)*(S*x + B{mode_i}'*p);
    dxdt = A{mode_i}*x + B{mode_i}*u;
    dpdt = -(Q*x + S'*u + A{mode_i}'*p);
    dxpdt = [dxdt; dpdt];
end