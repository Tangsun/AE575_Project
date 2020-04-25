function dxpdt = swlin_xp(t, xp, nx)
    global A Q
    x = xp(1: nx);
    p = xp(nx+1: end);
    cost_i = zeros(1, 2);
    for i = 1: 2
        cost_i(i) = p'*A{i}*x;
    end
    [~, mode_i] = min(cost_i);
    dxdt = A{mode_i}*x;
    dpdt = -(Q*x + A{mode_i}'*p);
    dxpdt = [dxdt; dpdt];
end