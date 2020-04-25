function e = single_shoot(p_0)
    global x0 nx
    [~, xp] = ode45(@(t, xp) swlin_xp(t, xp, nx), [0, 2], [x0; p_0]);
    p_f = xp(end, nx+1: end)';
    e = p_f - 0;
end