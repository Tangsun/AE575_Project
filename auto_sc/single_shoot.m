function e = single_shoot(p_0)
    global x0 T Ad Qd
    x = x0; p = p_0;
    for i = 1: length(T)-1
        cost_i = zeros(2, 1);
        for j = 1: 2
            cost_i(j) = p'*Ad{j}*x;
        end
        [~, mode] = min(cost_i);
        p = -Qd*x + Ad{mode}'*p;
        x = Ad{mode}*x;
    end
    e = p;
end