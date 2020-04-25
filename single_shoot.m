function e = single_shoot(p_0)
    global Ad Bd B_pd Q R S x0 T
    
    p = p_0; x = x0;
    for i = 1: length(T)-1
        cost_i = zeros(3, 1);
        for j = 1: 3
            cost_i(j) = p'*((Ad{j} - Bd{j}*inv(R)*S)*x - 0.5*Bd{j}*inv(R)*Bd{j}'*p);
        end
        [~, mode] = min(cost_i);
        u = -inv(R)*(S*x + Bd{mode}'*p);
        p = -B_pd{mode}*[x; u] - Ad{mode}'*p;
        x = Ad{mode}*x + Bd{mode}*u;
    end
    
    e = p;
end