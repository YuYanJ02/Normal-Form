function [x_L, lambda, omega_y, omega_z, c_coeffs] = compute_libration_point_params(mu, point,max_order)
    % 计算平动点位置和线性化参数
    % 参考公式(2.41)-(2.43)
    
    if point == 1  % L1
        % 解方程(2.42)
        gamma = fzero(@(g) g^5 - (3-mu)*g^4 + (3-2*mu)*g^3 - mu*g^2 + 2*mu*g - mu, 0.1);
        x_L = 1 - mu - gamma;
    elseif point == 2  % L2
        % 解方程(2.43)
        gamma = fzero(@(g) g^5 + (3-mu)*g^4 + (3-2*mu)*g^3 - mu*g^2 - 2*mu*g - mu, 0.1);
        x_L = 1 - mu + gamma;
    elseif point == 3  % L3
        % 解方程(2.41)
        gamma = fzero(@(g) g^5 + (2+mu)*g^4 + (1+2*mu)*g^3 - (1-mu)*g^2 - 2*(1-mu)*g - (1-mu), 0.1);
        x_L = -mu - gamma;
    end
    
    % 计算c2系数 (公式(2.46)之后)
    if point == 1
        c2 = 1/gamma^3 * (mu + (1-mu)*(gamma/(1-gamma))^3);
    elseif point == 2
        c2 = 1/gamma^3 * (mu + (1-mu)*(gamma/(1+gamma))^3);
    elseif point == 3
        c2 = 1/gamma^3 * (1 - mu + mu*(gamma/(1+gamma))^3);
    end
    
    % 计算特征值 (公式(2.49))
    eta = (c2 - 2 + sqrt(9*c2^2 - 8*c2)) / 2;
    lambda = sqrt(eta);      % 鞍点方向特征值
    omega_y = sqrt(-(c2 - 2 - sqrt(9*c2^2 - 8*c2)) / 2); % 平面中心频率
    omega_z = sqrt(c2);       % 垂直中心频率
    
    % 计算高阶c_n系数
    % max_order = 3;
    c_coeffs = zeros(max_order, 1);
    c_coeffs(2) = c2;
    
    for n = 3:max_order
        if point == 1
            c_coeffs(n) = 1/gamma^3 * ((-1)^n*mu + (1-mu)*(gamma/(1-gamma))^(n+1));
        elseif point == 2
            c_coeffs(n) = 1/gamma^3 * (mu + (1-mu)*(gamma/(1+gamma))^(n+1));
        elseif point == 3
            c_coeffs(n) = (-1)^n/gamma^3 * (1-mu + mu*(gamma/(1+gamma))^(n+1));
        end
    end
end