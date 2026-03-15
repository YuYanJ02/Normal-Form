function [action_angle, H_actions] = transform_to_action_angle(K_normal, lambda, omega_y, omega_z)
    % 变换到作用-角变量 - 公式(4.8)
    
    syms q1 q2 q3 p1 p2 p3 complex_
    syms I1 I2 I3 phi1 phi2 phi3 real
    
    % 作用-角变量定义
    q1_expr = sqrt(I1) * exp(phi1);
    p1_expr = sqrt(I1) * exp(-phi1);
    
    q2_expr = sqrt(I2) * exp(1i*phi2);
    p2_expr = -1i * sqrt(I2) * exp(-1i*phi2);
    
    q3_expr = sqrt(I3) * exp(1i*phi3);
    p3_expr = -1i * sqrt(I3) * exp(-1i*phi3);
    
    % 将正规形哈密顿量中的变量替换
    H_actions = subs(K_normal, [q1, q2, q3, p1, p2, p3], ...
        [q1_expr, q2_expr, q3_expr, p1_expr, p2_expr, p3_expr]);
    
    H_actions = simplify(expand(H_actions));
    
    % 存储变换关系
    action_angle.q1 = q1_expr;
    action_angle.p1 = p1_expr;
    action_angle.q2 = q2_expr;
    action_angle.p2 = p2_expr;
    action_angle.q3 = q3_expr;
    action_angle.p3 = p3_expr;
    action_angle.I1 = I1;
    action_angle.I2 = I2;
    action_angle.I3 = I3;
    action_angle.phi1 = phi1;
    action_angle.phi2 = phi2;
    action_angle.phi3 = phi3;
end