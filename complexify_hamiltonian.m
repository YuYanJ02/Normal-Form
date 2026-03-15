function [H_complex, T] = complexify_hamiltonian(H_poly, lambda, omega_y, omega_z)
    % 复数化变换 - 公式(4.2)
    
    % 定义实坐标
    syms x_tilde y_tilde z_tilde pxt pyt pzt real
    
    % 定义复坐标
    syms q1 q2 q3 p1 p2 p3 complex_
    
    % 复数化变换矩阵
    Q = 1/sqrt(2) * [1, 1i; 1i, 1];
    
    % 变换关系
    x_tilde_expr = q1;
    pxt_expr = p1;
    
    % y方向变换
    y_vec = Q * [q2; p2];
    y_tilde_expr = y_vec(1);
    pyt_expr = y_vec(2);
    
    % z方向变换
    z_vec = Q * [q3; p3];
    z_tilde_expr = z_vec(1);
    pzt_expr = z_vec(2);
    
    % 将实哈密顿量中的变量替换为复变量表达式
    H_complex = subs(H_poly, [x_tilde, y_tilde, z_tilde, pxt, pyt, pzt], ...
        [x_tilde_expr, y_tilde_expr, z_tilde_expr, pxt_expr, pyt_expr, pzt_expr]);
    
    H_complex = expand(H_complex);
    
    % 存储变换关系
    T.q1 = x_tilde_expr;
    T.p1 = pxt_expr;
    T.q2 = q2;
    T.p2 = p2;
    T.q3 = q3;
    T.p3 = p3;
end