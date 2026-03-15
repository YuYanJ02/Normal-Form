function [qp, diag_state, canon_state] = crtbp_state_to_complex( ...
    mu, libration_point, order_N, x, y, z, xd, yd, zd)
%CRTBP_STATE_TO_COMPLEX
%  将 CRTBP 共转系的物理状态 (x,y,z,xd,yd,zd)
%  变换到复数哈密顿坐标 (q1,q2,q3,p1,p2,p3)。
%
%  输入可以是标量，也可以是行/列向量（同长度）。
%
%  输出：
%   qp          : [N x 6]，每行 [q1 q2 q3 p1 p2 p3]
%   diag_state  : [N x 6]，对角化坐标 [x_tilde y_tilde z_tilde pxt pyt pzt]
%   canon_state : [N x 6]，正则坐标  [x y z px py pz]

    % 确保都是列向量
    x  = x(:);  y  = y(:);  z  = z(:);
    xd = xd(:); yd = yd(:); zd = zd(:);
    N  = numel(x);

    % ---------- 1. 物理 -> 正则哈密顿坐标 ----------
    px = xd - y;
    py = yd + x;
    pz = zd;

    canon_state = [x y z px py pz];

    % ---------- 2. 构造 H_poly & 对角化，得到 T_diag ----------
    [x_L, lambda, omega_y, omega_z, c_coeffs] = ...
        compute_libration_point_params(mu, libration_point, order_N);

    H_poly = expand_hamiltonian(mu, c_coeffs, order_N);

    [H_diag, T_diag, xi_tilde] = ...
        hamiltonian_after_diagonalization(H_poly, mu, libration_point, order_N);

    % 提取 forward 变换表达式
    syms xs ys zs pxx pyy pzz real
    vars_old = {xs, ys, zs, pxx, pyy, pzz};

    f_x_tilde = matlabFunction( ...
        subs(T_diag.forward.x_tilde, {'x','y','z','pxx','pyy','pzz'}, vars_old), ...
        'Vars', vars_old);
    f_y_tilde = matlabFunction( ...
        subs(T_diag.forward.y_tilde, {'x','y','z','pxx','pyy','pzz'}, vars_old), ...
        'Vars', vars_old);
    f_z_tilde = matlabFunction( ...
        subs(T_diag.forward.z_tilde, {'x','y','z','pxx','pyy','pzz'}, vars_old), ...
        'Vars', vars_old);
    f_pxt = matlabFunction( ...
        subs(T_diag.forward.pxt, {'x','y','z','pxx','pyy','pzz'}, vars_old), ...
        'Vars', vars_old);
    f_pyt = matlabFunction( ...
        subs(T_diag.forward.pyt, {'x','y','z','pxx','pyy','pzz'}, vars_old), ...
        'Vars', vars_old);
    f_pzt = matlabFunction( ...
        subs(T_diag.forward.pzt, {'x','y','z','pxx','pyy','pzz'}, vars_old), ...
        'Vars', vars_old);

    % 数值对角化坐标
    x_tilde = zeros(N,1);
    y_tilde = zeros(N,1);
    z_tilde = zeros(N,1);
    pxt     = zeros(N,1);
    pyt     = zeros(N,1);
    pzt     = zeros(N,1);

    for k = 1:N
        x_tilde(k) = f_x_tilde(x(k), y(k), z(k), px(k), py(k), pz(k));
        y_tilde(k) = f_y_tilde(x(k), y(k), z(k), px(k), py(k), pz(k));
        z_tilde(k) = f_z_tilde(x(k), y(k), z(k), px(k), py(k), pz(k));
        pxt(k)     = f_pxt    (x(k), y(k), z(k), px(k), py(k), pz(k));
        pyt(k)     = f_pyt    (x(k), y(k), z(k), px(k), py(k), pz(k));
        pzt(k)     = f_pzt    (x(k), y(k), z(k), px(k), py(k), pz(k));
    end

    diag_state = [x_tilde y_tilde z_tilde pxt pyt pzt];

    % ---------- 3. 对角化坐标 -> 复数坐标 ----------
    Qinv = 1/sqrt(2) * [1, -1i; -1i, 1];

    q1 = x_tilde;
    p1 = pxt;

    qp = zeros(N, 6);   % [q1 q2 q3 p1 p2 p3]

    for k = 1:N
        % y 方向
        yp_vec = Qinv * [y_tilde(k); pyt(k)];
        q2 = yp_vec(1);
        p2 = yp_vec(2);

        % z 方向
        zp_vec = Qinv * [z_tilde(k); pzt(k)];
        q3 = zp_vec(1);
        p3 = zp_vec(2);

        qp(k, :) = [q1(k) q2 q3 p1(k) p2 p3];
    end
end