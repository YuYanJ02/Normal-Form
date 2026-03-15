function [H_diag, T_diag, xi_tilde] = hamiltonian_after_diagonalization(H_poly, mu, libration_point,max_order)
    % 将对角化后的哈密顿量进行变换
    syms x y z pxx pyy pzz real
    
    % 1. 对角化线性系统
    [C, lambda_x, omega_y, omega_z, T_diag, xi_tilde] = diagonalize_linear_system(mu, libration_point,max_order);
    
    % 2. 将哈密顿量变换到对角化坐标
    % H_diag(x_tilde, y_tilde, z_tilde, pxt, pyt, pzt) = H_poly(x, y, z, px, py, pz)
    
    % 获取变换关系
    x_expr = T_diag.inverse.x;
    y_expr = T_diag.inverse.y;
    z_expr = T_diag.inverse.z;
    px_expr = T_diag.inverse.px;
    py_expr = T_diag.inverse.py;
    pz_expr = T_diag.inverse.pz;
    
    % 将H_poly中的变量替换为对角化坐标表达式
    H_diag = subs(H_poly, [x, y, z, pxx, pyy, pzz], ...
        [x_expr, y_expr, z_expr, px_expr, py_expr, pz_expr]);
    
    H_diag = simplify(expand(H_diag));
    
    % 3. 验证线性部分已对角化
    verify_diagonalization(H_diag, lambda_x, omega_y, omega_z, xi_tilde);
    
    fprintf('哈密顿量对角化完成\n');
end

function verify_diagonalization(H_diag, lambda_x, omega_y, omega_z, xi_tilde)
    % 验证哈密顿量的二次部分已对角化
    
    % 提取二次项
    H2 = extract_homogeneous_part(H_diag, 2);
    
    % 期望的对角化形式: H2 = lambda_x * x_tilde * pxt + 
    %                   omega_y/2 * (y_tilde^2 + pyt^2) + 
    %                   omega_z/2 * (z_tilde^2 + pzt^2)
    
    % 检查交叉项是否为零
    syms x_tilde y_tilde z_tilde pxt pyt pzt real
    
    % 检查x_tilde * y_tilde项
    coeff_xy = diff(diff(H2, x_tilde), y_tilde);
    if coeff_xy ~= 0
        fprintf('警告: 发现非零交叉项 x_tilde*y_tilde: %.2e\n', double(coeff_xy));
    end
    
    % 检查x_tilde * pyt项
    coeff_xpy = diff(diff(H2, x_tilde), pyt);
    if coeff_xpy ~= 0
        fprintf('警告: 发现非零交叉项 x_tilde*pyt: %.2e\n', double(coeff_xpy));
    end
    
    % 提取期望项的系数
    coeff_xp = diff(diff(H2, x_tilde), pxt);
    coeff_yy = diff(diff(H2, y_tilde), y_tilde)/2;
    coeff_pp = diff(diff(H2, pyt), pyt)/2;
    coeff_zz = diff(diff(H2, z_tilde), z_tilde)/2;
    coeff_pzpz = diff(diff(H2, pzt), pzt)/2;
    
    fprintf('对角化验证:\n');
    fprintf('  λ_x项系数: %.6f (期望: %.6f)\n', double(coeff_xp), lambda_x);
    fprintf('  ω_y项系数(y^2): %.6f (期望: %.6f)\n', double(coeff_yy), omega_y/2);
    fprintf('  ω_y项系数(py^2): %.6f (期望: %.6f)\n', double(coeff_pp), omega_y/2);
    fprintf('  ω_z项系数(z^2): %.6f (期望: %.6f)\n', double(coeff_zz), omega_z/2);
    fprintf('  ω_z项系数(pz^2): %.6f (期望: %.6f)\n', double(coeff_pzpz), omega_z/2);
    
    tolerance = 1e-8;
    if abs(double(coeff_xp) - lambda_x) < tolerance && ...
       abs(double(coeff_yy) - omega_y/2) < tolerance && ...
       abs(double(coeff_pp) - omega_y/2) < tolerance && ...
       abs(double(coeff_zz) - omega_z/2) < tolerance && ...
       abs(double(coeff_pzpz) - omega_z/2) < tolerance
        fprintf('对角化验证通过!\n');
    else
        warning('对角化验证未通过');
    end
end






function [C, lambda_x, omega_y, omega_z, T_diag, xi_tilde] = diagonalize_linear_system(mu, libration_point,max_order)
    % 对角化CRTBP线性系统
    % 输入:
    %   mu: 质量参数
    %   libration_point: 平动点 (1:L1, 2:L2, 3:L3)
    % 输出:
    %   C: 对角化变换矩阵 (6x6)
    %   lambda_x: 鞍点方向特征值
    %   omega_y: 平面中心方向频率
    %   omega_z: 垂直中心方向频率
    %   T_diag: 对角化变换的符号表达式
    %   xi_tilde: 对角化坐标符号变量
    
    fprintf('开始对角化线性系统...\n');
    
    % 1. 计算平动点位置和线性化系数c2
    [x_L, lambda_x, omega_y, omega_z, c_coeffs] = compute_libration_point_params(mu, libration_point,max_order);
    c2 = c_coeffs(2);
    
    fprintf('线性化参数:\n');
    fprintf('  c2 = %.6f\n', c2);
    fprintf('  λ = %.6f\n', lambda_x);
    fprintf('  ωy = %.6f\n', omega_y);
    fprintf('  ωz = %.6f\n', omega_z);
    
    % 2. 计算缩放因子 (公式2.49后的描述)
    [s_lambda, s_omega_y] = compute_scaling_factors(lambda_x, omega_y, c2);
    
    % 3. 构造对角化变换矩阵C (公式2.54, 第19页矩阵)
    C = construct_diagonalization_matrix(lambda_x, omega_y, omega_z, c2, s_lambda, s_omega_y);
    J = [ 0 0 0 1 0 0;
          0 0 0 0 1 0;
          0 0 0 0 0 1;
         -1 0 0 0 0 0;
         0 -1 0 0 0 0;
         0 0 -1 0 0 0
        ];
    A=C'*J*C
    
    
    % 4. 创建符号变换
    [T_diag, xi_tilde] = create_symbolic_diagonalization(C, lambda_x, omega_y, omega_z);
    
    fprintf('对角化完成!\n');
end

function [s_lambda, s_omega_y] = compute_scaling_factors(lambda_x, omega_y, c2)
    % 计算缩放因子s_lambda和s_omega_y (公式第18页)
    
    % 公式: s_lambda = sqrt(2*lambda_x*((4+3*c2)*lambda_x^2 + 4+5*c2-6*c2^2))
    s_lambda = sqrt(2 * lambda_x * ((4 + 3*c2) * lambda_x^2 + 4 + 5*c2 - 6*c2^2));

    % 公式: s_omega_y = sqrt(omega_y*((4+3*c2)*omega_y^2 - 4 - 5*c2 + 6*c2^2))
    s_omega_y = sqrt(omega_y * ((4 + 3*c2) * omega_y^2 - 4 - 5*c2 + 6*c2^2));
    
    fprintf('缩放因子:\n');
    fprintf('  s_lambda = %.6f\n', s_lambda);
    fprintf('  s_omega_y = %.6f\n', s_omega_y);
end

function C = construct_diagonalization_matrix(lambda_x, omega_y, omega_z, c2, s_lambda, s_omega_y)
    % 构造对角化变换矩阵C (第19页矩阵)
    % 注意: 这里使用数值矩阵，保持高精度
    
    % 初始化6x6矩阵
    C = zeros(6, 6);
    
    % 第1列: 对应x方向 (鞍点)
    C(1, 1) = 2 * lambda_x / s_lambda;
    C(2, 1) = (lambda_x^2 - 2*c2 - 1) / s_lambda;
    C(3, 1) = 0;
    C(4, 1) = (lambda_x^2 + 2*c2 + 1) / s_lambda;
    C(5, 1) = (lambda_x^3 + (1 - 2*c2) * lambda_x) / s_lambda;
    C(6, 1) = 0;
    
    % 第2列: 对应y方向 (平面中心)
    C(1, 2) = 0;
    C(2, 2) = (-omega_y^2 - 2*c2 - 1) / s_omega_y;
    C(3, 2) = 0;
    C(4, 2) = (-omega_y^2 + 2*c2 + 1) / s_omega_y;
    C(5, 2) = 0;
    C(6, 2) = 0;
    
    % 第3列: 对应z方向 (垂直中心)
    C(1, 3) = 0;
    C(2, 3) = 0;
    C(3, 3) = 1 / sqrt(omega_z);
    C(4, 3) = 0;
    C(5, 3) = 0;
    C(6, 3) = 0;
    
    % 第4列: 对应px方向 (鞍点共轭)
    C(1, 4) = -2 * lambda_x / s_lambda;
    C(2, 4) = (lambda_x^2 - 2*c2 - 1) / s_lambda;
    C(3, 4) = 0;
    C(4, 4) = (lambda_x^2 + 2*c2 + 1) / s_lambda;
    C(5, 4) = -(lambda_x^3 + (1 - 2*c2) * lambda_x) / s_lambda;
    C(6, 4) = 0;
    
    % 第5列: 对应py方向 (平面中心共轭)
    C(1, 5) = 2 * omega_y / s_omega_y;
    C(2, 5) = 0;
    C(3, 5) = 0;
    C(4, 5) = 0;
    C(5, 5) = (-omega_y^3 + (1 - 2*c2) * omega_y) / s_omega_y;
    C(6, 5) = 0;
    
    % 第6列: 对应pz方向 (垂直中心共轭)
    C(1, 6) = 0;
    C(2, 6) = 0;
    C(3, 6) = 0;
    C(4, 6) = 0;
    C(5, 6) = 0;
    C(6, 6) = sqrt(omega_z);
    
    % 显示矩阵条件数
    cond_C = cond(C);
    fprintf('变换矩阵C的条件数: %.2e\n', cond_C);
    if cond_C > 1e10
        warning('变换矩阵C条件数很大，可能数值不稳定');
    end
end

function [T_diag, xi_tilde] = create_symbolic_diagonalization(C, lambda_x, omega_y, omega_z)
    % 创建符号对角化变换
    
    % 定义原始坐标符号变量 (平动点中心坐标)
    syms x y z pxx pyy pzz real
    
    % 定义对角化后坐标符号变量
    syms x_tilde y_tilde z_tilde pxt pyt pzt real
    
    % 创建坐标向量
    xi_original = [x; y; z; pxx; pyy; pzz];
    xi_tilde_vec = [x_tilde; y_tilde; z_tilde; pxt; pyt; pzt];
    
    % 使用数值矩阵C进行变换
    % 注意: 由于C是数值矩阵，我们需要将其转换为符号矩阵
    C_sym = sym(C);
    
    % 正向变换: xi_tilde = C^{-1} * xi_original
    % 逆向变换: xi_original = C * xi_tilde
    
    % 计算正向变换
    xi_tilde_expr = C_sym \ xi_original;
    
    % 计算逆向变换
    xi_original_expr = C_sym * xi_tilde_vec;
    
    % 提取各个分量
    T_diag.forward.x_tilde = xi_tilde_expr(1);
    T_diag.forward.y_tilde = xi_tilde_expr(2);
    T_diag.forward.z_tilde = xi_tilde_expr(3);
    T_diag.forward.pxt = xi_tilde_expr(4);
    T_diag.forward.pyt = xi_tilde_expr(5);
    T_diag.forward.pzt = xi_tilde_expr(6);
    
    T_diag.inverse.x = xi_original_expr(1);
    T_diag.inverse.y = xi_original_expr(2);
    T_diag.inverse.z = xi_original_expr(3);
    T_diag.inverse.px = xi_original_expr(4);
    T_diag.inverse.py = xi_original_expr(5);
    T_diag.inverse.pz = xi_original_expr(6);
    
    % 存储符号变量
    xi_tilde.x_tilde = x_tilde;
    xi_tilde.y_tilde = y_tilde;
    xi_tilde.z_tilde = z_tilde;
    xi_tilde.pxt = pxt;
    xi_tilde.pyt = pyt;
    xi_tilde.pzt = pzt;
    
    % 简化表达式
    fields = fieldnames(T_diag.forward);
    fields_inv = fieldnames(T_diag.inverse);
    for i = 1:length(fields)
        T_diag.forward.(fields{i}) = simplify(T_diag.forward.(fields{i}));
        T_diag.inverse.(fields_inv{i}) = simplify(T_diag.inverse.(fields_inv{i}));
    end
    
    fprintf('符号对角化变换创建完成\n');
end


function H_part = extract_homogeneous_part(H, degree, vars)
    % 从多项式H中提取指定次数的齐次部分
    % 输入:
    %   H: 符号多项式表达式
    %   degree: 要提取的次数
    %   vars: (可选) 变量列表，如果未提供则自动检测
    % 输出:
    %   H_part: 次数为degree的齐次部分
    
    % 如果没有提供变量列表，自动检测
    if nargin < 3
        vars = symvar(H);
        if isempty(vars)
            % 如果是常数，检查是否为零次
            if degree == 0
                H_part = H;
            else
                H_part = sym(0);
            end
            return;
        end
    end
    
    % 确保H是展开的多项式
    H_expanded = expand(H);
    expr_double = vpa(H_expanded,2)
    % 获取系数和单项式
    [coeffs_, terms] = coeffs(H_expanded, vars);
    
    % 初始化结果
    H_part = sym(0);
    
    % 遍历所有项
    for i = 1:length(terms)
        term = terms(i);
        
        % 计算该项的总次数
        total_degree = compute_monomial_degree(term, vars);
        
        % 如果次数匹配，则添加到结果中
        if total_degree == degree
            H_part = H_part + coeffs_(i) * term;
        end
    end
    
    % 简化结果
    H_part = simplify(H_part);
end


function deg = compute_monomial_degree(term, vars)
    % 计算单项式的总次数
    
    deg = 0;
    
    % 对每个变量计算指数并求和
    for j = 1:length(vars)
        var = vars(j);
        
        % 计算变量var在term中的指数
        % 方法：如果term包含var，则指数为diff(log(term), var)*var
        % 但更稳健的方法是使用符号计算
        
        % 将term表示为var的多项式
        [coeffs_var, exps] = coeffs(term, var);
        
        % 找到包含var的项
        for k = 1:length(exps)
            % 检查指数是否为var的幂
            if has(exps(k), var)
                % 计算var的指数
                exponent = diff(log(exps(k)), var) * var;
                deg = deg + double(exponent);
                break;
            end
        end
    end
end