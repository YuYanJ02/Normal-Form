function [K_normal, G_all, transforms] = lie_transformation(H_complex, nu, max_order)
    % 执行Dragt-Finn Lie级数变换
    % nu = [lambda, 1i*omega_y, 1i*omega_z]
    
    syms q1 q2 q3 p1 p2 p3 complex_
    
    % 初始化
    K_current = H_complex;
    K_normal = 0;
    G_all = cell(max_order, 1);
    transforms.forward = cell(max_order, 1);
    transforms.inverse = cell(max_order, 1);
    
    % 提取二次部分
    K0 = extract_homogeneous_part(K_current, 2);
    K_normal = K0;
    
    fprintf('开始Lie级数变换 (N = %d)...\n', max_order);
    
    for n = 3:max_order
        fprintf('  处理 %d 阶项...\n', n);
        
        % 提取当前阶数的项
        Hn = extract_homogeneous_part(K_current, n);
        
        if ~isempty(Hn)
            % 计算生成函数Gn
            Gn = compute_generator(Hn, nu, n);
            G_all{n} = Gn;
            
            % 应用变换: K_new = exp(L_Gn) K_current
            K_current = apply_lie_transform(K_current, Gn, max_order,n);
            
            % 更新正规形 (保留共振项)
            Kn_resonant = extract_resonant_terms(K_current, nu, n);
            K_normal = K_normal + Kn_resonant;
            
            % 记录变换
            transforms.forward{n} = @(coord) apply_forward_transform(coord, Gn);
            transforms.inverse{n} = @(coord) apply_inverse_transform(coord, Gn);
        end
    end
    
    fprintf('Lie级数变换完成\n');
end

function Gn = compute_generator(Hn, nu, order_n)
    % 计算生成函数Gn - 公式(4.6)
    
    syms q1 q2 q3 p1 p2 p3 complex_
    
    % 提取Hn中的单项式
    [coeffs_, terms] = coeffs(Hn, [q1, q2, q3, p1, p2, p3]);
    Gn = 0;
    
    for i = 1:length(terms)
        term = terms(i);
        coeff = coeffs_(i);
        
        % 获取指数向量
        q1_val = 1;
        q2_val = 1;
        q3_val = 1;
        p1_val = 1;
        p2_val = 1;
        p3_val = 1;

        kq = [diff(term, q1)/q1, diff(term, q2)/q2, diff(term, q3)/q3];
        kp = [diff(term, p1)/p1, diff(term, p2)/p2, diff(term, p3)/p3];
        kq = subs(kq, {q1, q2, q3,p1, p2, p3}, {q1_val, q2_val, q3_val,p1_val, p2_val, p3_val});
        kp = subs(kp, {q1, q2, q3,p1, p2, p3}, {q1_val, q2_val, q3_val,p1_val, p2_val, p3_val});

        % 转换为数值
        kq_num = double(kq);
        kp_num = double(kp);
        
        % 检查是否为共振项 (kp == kq)
        if ~all(kq_num == kp_num)
            % 计算分母 <kp - kq, nu>
            denominator = sum((kp_num - kq_num) .* nu);
            
            if abs(denominator) > 1e-10
                % 非共振项，可以消除
                gn_coeff = -coeff / denominator;
                Gn = Gn + gn_coeff * term;
            end
        end
    end
    
    Gn = simplify(Gn);
end

function K_new = apply_lie_transform(K, Gn, max_order, n)
    % 应用 Lie 变换: exp(L_Gn) K
    % K: 当前哈密顿量
    % Gn: 第 n 阶生成函数（齐次多项式）
    % max_order: 目标正规形最高阶
    % n: Gn 的齐次阶数

    syms q1 q2 q3 p1 p2 p3 complex_
    vars = [q1 q2 q3 p1 p2 p3];

    % 起始：L_G^0 K = K
    K_new   = K;
    L_power = K;

    % 估计需要的最高 Lie 幂次：
    %   每次作用 L_G (G 为 n 阶) 约把阶数提升 (n-2),
    %   从 2 阶开始，直到不超过 max_order:
    %       2 + k_max*(n-2) <= max_order
    %   =>  k_max = floor((max_order-2)/(n-2))
    if n <= 2
        k_max = 1;   % 线性/二次生成函数，用一次修正就够
    else
        k_max = floor((max_order - 2) / (n - 2));
        if k_max < 1
            k_max = 1;
        end
    end

    for k = 1:k_max
        % 计算 L_G^k K
        L_power = poisson_bracket(L_power, Gn);

        % 对应级数项 1/k! * L_G^k K
        term = (1/factorial(k)) * L_power;

        % 累加并按 max_order 截断高阶项
        K_new = truncate_polynomial(K_new + term, max_order);
    end
end

function PB = poisson_bracket(f, g)
    % 计算泊松括号 {f, g}
    
    syms q1 q2 q3 p1 p2 p3 complex_
    
    PB = 0;
    for j = 1:3
        qj = eval(['q' num2str(j)]);
        pj = eval(['p' num2str(j)]);
        PB = PB + diff(f, qj)*diff(g, pj) - diff(f, pj)*diff(g, qj);
    end
    
    PB = simplify(PB);
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



function H_truncated = truncate_polynomial(H, max_order)
    % 使用MATLAB内置函数进行截断
    
    % 获取变量
    vars = symvar(H);
    
    if isempty(vars)
        H_truncated = H;
        return;
    end
    
    % 展开多项式
    H_expanded = expand(H);
    
    % 获取所有项
    [coeffs_, terms] = coeffs(H_expanded, vars);
    
    H_truncated = sym(0);
    
    for i = 1:length(terms)
        % 使用polynomialDegree计算次数
        term_degree = polynomialDegree(terms(i), vars);
        
        if term_degree <= max_order
            H_truncated = H_truncated + coeffs_(i) * terms(i);
        end
    end
end



function H_resonant = extract_resonant_terms(H, nu, order_n)
    % 使用符号方法提取共振项（更准确）
    
    % 定义符号变量
    syms q1 q2 q3 p1 p2 p3
    
    vars = [q1, q2, q3, p1, p2, p3];
    
    % 提取指定阶数的项
    H_n = extract_homogeneous_part(H, order_n, vars);
    
    if H_n == 0
        H_resonant = sym(0);
        return;
    end
    
    % 使用多项式分解
    [coeffs_, monomials] = coeffs(expand(H_n), vars);
    
    H_resonant = sym(0);
    
    for i = 1:length(monomials)
        monomial = monomials(i);
        coeff = coeffs_(i);
        
        % 计算指数向量
        exponents = compute_exponents_symbolic(monomial, vars);
        
        % 分离q和p的指数
        kq = exponents(1:3);
        kp = exponents(4:6);
        
        % 计算内积
        inner_product = dot(kp - kq, nu);
        
        % 检查是否为零
        % if simplify(inner_product) == 0
        if inner_product == 0
            H_resonant = H_resonant + coeff * monomial;
        end
    end
end

function exponents = compute_exponents_symbolic(monomial, vars)
    % 使用符号方法计算指数
    
    exponents = zeros(1, length(vars));
    
    for i = 1:length(vars)
        var = vars(i);
        
        % 计算变量在单项式中的指数
        % 方法：log(monomial)对var求导
        try
            % 先检查是否包含该变量
            if has(monomial, var)
                % 使用多项式次数函数
                exponents(i) = polynomialDegree(monomial, var);
            else
                exponents(i) = 0;
            end
        catch
            % 如果失败，使用字符串方法
            monomial_str = char(monomial);
            var_str = char(var);
            exponents(i) = extract_exponent(monomial_str, var_str);
        end
    end
end