function [H_poly, c_coeffs, T_cell] = expand_hamiltonian(mu, c_coeffs, max_order)
    % 根据 Jorba & Masdemont (1998) 展开共线平动点附近的 Hamiltonian 为多项式
    % 输入：
    %   mu                - 质量参数 (0 < mu <= 0.5)
    %   libration_point  - 1,2,3 对应 L1, L2, L3
    %   max_order        - 展开的最高阶数（包括二次项）
    % 输出：
    %   H_poly    - 符号多项式 Hamiltonian H(x,y,z,px,py,pz)
    %   c_coeffs  - 计算得到的 c_n 系数向量（索引从2开始）
    %   T_cell    - 齐次多项式 T_n(x,y,z) 的元胞数组

    % --- 1. 计算平动点距离参数 gamma 和系数 c_n ---
   
    % --- 2. 定义符号变量 ---
    syms x y z pxx pyy pzz real

    % --- 3. 二次部分（公式 (7)）---
    c2 = c_coeffs(2);
    H2 = 1/2*(pxx^2 + pyy^2 + pzz^2) + y*pxx - x*pyy ...
         - c2*x^2 + c2/2*(y^2 + z^2);

    % --- 4. 递归生成齐次多项式 T_n(x,y,z) 并累加高阶项 ---
    %     T_n = rho^n * P_n(x/rho) 是 x,y,z 的 n 次齐次多项式
    %     递归： T0 = 1, T1 = x
    %            Tn = (2*n-1)/n * x * T_{n-1} - (n-1)/n * (x^2+y^2+z^2) * T_{n-2}
    T_cell = cell(1, max_order);
    T_cell{1} = sym(1);        % T0 (实际不使用，但递归需要)
    T_cell{2} = x;             % T1
    for n = 2:max_order
        Tn = (2*n-1)/n * x * T_cell{n} - (n-1)/n * (x^2+y^2+z^2) * T_cell{n-1};
        T_cell{n+1} = expand(Tn);   % 确保展开为多项式
    end

    % 累加高阶势能项： - sum_{n=2}^{max_order} c_n * T_n
    % 注意：论文中 n=2 项已经在 H2 中单独处理，因此这里 n 从 3 开始
    H_high = sym(0);
    for n = 3:max_order
        if n <= length(c_coeffs) && c_coeffs(n) ~= 0
            H_high = H_high - c_coeffs(n) * T_cell{n+1}; % T_n 对应索引 n+1
        end
    end

    % --- 5. 合并并展开最终多项式 ---
    H_poly = expand(H2 + H_high);
end

% ----------------------------------------------------------------------
