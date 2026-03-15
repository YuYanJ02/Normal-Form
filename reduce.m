function new_expr_collected = reduce(expr,threshold)
% 定义符号变量
syms p1 p2 p3 q1 q2 q3

% 展开表达式，确保为单项式的和
expr_expanded = expand(expr);

% 提取系数和对应的单项式
[coeffs_list, terms] = coeffs(expr_expanded, [p1, p2, p3, q1, q2, q3]);

% 设置阈值

new_expr = 0;

% 循环处理每一项
for i = 1:length(coeffs_list)
    coeff = coeffs_list(i);
    coeff_val = vpa(coeff);          % 转换为高精度数值
    if abs(coeff_val) >= threshold
        new_expr = new_expr + coeff * terms(i);  % 保留该项
    end
end

% 简化新表达式，并按变量整理
new_expr = simplify(new_expr);
new_expr_collected = collect(new_expr, [p1, p2, p3, q1, q2, q3]);

