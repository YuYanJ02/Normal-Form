function [coord_forward, coord_inverse] = build_coordinate_transforms(G_all, max_order)
    % 构建坐标变换
    
    syms q1 q2 q3 p1 p2 p3 complex_


    % 如果只需要双精度浮点数
    expr_double = vpa(G_all{3},3);

    % 初始化为恒等变换
    coord_forward.q1 = q1;
    coord_forward.q2 = q2;
    coord_forward.q3 = q3;
    coord_forward.p1 = p1;
    coord_forward.p2 = p2;
    coord_forward.p3 = p3;
    
    coord_inverse.q1 = q1;
    coord_inverse.q2 = q2;
    coord_inverse.q3 = q3;
    coord_inverse.p1 = p1;
    coord_inverse.p2 = p2;
    coord_inverse.p3 = p3;
    
    % 应用所有生成函数的变换
    for n = 3:max_order
        if ~isempty(G_all{n})
            % 正向变换
            coord_forward = apply_lie_transform_to_coords(coord_forward, G_all{n}, 1);
            
            % 逆向变换 (使用-Gn)
            coord_inverse = apply_lie_transform_to_coords(coord_inverse, -G_all{n}, 1);
        end
    end
    
    % 化简表达式
    fields_forward = fieldnames(coord_forward);
    fields_inverse = fieldnames(coord_inverse);
    for i = 1:length(fields_inverse)
        coord_forward.(fields_forward{i}) = simplify(coord_forward.(fields_forward{i}));
        coord_inverse.(fields_inverse{i}) = simplify(coord_inverse.(fields_inverse{i}));
    end
end

function coords_out = apply_lie_transform_to_coords(coords, Gn, tau)
    % 将Lie变换应用到坐标上
    
    syms q1 q2 q3 p1 p2 p3 complex_
    
    % 提取坐标场
    q_fields = {'q1', 'q2', 'q3', 'p1', 'p2', 'p3'};
    coords_out = coords;
    for i = 1:length(q_fields)
        field = q_fields{i};
        coord_expr = coords.(field);
        
        % 计算Lie级数: exp(tau * L_Gn) coord
        L1 = poisson_bracket(coord_expr, Gn);
        L2 = poisson_bracket(L1, Gn);
        L3 = poisson_bracket(L2, Gn);
        
        % 应用变换
        coord_transformed = coord_expr + tau * L1 + (tau^2/2) * L2 + (tau^3/6) * L3;
        coords_out.(field) = simplify(coord_transformed);
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

