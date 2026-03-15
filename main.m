%% CRTBP Normal Form Computation
clear; clc;

%% 生成对角复数化后的哈密尔顿函数表达式
% 1. 参数设置
mu = 0.0121505856;
libration_point = 1;
order_N = 10;  % 从低阶开始测试

% 2. 计算平动点参数
[x_L, lambda, omega_y, omega_z, c_coeffs] = compute_libration_point_params(mu, libration_point,order_N);

% 3. 展开哈密顿量（在平动点中心坐标）
[H_poly] = expand_hamiltonian(mu, c_coeffs, order_N);

% 4. 对角化变换
[H_diag, T_diag, xi_tilde] = hamiltonian_after_diagonalization(H_poly, mu, libration_point,order_N);

% 5. 复数化变换（在对角化坐标基础上）
nu = [lambda, 1i*omega_y, 1i*omega_z];
[H_complex, T_complex] = complexify_hamiltonian(H_diag, lambda, omega_y, omega_z);

% 6. 写到文本文件（给 C++ 用）
export_H_complex(H_complex, 'H_complex_terms.txt');


%% 生成对角复数化后的坐标转换表达式
% 1. CRTBP 轨道数据 
% 输入API网址,获取所需轨道初值
url = 'https://ssd-api.jpl.nasa.gov/periodic_orbits.api?sys=earth-moon&family=dro';
orbit = webread(url);
total = numel(orbit.data);
%XX = zeros(total,6);
disp(['dro总轨道数：', num2str(total)]);
k=8500;
i=1;
% 获取第k个dro初值
while k <= total
    numeric_cell = cellfun(@str2double,orbit.data{k,1} , 'UniformOutput', false);
    x0(i,:) = cell2mat(numeric_cell);
    k = k + 1;
    i = i + 1;
end
% % 2:1DRO
% x0 = [0.8089 ;0;0 ;0;0.5156 ; 0 ;0;3.1416]; 
% % Halo
% x0 = [1.14318259622366 ; 1.04e-27 ; 0.158428843658573 ; 5.78e-15 ; -0.222135177096744 ; -4.19e-15;0;3.13665420395033];
x0 = [1.02950089 ;0;-0.1868081 ;0;-0.11898 ; 0 ;0;1.609]; 
[T,XX] = main_Orbit_(x0');
x  = XX(:,1);
y  = XX(:,2);
z  = XX(:,3);
xd = XX(:,4);
yd = XX(:,5);
zd = XX(:,6);


% 2. 共转 -> 正则 -> 对角化 -> 复数   qp(k,:) = [q1 q2 q3 p1 p2 p3]
[qp, diag_state, canon_state] = crtbp_state_to_complex( ...
    mu, libration_point, order_N, x, y, z, xd, yd, zd);


% 3. 写到文本文件（给 C++ 用）
outFile = 'qp_points.txt';
fid = fopen(outFile, 'w');
for k = 1:size(qp,1)
    q1 = qp(k,1); q2 = qp(k,2); q3 = qp(k,3);
    p1 = qp(k,4); p2 = qp(k,5); p3 = qp(k,6);
    % 实部 虚部 形式，便于 C++ 读取
    fprintf(fid, '%.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e %.16e\n', ...
        real(q1), imag(q1), real(q2), imag(q2), real(q3), imag(q3), ...
        real(p1), imag(p1), real(p2), imag(p2), real(p3), imag(p3));
end
fclose(fid);





% % 6. Lie级数变换
% [K_normal, G_all, transforms] = lie_transformation(H_complex, nu, order_N);
% 
% K_normal = reduce(K_normal,1e-15);
% 
% % 7. 构建坐标变换
% [coord_forward, coord_inverse] = build_coordinate_transforms(G_all, order_N);
% % p1_forward = vpa(reduce(coord_forward.p1,1e-15),2)
% % p1_inverse = vpa(reduce(coord_inverse.p1,1e-15),2)
% 
% 
% 
% % 8. 变换到作用-角变量
% [action_angle, H_actions] = transform_to_action_angle(K_normal, lambda, omega_y, omega_z);
% 
% % 9. 完整的坐标变换链
% fprintf('\n完整的坐标变换链:\n');
% fprintf('1. 平动点中心坐标 -> 对角化坐标\n');
% fprintf('2. 对角化坐标 -> 复数坐标\n');
% fprintf('3. 复数坐标 -> 正规形坐标 (Lie变换)\n');
% fprintf('4. 正规形坐标 -> 作用-角变量\n');
% 
% expr_double = vpa(K_normal,9)
% expr_double = vpa(H_actions,9)
