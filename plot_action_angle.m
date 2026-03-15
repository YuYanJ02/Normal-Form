% ===== 1) 读文件 =====
fname = 'action_angle.txt';   % 改成你的路径/文件名
A = readmatrix(fname);

ReI1 = A(:,1);  ImI1 = A(:,2);
ReI2 = A(:,3);  ImI2 = A(:,4);
ReI3 = A(:,5);  ImI3 = A(:,6);
phi1 = A(:,7);  phi2 = A(:,8);  phi3 = A(:,9);

% 如果你有对应时间向量 t（和点数一致），用它替换 k
k = (1:size(A,1)).';

% ===== 2) 画作用量（通常看实部即可）=====
figure; 
subplot(3,1,1); plot(k, ReI1, 'LineWidth', 1.2); grid on; ylabel('I1 (Re)'); title('Actions');
subplot(3,1,2); plot(k, ReI2, 'LineWidth', 1.2); grid on; ylabel('I2 (Re)');
subplot(3,1,3); plot(k, ReI3, 'LineWidth', 1.2); grid on; ylabel('I3 (Re)'); xlabel('index');

% 如需检查虚部（理论上应接近 0）
figure;
subplot(3,1,1); plot(k, ImI1, 'LineWidth', 1.2); grid on; ylabel('Im(I1)'); title('Imag parts (should be small)');
subplot(3,1,2); plot(k, ImI2, 'LineWidth', 1.2); grid on; ylabel('Im(I2)');
subplot(3,1,3); plot(k, ImI3, 'LineWidth', 1.2); grid on; ylabel('Im(I3)'); xlabel('index');

% ===== 3) 画角变量（建议 unwrap）=====
phi2u = unwrap(phi2);
phi3u = unwrap(phi3);

figure;
subplot(3,1,1); plot(k, phi1, 'LineWidth', 1.2); grid on; ylabel('\phi_1'); title('Angles');
subplot(3,1,2); plot(k, phi2u, 'LineWidth', 1.2); grid on; ylabel('\phi_2 (unwrap)');
subplot(3,1,3); plot(k, phi3u, 'LineWidth', 1.2); grid on; ylabel('\phi_3 (unwrap)'); xlabel('index');

% ===== 4) 可选：看 I 的漂移量（稳定性指标）=====
I1d = ReI1 - ReI1(1);
I2d = ReI2 - ReI2(1);
I3d = ReI3 - ReI3(1);

figure;
plot(k, I1d, 'LineWidth', 1.2); hold on;
plot(k, I2d, 'LineWidth', 1.2);
plot(k, I3d, 'LineWidth', 1.2);
grid on; legend('\Delta I1','\Delta I2','\Delta I3');
title('Action drift');
xlabel('index'); ylabel('\Delta I');