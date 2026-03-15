function [T,xx] = main_Orbit_(x0)

mu = 1.21506683e-2; 
RelTol = 1e-13; 
AbsTol = 1e-13; 
MaxNum = 25; 
Maxtf = x0(8); % 飞行时间

auxdata.mu = mu;
auxdata.RelTol = RelTol;
auxdata.AbsTol = AbsTol;
auxdata.MaxNum = MaxNum;
auxdata.Maxtf = Maxtf;

options = odeset('RelTol',RelTol,'AbsTol',AbsTol);

[T , xx] = ode45(@fun_crtbp_ , [0 , Maxtf] , x0(1:6) , options , auxdata); %利用修正后的初值求解crtbp

% 画轨道
%figure('position' , [100 , 100 , 400 , 280]) % 绘制区域的位置和大小
%plot3(xx(: , 1) , xx(: , 2) , xx(: , 3) , 'k' , 'LineWidth' , 2); 
%xlabel('x/LU'); ylabel('y/LU'); zlabel('z/LU');
%axis equal;
%axis([-0.5 2 , -1 1 -0.1 0.1])


end
% -------------------------------------------------------------------------