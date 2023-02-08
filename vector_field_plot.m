% %绘制向量场图
% %例⼀
% clear all;clc;
% [X,Y] = meshgrid(-2:.2:2,-3:.2:3);
% Z = X.*exp(-X.^2 - Y.^2);
% [DX,DY] = gradient(Z); %Dx为⽔平⽅向上的梯度，第⼀列元素为原矩阵第⼆列与第⼀列元素之差，
% %第⼆列元素为原矩阵第三列与第⼀列元素之差除以2，以此类推
% [DDX,DDY]=gradient(Z,0.1,0.25);%这个就是将默认除以的2，x轴换为0.1，y轴换为0.25.
% %[FX，FY]=gradient(F,HX,HY) HX，HY参数表⽰各⽅向相邻两点的距离
% contour(X,Y,Z,7);       %7个等级的等值线图
% hold on,colormap hsv
% quiver(X,Y,DX,DY,0.8);  %绘制梯度（向量）场
% %1.2代表所画出的向量长度的缩放系数
% hold off
% 
% %例⼆
% n=-2:.15:2;
% [X,Y,Z]=peaks(n);
% contour(X,Y,Z,10);colormap autumn
% [U,V]=gradient(Z);hold on
% quiver(X,Y,U,V);hold off
% legend('等值线','向量场');

%例三：三维向量场图绘制
%本例是绘制三维参数曲线的很好的例⼦
vz = 10;            % Velocity
a = -32;            % Acceleration
t = 0:.1:pi/2;
x = 2*cos(t);y = 3*sin(t);
z = vz*t + 1/2*a*t.^2;
plot3(x,y,z,'r');hold on%绘制曲线图像
u = gradient(x);v = gradient(y);w = gradient(z);
quiver3(x,y,z,u,v,w,0.2,'b');   %绘制切向量

% gradient 和diff区别
% gradient对矩阵来说是求数值梯度，
% diff对于矩阵来说是求差分
% 当然对于符号函数diff(func,x)用于求对于func对于x的偏导数。