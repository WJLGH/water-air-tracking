%% 数据准备阶段
generate_trajectory       
% size(trajectory);
padding = 0.2;
% 采样次数
Tc = 1;
trace_ints = zeros(length(t)*Tc,2);
tpos_ints = zeros(length(t),2);
period_ints = zeros(Tc,2);
for i = 1:length(t)
    for j=1:Tc
        pos = trajectory(i,1:2);
        trace_ints(Tc*i+j,:) = get_intensity(pos,t(i),padding);   
    end
    tpos_ints(i,:) = mean(trace_ints(Tc*i:Tc*(i+1),:));
end

ha = 1;
hw = 1;
n = 1.34;
% 波束发散角
theta = pi/2;
% 水中衰减系数
c = 0.15;
fov = pi/2;
Ar = 0.01 ^ 2;
bw = -log(2) / log(cos(theta / 2));
tau = 0.9;
% 接收器的光电效应比
R = 0.6;
P_tx = 0.03;

%状态方程
syms px py d vx vy;
sy = [ px py vx vy];
sy = [ px py];
% 状态方程
f = [px+vx; py+vy;vx;vy ];
f = [px; py; ];
%分别求取两个PD测量得到的光强
r = sqrt(px^2+py^2);
s = solve(n*d/sqrt(hw^2+d^2)-(r-d)/sqrt((r-d)^2+ha^2),d);
mid = vpa(s(1));
dw = sqrt(mid^2+hw^2);
da = sqrt((r-mid)^2+ha^2);
a_t = atan(mid/hw);
a_r = atan((r-mid)/ha);
Is = P_tx * (bw+1)/(2*pi*dw)*cos(a_t)^bw;
la = 1/(2*theta*da^2);
lw = exp(-c*dw);
A_eff = Ar*cos(a_r);
intensity = R*tau*la*lw*A_eff*Is;

r1 = sqrt(px^2+(py+padding)^2);
s1 = solve(n*d/sqrt(hw^2+d^2)-(r1-d)/sqrt((r1-d)^2+ha^2),d);
mid1 = vpa(s1(1));
dw1 = sqrt(mid1^2+hw^2);
da1 = sqrt((r1-mid1)^2+ha^2);
a_t1 = atan(mid1/hw);
a_r1 = atan((r1-mid1)/ha);
Is1 = P_tx * (bw+1)/(2*pi*dw1)*cos(a_t1)^bw;
la1 = 1/(2*theta*da1^2);
lw1 = exp(-c*dw1);
A_eff1 = Ar*cos(a_r1);
intensity1 = R*tau*la1*lw1*A_eff1*Is1;
h = [
    %RSS
    intensity;    intensity1;
    ];
F_ = jacobian(f,sy);
H_ = jacobian(h,sy);
% 根据模型预测的两个PD的光强
mints = [];
for i = 1:length(trajectory)
    pos = trajectory(i,1:2);
    mints = [
        mints;
        double(subs(h,[px,py],pos))';
        ];  
end

tpos_ints(1,:) = mints(1,:);
figure()
hold on
plot(tpos_ints(:,1),'*');
plot(mints(:,1),"*");
plot(tpos_ints(:,2),'--');
plot(mints(:,2),'--');
legend('p1','esp1','p2','es p2')
% legend('p2','es p2')
hold off

%%
% eye返回单位矩阵
n = length(f);
m = length(h);
P = eye(n);
q = 1e-10; % 过程标准差
Q = 1e-3*eye(n);   % 过程方差
% 测量值的方差
R = diag([1e-17 1e-17]);
% 初始化估计状态 X_{0|0}
xhat = [pos';0;0]; 
xhat = [pos']; 
% 最优估计值
xV = zeros(n,length(tpos_ints));

F = double(subs(F_ ,sy,xhat'));
xM = zeros(n,length(trajectory));
zV = trace_ints';
% 进行移动点的位置估计
for k = 1:length(trace_ints)   
    z = zV(:,k);
%     % 计算f的雅可比矩阵，其中x1对应黄金公式line2
%     % 控制变量u，表示系统的输入
    x1 = double(subs(f,sy,xhat'));
%     F = double(subs(F_ ,sy,xhat'));
    % 过程方差预测，对应line3
    P = F*P*F'+Q;
%     % 计算h的雅可比矩阵
    z1 = double(subs(h,sy,x1'));
    H = double(subs(H_,sy,x1'));    
    % 卡尔曼增益，对应line4
    % inv返回逆矩阵
    S = H*P*H'+R;
    K = P*H'*pinv(S);
%     % 状态EKF估计值，对应line5
    xhat = x1+K*(z-z1);
%     % EKF方差，对应line6
    P = P-K*H*P; 
    if mod(k,Tc) ==0
        xM(:,k/Tc) = xhat;
    end 
end
figure
plot(xM(1,2:length(trajectory))',xM(2,2:length(trajectory))')
hold on 
plot(trajectory(:,1),trajectory(:,2));
hold off
xlabel('x/m');
ylabel('y/m');
legend('EKF estimated trajectory','AUV real trajectory')
grid on



figure
plot(xM(1,2:length(trajectory))')
hold on 
plot(trajectory(:,1));
plot(xM(2,2:length(trajectory))');
plot(trajectory(:,2));
hold off
legend('x','tx','y','ty')
% xM(:,end)
mse = zeros(1,length(trajectory));
for i = 1:length(trajectory)
    mse(i) = norm(trajectory(i,:) - xM(:,i)');
end 
figure();
plot(mse);
mean(mse);
xlabel('time/slot')
ylabel('mean square error/m')
grid on
xlim([0,25]);