%% �켣����
tic
clear;

pos = [0.8,-0.8];
dt = 0.006;
T = 0.6;
t = 0:dt:T;
% �켣��(waypoint)��λ�ú��ٶ�[ px py vx vy]
trajectory = zeros(length(t),4);
lv = 3; % m/s
yaw = -80;
v = [sind(yaw),cosd(yaw)]*lv;
movement = v*dt;
dist = zeros(length(t),1);

for i = 1:length(t)
%     yaw = -45; % -45�� 
    pos = pos + movement;
    dist(i) = norm(pos);
    trajectory(i,:) = [pos,v];
end 

% plot(trajectory(:,1),trajectory(:,2),"*")
% plot(dist);
% step = norm(trajectory(1,:)-trajectory(2,:))
  

% generate_circle_traj
% ǿ�����ݻ�ȡ
padding = 0.3;
% ��������
Tc = 1;
trace_ints = zeros(length(t)*Tc,2);
tmp = zeros(Tc,2);
tpos_ints = zeros(length(t),2);

for i = 1:length(t)   
    pos = trajectory(i,1:2);
    for j=1:Tc
        tmp(j,:) = get_intensity(pos,t(i),padding);   
        trace_ints(Tc*i+j,:) = tmp(j,:);
    end
    tpos_ints(i,:) = mean(tmp,1);
end

ha = 1;
hw = 1;
n = 1.34;
% ������ɢ��
theta = pi/2;
% ˮ��˥��ϵ��
c = 0.15;
fov = pi/2;
Ar = 0.01 ^ 2;
bw = -log(2) / log(cos(theta / 2));
tau = 0.9;
% �������Ĺ��ЧӦ��
R = 0.6;
P_tx = 0.03;

%״̬����
syms px py d vx vy;
sy = [ px py vx vy];
% ״̬����
f = [px+vx; py+vy;vx;vy];
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

% mints = zeros(size(tpos_ints));
% for i = 1:length(trajectory)
%     pos = trajectory(i,1:2);
%     mints(i,:)= double(subs(h,[px,py],pos))';
% end

tpos_ints(1,:) = mints(1,:);

% figure()
% hold on
% plot(tpos_ints(:,1),'*');
% plot(mints(:,1),"*");
% plot(tpos_ints(:,2),'--');
% plot(mints(:,2),'--');
% legend('p1','esp1','p2','es p2')
% % legend('p2','es p2')
% hold off
toc

%% EKF Ԥ��켣
% eye���ص�λ����
tic
n = length(f);
m = length(h);
P = eye(n);
q = 1e-10; % ���̱�׼��
Q = 1e-3*eye(n);   % ���̷���
% ����ֵ�ķ���
R = diag(std(tpos_ints)).^2;
R = eye(2)*10^(-9); %  1.6303
R = eye(2)*10^(-14); %    0.1726
R = eye(2)*10^(-15); %   0.3797
R = eye(2)*10^(-13); %     0.2629
R = eye(2)*10^(-14.5); %     0.1774
R = eye(2)*10^(-13.5); %      0.1733 *
R = eye(2)*10^(-14.3); %    0.1746
R = eye(2)*10^(-13.8); %    0.1661
%  e-18 0.3321
%  e-17.5 0.0536
%  e-16 0.18
% ��ʼ��״̬
xhat = [1;1;0;0]; 
% ���Ź���ֵ
xV = zeros(n,length(tpos_ints));
F = double(subs(F_ ,sy,xhat'));
xM = zeros(n,length(trajectory));
zV = tpos_ints';
% �����ƶ����λ�ù���
for k = 1:length(tpos_ints)   
    z = zV(:,k);  
%     % ����f���ſɱȾ�������x1��Ӧ�ƽ�ʽline2
%     % ���Ʊ���u����ʾϵͳ������
    x1 = double(subs(f,sy,xhat'));
%     F = double(subs(F_ ,sy,xhat'));
    % ���̷���Ԥ�⣬��Ӧline3
    P = F*P*F'+Q;
%     % ����h���ſɱȾ���
    z1 = double(subs(h,[px,py],x1(1:2)'));
    H = double(subs(H_,sy,x1'));    
    % ���������棬��Ӧline4
    % inv���������
    S = H*P*H'+R;
    K = P*H'*pinv(S);
%     % ״̬EKF����ֵ����Ӧline5
    xhat = x1+K*(z-z1);
%     % EKF�����Ӧline6
    P = P-K*H*P; 
    if mod(k,Tc) == 0
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
plot(xM(2,2:length(trajectory))',"*");
plot(trajectory(:,2),"*");
hold off
legend('x','tx','y','ty')
figure
plot(xM(3,2:length(trajectory))')
hold on 
plot(trajectory(:,3));
plot(xM(4,2:length(trajectory))',"*");
plot(trajectory(:,4),"*");
hold off
legend('vx','est_vx','vy','est_vy')
% xM(:,end)
mse = zeros(1,length(trajectory));
for i = 1:length(trajectory)
    mse(i) = norm(trajectory(i,1:2) - xM(1:2,i)');
end 
figure();
plot(mse);
mean(mse);
xlabel('time/slot')
ylabel('mean square error/m')
grid on
% xlim([0,25]);
mean(mse)
toc