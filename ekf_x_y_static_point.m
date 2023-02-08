pos = [0.51,-0.5];
pos = [0.2,0.4];
padding = 0.2;
% pos = [1, 0.5];
t = 0:0.2:10;
% y = -1.5:0.1:1.5;
% y = -0.5:1:0.5;
ints = zeros(length(t),2);
for i = 1:length(t)
    ints(i,:) = get_intensity(pos,t(i),padding);
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
syms px py d;
sy = [ px py];
% ״̬����
f = [px; py ];
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

mints = double(subs(h,[px,py],pos))';
mints = repmat(mints,length(ints),1);
figure()
hold on
plot(ints(:,1),'*');
plot(mints(:,1),"*");
plot(ints(:,2),'--');
plot(mints(:,2),'--');
legend('p1','esp1','p2','es p2')
% legend('p2','es p2')
hold off

%%
%��ֵƫ�������Rֵ
%�����ٶ���������Rֵ
% eye���ص�λ����
n = length(f);
m = length(h);
P = eye(n);
q = 1e-10; % ���̱�׼��
Q = 1e-3*eye(n);   % ���̷���
% ����ֵ�ķ���
R = diag([1e-17 1e-17]);
R = diag(std(ints)).^2;
R = diag([1e-15 1e-15]);
% ״̬����ֵ
zV = ints';
% ��ʼ��״̬
xhat = [1.2;0.2]; 
% ���Ź���ֵ
xV = zeros(n,length(ints));
F = double(subs(F_ ,sy,xhat'));

for k = 1:length(ints)   
    z = zV(:,k);
%     % ����f���ſɱȾ�������x1��Ӧ�ƽ�ʽline2
%     % ���Ʊ���u����ʾϵͳ������
    x1 = double(subs(f,sy,xhat'));
%     F = double(subs(F_ ,sy,xhat'));
    % ���̷���Ԥ�⣬��Ӧline3
    P = F*P*F'+Q;
%     % ����h���ſɱȾ���
    z1 = double(subs(h,sy,x1'));
    H = double(subs(H_,sy,x1'));    
    % ���������棬��Ӧline4
    % inv���������
    S = H*P*H'+R;
    K = P*H'*pinv(S);
%     % ״̬EKF����ֵ����Ӧline5
    xhat = x1+K*(z-z1);
%     % EKF�����Ӧline6
    P = P-K*H*P; 
    xV(:,k) = xhat;
end
figure
plot(xV(1,:)')
hold on 
plot(xV(2,:)');
hold off
legend('x','y')
xV(:,end)
mse = zeros(1,length(ints));
for i = 1:length(ints)
    mse(i) = norm(pos' - xV(:,i));
end 
figure();
plot(mse);
mean(mse);