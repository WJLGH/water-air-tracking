%生成一个圆形轨迹
clear;
pos = [ 1 -0.5];
dt = 0.003;
T = 0.5;
t = 0:dt:T;
trajectory = zeros(length(t),3);
dist = zeros(length(t),1);
refv = 3; % m/s
refv = 10; % m/s
for i = 1:length(t)
    yaw = t(i)/T*360;
    dist(i) = norm(pos);
    movement = [sind(yaw),cosd(yaw)]*refv*dt;
    pos = pos + movement;
    trajectory(i,:) = [pos yaw];
end 

% figure
% plot(dist);