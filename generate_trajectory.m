
pos = [0.8,-0.54]+0.3;
dt = 0.01;
T = 0.6;
t = 0:dt:T;
% 轨迹点(waypoint)的位置和速度[ px py vx vy]
trajectory = zeros(length(t),3);
lv = 10; % m/s
yaw = -45;
v = [sind(yaw),cosd(yaw)]*lv;
movement = v*dt;
dist = zeros(length(t),1);

for i = 1:length(t)
%     yaw = -45; % -45° 
    pos = pos + movement;
    dist(i) = norm(pos);
%     trajectory(i,:) = [pos,movement];
    trajectory(i,:) = [pos,yaw];
end 

% subplot(1,2,1)
% plot(trajectory(:,1),trajectory(:,2),"*")
% subplot(1,2,2)
% plot(dist);
% step = norm(trajectory(1,:)-trajectory(2,:))
