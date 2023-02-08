% 首先将光线看做直线，水面看做平面找到，光线与平面的交点，
% 然后在这个交点附近计算最契合折射定律的点  (暴力枚举,求其差异最小点或者小于阈值的点)
function [ipos3,err] = find_incident_point(rpos,X,Y,ha,hw,h,hx,hy)
% rpos是一个二维向量，ha:平均水面高度，hw平均水深,h为海面的高度分布
%   
%  初始猜测的入射点
    ipos = abs(hw)/abs(ha+hw)*rpos;
    % 初始猜测的入射点在网格中的索引
    [xmax,ymax] = size(X);
    x = round( (rpos(1)- X(1,1))/(X(1,xmax)-X(1,1)) *xmax);
    y = round( (rpos(2)- Y(1,1))/(Y(ymax,1)-Y(1,1)) *ymax);
    span = 50;
    lx = max(1,x-span);
    rx = min(xmax,x+span);
    uy = max(1,y-span);
    dy = min(ymax,y+span);
    xspan = lx:1:rx;
    yspan = uy:1:dy;
    tx = ipos(1);
    ty = ipos(2);
    tz = hw;
    rpos3 = [ rpos ha+hw ];
    n1 = 1.34;
    n2 = 1;
    tarcos = -inf;
%     较小的区域导致算法无法找到最优的入射点，从而使得与实际估计差距很大，需要将整个波面的采样平面中搜索。
    for i = 1:xmax
%     for i = xspan
        for j = 1:ymax
%         for j = yspan
            ix = X(i,j);
            iy = Y(i,j);
            iz = h(i,j)+hw;
            idx = hx(i,j);
            idy = hy(i,j);
            N = [ idx idy -1];
            I = [ix iy iz];
            T_r = rpos3 - I;

            N = N/norm(N);
            I = I/norm(I);
            dNT = -dot(N,I);
            
            C = n1/n2;
            tmp = 1-C^2*(1 - dNT^2);
            if tmp < 0
                ctar = -inf;
            else 
                T = (C*dNT-sqrt(tmp))*N + C*I;
                ctar = dot(T_r,T)/norm(T_r)/norm(T);
            end
            if tarcos<ctar
                tx = ix;
                ty = iy;
                tz = iz;
                tarcos = ctar;
            end
        end
    end
    ipos3 =  [tx ty tz];
    err = abs(1-tarcos);
    if err == inf
        ipos3 = inf;
    end
%     T_r = rpos3 - ipos3;
%     da = norm(T_r);
%     dw = norm(ipos3);
%     theta_i = acos(tz/norm(ipos3));
%     theta_r = acos(T_r(3)/norm(T_r));
end
    

    
