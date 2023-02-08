function  slope = difference_slop(A,d,padding)
% 使用差分近似计算梯度
% A为高度值，d为计算的方向，是X方向还y方向
% padding相邻两点的距离
    [n,m] = size(A);
    if d == 1
        slope = (A(:,2:m)-A(:,1:m-1)) / padding;
    else
        slope = (A(2:n,:)-A(1:n-1,:)) / padding;
    end 
end

