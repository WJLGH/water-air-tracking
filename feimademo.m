xmax = 5;
a = 1:0.2:xmax;
n = length(a);
t = zeros(n,2);

t(:,1) = a;
t(:,2) = 2;
a = t;
n2 = 1.34;
rpos = zeros(n,2);
syms x
for i=1:n
    r= norm(a(i,1));
    xp = vpasolve(n2*x/sqrt(1+x^2)-(r-x)/sqrt((r-x)^2+1));
    rpos(i,:) = [xp 1];
end

scatter(rpos(:,1),rpos(:,2));
hold on;
scatter(a(:,1),a(:,2));
hold off;
xlim([0,20])
ylim([0,2.2]);
% 结论，当入射角接近临界角的时候，可以折射到任意远的位置
% 在用solve解高次方程时，得到的解是一个带z的方程，
%   使用double函数将解变为实数解，
%   使用vpa获得相应的数值解
% s = solve(n*x/sqrt(1+x^2)-(r-x)/sqrt((r-x)^2+1))
% vpa(s)
% 上面的两步可以通过一下的 vpasolve 函数进行代替
% vpasolve(n2*x/sqrt(1+x^2)-(r-x)/sqrt((r-x)^2+1))