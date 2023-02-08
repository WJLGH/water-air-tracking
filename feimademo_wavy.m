t= 0.01:0.01:0.3;

n = 1.34;
ha = 1;
hw = 1;
t_s = 0.1;
syms x
% 在用solve解高次方程时，得到的解是一个带z的方程，
%   使用double函数将解变为实数解，
%   使用vpa获得相应的数值解
% s = solve(n*x/sqrt(1+x^2)-(r-x)/sqrt((r-x)^2+1))
% vpa(s)
% 上面的两步可以通过一下的 vpasolve 函数进行代替
% vpasolve(n*x/sqrt(ha+x^2)-(r-x)/sqrt((r-x)^2+ha))
err = [];
ipos = [];
for r = t
    mret =vpasolve(sin(atan((r-x)/ha)-t_s)-n*sin(atan(x/hw)-t_s));
    appret = (r*hw+ha*hw*(n-1)*t_s) / (n*ha+hw);
    ipos = [ipos mret];
    err = [err abs(appret - mret)  ];
end
% plot(err); 误差最大小于10^-4，满足了近似的要求
h = ones(1,length(t));
scatter(ipos,h);
hold on;
scatter(t,h*2)
hold off;
ylim([0,2])