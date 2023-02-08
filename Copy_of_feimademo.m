r=  1.077;
n = 1.34;
syms x 
% 在用solve解高次方程时，得到的解是一个带z的方程，
%   使用double函数将解变为实数解，
%   使用vpa获得相应的数值解
% s = solve(n*x/sqrt(1+x^2)-(r-x)/sqrt((r-x)^2+1),x);
% s = subs(s,1.2);
% vpa(s);
% 上面的两步可以通过一下的 vpasolve 函数进行代替
vpasolve(n*x/sqrt(1+x^2)-(r-x)/sqrt((r-x)^2+1))