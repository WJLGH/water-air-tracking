padding = .2;
t = 0.1:padding:1;
a = 3*t;
% 不指定间隔大小默认为1
gradient(a)
% ans =  0.6000    0.6000    0.6000    0.6000    0.6000
% 指定两个点之间的间隔为padding 相当于gradient(a)/padding
gradient(a,padding)
% ans = 3.0000    3.0000    3.0000    3.0000    3.0000

x = -2:0.2:2;
y = x';
z = x .* exp(-x.^2 - y.^2);
% 必须设置两个点之间的间隔，否则求出的梯度结果不对。
[px,py] = gradient(z,0.2);