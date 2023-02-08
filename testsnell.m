xmax = 5;
a = 1:0.2:xmax;
n = length(a);
t = zeros(n,2);
t(:,1) = a;
t(:,2) = 2;
a = t;


n1 = 1.34;
n2 = 1;
errors = zeros(n);
for i=1:n
    phi = rand()*2*pi;
%     phi = 0
%     phi = pi/4;
%     phi = pi*2/4;
%     phi = pi*3/4;
%     phi = pi*5/4;
%     phi = pi*6/4;
%     phi = pi*7/4;

    alpha = cos(phi);
    beta = sin(phi);
    rpos3 = [a(i,1)*alpha a(i,1)*beta 2];
    I = [rpos(i,1)*alpha rpos(i,1)*beta 1];
    N = [0  0 -1];
    T_r = rpos3 - I;
    N = N/norm(N);
    I = I/norm(I);
    dNT = -dot(N,I);
    C = n1/n2;
    tmp = 1-C^2*(1 - dNT^2);
    T = (C*dNT-sqrt(tmp))*N + C*I;
    T = T/norm(T);
    T_r = T_r/norm(T_r);
    errors(i) = abs(acos(dot(T,T_r)));
end
% 与真实的折射定律相差很小验证了向量折射定律的正确性
plot(errors)