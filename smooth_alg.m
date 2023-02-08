%Matlab
data_length = 201;
noise = wgn(1,data_length,2);  
for n=1:1:data_length  
    s(1,n) = 3*sin(0.2*n)+3*sin(0.4*n);  
end  
x = s + noise;
h=[-1,  3, -3,  1;
    3, -6,  3,  0;
   -3,  0,  3,  0;
    1,  4,  1,  0];
for n=1:1:data_length
    t = n*0.0001;
    if n<4
        y(1,n) = x(1,n);
    else
    temp = [x(1,n), x(1,n-1), x(1,n-2), x(1,n-3)];
    y(1,n) = [t*t*t, t*t, t, 1] * h * temp';
    end
end    
plot(x,'b');  
hold on;
plot(smooth(x),'g');  
hold on;
plot(y,'r');   
hold off;