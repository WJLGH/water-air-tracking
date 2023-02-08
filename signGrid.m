function sgn = signGrid(n)
% return a matrix the size of n with alternate sign for every indice
% ex:     sgn = signGrid(3) ;
%         sgn =
%             -1     1    -1
%              1    -1     1
%             -1     1    -1

    [x,y] = meshgrid(1:n,1:n) ;
    sgn = ones( n ) ;
    sgn(mod(x+y,2)==0) = -1 ;
end

% 你的神秘函数 signCor (ifft2(Ht) ，meshSize) ，简单地将 Ht 的其他元素的符号颠倒过来。有一个更快的方法可以达到这个目的: 
% 简单地用一个同样大小的矩阵(Grid _ sign)乘以 Ht，这个矩阵是交替的 + 1-1... ... 等等。