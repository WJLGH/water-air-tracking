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

% ������غ��� signCor (ifft2(Ht) ��meshSize) ���򵥵ؽ� Ht ������Ԫ�صķ��ŵߵ���������һ������ķ������Դﵽ���Ŀ��: 
% �򵥵���һ��ͬ����С�ľ���(Grid _ sign)���� Ht����������ǽ���� + 1-1... ... �ȵȡ�