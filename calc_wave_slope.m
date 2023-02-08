function [Zx,Zy] = calc_wave_slope(H0,W,time,Grid_Sign,Kx,Ky)
% Z = calc_wave( H0,W,time,Grid_Sign )
%
% This function calculate the wave height based on the wave coefficients H0
% and W, for a given "time". Default time=0 if not supplied.
% Fourth output parameter is optional (easy to recalculate anyway)

    % recalculate the grid sign if not supplied in input
    if nargin < 5
        Grid_Sign = signGrid( param.meshsize ) ;
    end
    % Assign time=0 if not specified in input
    if nargin < 3 ; time = 0 ; end
% Ht = H0 .* exp(1i .* W .* (t * timeStep)) + \
%     conj(flip(flip(H0,1),2)) .* exp(-1i .* W .* (t * timeStep));
% wt = exp(1i .* W .* time ) ;
% Ht = H0 .* wt + conj(rot90(H0,2)) .* conj(wt) ;  
% Z = real( ifft2(Ht) .* Grid_Sign ) ;
% conj()计算一个复数的共轭
% =rot90（A，k）矩阵A按照逆时针旋转90*k度
    wt = exp(1i .* W .* time ) ;
    Ht = H0 .* wt + conj(rot90(H0,2)) .* conj(wt) ;  
    Hkx = Ht.*1i.*Kx;
    Hky = Ht.*1i.*Ky;
    Zx = real( ifft2(Hkx) .* Grid_Sign ) ;
    Zy = real( ifft2(Hky) .* Grid_Sign ) ;
    n = sqrt(Zx.^2+Zy.^2+1);
%   向量归一化处理
    Zx = Zx./n;
    Zy = Zy./n;
%     进行傅里叶反变换，然后将正负间隔设置？
%     Z = real( ifft2(Ht) .* Grid_Sign ) ;
end
