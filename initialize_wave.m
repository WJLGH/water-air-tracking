function [H0, W,Kx,Ky, Grid_Sign] =  initialize_wave( param )
% function [H0, W, Grid_Sign] =  initialize_wave( param )
%
% This function return the wave height coefficients H0 and W for the
% parameters given in input. These coefficients are constants for a given
% set of input parameters.
% Third output parameter is optional (easy to recalculate anyway)
% 第一个函数叫做 initialize _ wave。M 这里计算的所有东西以后都会是常量(当您以后使波动具有动画效果时，
% 它不会随时间而变化) ，所以把它放到一个块中是有意义的。
rng(param.rng);  %// setting seed for random numbers

gridSize = param.meshsize * [1 1] ;

meshLim = pi * param.meshsize / param.patchsize ;
N = linspace(-meshLim , meshLim , param.meshsize ) ;
M = linspace(-meshLim , meshLim , param.meshsize ) ;
[Kx,Ky] = meshgrid(N,M) ;

K = sqrt(Kx.^2 + Ky.^2);    %// ||K||
W = sqrt(K .* param.g);     %// deep water frequencies (empirical parameter)

% 该函数用于把极坐标（柱坐标）转换为笛卡尔坐标。
% 初始的 windr 参数现在用一个标量值表示，该值代表风的“方位角”(以度数为单位)(从0到360)
% 由于函数 pol2cart，它后来被转换为 x 和 y 组件。保证了最终风速向量的模总是1
[windx , windy] = pol2cart( deg2rad(param.winddir) , 1) ;
P = phillips(Kx, Ky, [windx , windy], param.windSpeed, param.A, param.g) ;
H0 = 1/sqrt(2) .* (randn(gridSize) + 1i .* randn(gridSize)) .* sqrt(P); % height field at time t = 0
% 如果结果向量有三个
if nargout == 5
    Grid_Sign = signGrid( param.meshsize ) ;
end