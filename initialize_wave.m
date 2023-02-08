function [H0, W,Kx,Ky, Grid_Sign] =  initialize_wave( param )
% function [H0, W, Grid_Sign] =  initialize_wave( param )
%
% This function return the wave height coefficients H0 and W for the
% parameters given in input. These coefficients are constants for a given
% set of input parameters.
% Third output parameter is optional (easy to recalculate anyway)
% ��һ���������� initialize _ wave��M �����������ж����Ժ󶼻��ǳ���(�����Ժ�ʹ�������ж���Ч��ʱ��
% ��������ʱ����仯) �����԰����ŵ�һ��������������ġ�
rng(param.rng);  %// setting seed for random numbers

gridSize = param.meshsize * [1 1] ;

meshLim = pi * param.meshsize / param.patchsize ;
N = linspace(-meshLim , meshLim , param.meshsize ) ;
M = linspace(-meshLim , meshLim , param.meshsize ) ;
[Kx,Ky] = meshgrid(N,M) ;

K = sqrt(Kx.^2 + Ky.^2);    %// ||K||
W = sqrt(K .* param.g);     %// deep water frequencies (empirical parameter)

% �ú������ڰѼ����꣨�����꣩ת��Ϊ�ѿ������ꡣ
% ��ʼ�� windr ����������һ������ֵ��ʾ����ֵ�����ġ���λ�ǡ�(�Զ���Ϊ��λ)(��0��360)
% ���ں��� pol2cart����������ת��Ϊ x �� y �������֤�����շ���������ģ����1
[windx , windy] = pol2cart( deg2rad(param.winddir) , 1) ;
P = phillips(Kx, Ky, [windx , windy], param.windSpeed, param.A, param.g) ;
H0 = 1/sqrt(2) .* (randn(gridSize) + 1i .* randn(gridSize)) .* sqrt(P); % height field at time t = 0
% ����������������
if nargout == 5
    Grid_Sign = signGrid( param.meshsize ) ;
end