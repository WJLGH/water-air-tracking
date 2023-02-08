
param.meshsize  = 128 ;     %// main grid size
param.patchsize = 200 ;     
param.meshsize  = 8 ;     %// main grid size
param.patchsize = 8;    
param.windSpeed = 100  ;    %// what unit ? [m/s] ??
% 注意，初始的 windr 参数现在用一个标量值表示，该值代表风的“方位角”(以度数为单位)(从0到360)。
param.winddir   = 90   ;    %// Azimuth
param.rng = 13 ;            %// setting seed for random numbers
param.A         = 1e-7 ;    %// Scaling factor
param.g         = 9.81 ;    %// gravitational constant

param.span = 10;
param.xLim = [-param.span param.span] ;     %// domain limits X
param.yLim = [-param.span param.span] ;     %// domain limits Y
param.zLim = [-1e-4 1e-4]*2 ;

gridSize = param.meshsize * [1 1] ;


x = linspace( param.xLim(1) , param.xLim(2) , param.meshsize ) ;
y = linspace( param.yLim(1) , param.yLim(2) , param.meshsize ) ;
[X,Y] = meshgrid(x, y);


[H0, W, K, Grid_Sign] =  initialize_wave( param ) ;

t0 = 0 ;
Z = calc_wave( H0 , W , t0 , Grid_Sign ) ;


timeStep = 1./25 ;
nSteps = 2000 ;
padding = param.span*2/param.meshsize;
idx = param.meshsize/2;
dx1= (Z(idx,idx)-Z(idx-1,idx))/padding;
dx2 = (Z(idx+1,idx)-Z(idx,idx))/padding;
slop = calc_slope( H0 , W , K,t0 , Grid_Sign ) ;
dx1
dx2
slop(idx,idx)