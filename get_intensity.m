function iarr = get_intensity(rpos,t,padding)


    % // Default parameters
    param.meshsize  = 128 ;     %// main grid size
    param.patchsize = 200 ;     
    param.windSpeed = 100  ;    %// what unit ? [m/s] ??
    % 注意，初始的 windr 参数现在用一个标量值表示，该值代表风的“方位角”(以度数为单位)(从0到360)。
    param.winddir   = 90   ;    %// Azimuth
    param.rng = 1 ;            %// setting seed for random numbers
    param.A         = 3 ;   %// Scaling factor
%     param.A         = 1e-10 ;   %// Scaling factor
    param.A         = 1e-1 ;   %// Scaling factor
    param.g         = 9.81 ;    %// gravitational constant
    
    param.xLim = [-2 2] ;     %// domain limits X
    param.yLim = [-2 2] ;     %// domain limits Y
    param.zLim = [-1e-4 1e-4]*2 ;
    
    gridSize = param.meshsize * [1 1] ;
    
    %  // Define the grid X-Y domain
    x = linspace( param.xLim(1) , param.xLim(2) , param.meshsize ) ;
    y = linspace( param.yLim(1) , param.yLim(2) , param.meshsize ) ;
    [X,Y] = meshgrid(x, y);
    
    %  // get the grid parameters which remain constants (not time dependent)
    [H0, W, Kx,Ky, Grid_Sign] =  initialize_wave( param ) ;
    
    % // calculate wave at t
    Z = calc_wave( H0 , W , t , Grid_Sign) ;
    % [Zx,Zy] = calc_wave_slope(H0,W,t0,Grid_Sign,Kx,Ky);
    padding_x =  abs( param.xLim(2) - param.xLim(1))/ param.meshsize ;
    padding_y = abs(param.yLim(1)- param.yLim(2) )/ param.meshsize;
    [gx,gy] = gradient(Z,padding_x,padding_y);
    
    % 需要两次独立的测量,rpos是一个二维的
%     两次测量的间隔距离越大，独立性越好，EKF算法也更容易收敛
    offset = [
        0 0; 
        0 padding;
        ];
    iarr = zeros(1,length(offset));
    INT_MIN = 10^(-10);
    ha = 1;
    hw = 1;

for i = 1:length(offset)
    pos = offset(i,:)+rpos;
%     pos
    [ipos3,err] = find_incident_point(pos,X,Y,ha,hw,Z,gx,gy);
%     ipos3
    if err == inf
        iarr(i) = INT_MIN;
        continue;
    end
    rpos3 = [pos ha+hw];
    iarr(i) =  cal_intensity(rpos3,ipos3);
end
