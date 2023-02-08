%% // clear workspace
clear all; close all; clc;

%% // Default parameters
param.meshsize  = 128 ;     %// main grid size
param.patchsize = 200 ;     
param.windSpeed = 100  ;    %// what unit ? [m/s] ??
% 注意，初始的 windr 参数现在用一个标量值表示，该值代表风的“方位角”(以度数为单位)(从0到360)。
param.winddir   = 90   ;    %// Azimuth
param.rng = 1 ;            %// setting seed for random numbers
param.A         = 3 ;   %// Scaling factor
param.g         = 9.81 ;    %// gravitational constant

param.xLim = [-2 2] ;     %// domain limits X
param.yLim = [-2 2] ;     %// domain limits Y
param.zLim = [-1e-4 1e-4]*2 ;

gridSize = param.meshsize * [1 1] ;

%% // Define the grid X-Y domain
x = linspace( param.xLim(1) , param.xLim(2) , param.meshsize ) ;
y = linspace( param.yLim(1) , param.yLim(2) , param.meshsize ) ;
[X,Y] = meshgrid(x, y);

%% // get the grid parameters which remain constants (not time dependent)
[H0, W, Kx,Ky, Grid_Sign] =  initialize_wave( param ) ;

%% // calculate wave at t0
t0 = 0;
Z = calc_wave( H0 , W , t0 , Grid_Sign) ;
% [Zx,Zy] = calc_wave_slope(H0,W,t0,Grid_Sign,Kx,Ky);
padding_x =  abs( param.xLim(2) - param.xLim(1))/ param.meshsize ;
padding_y = abs(param.yLim(1)- param.yLim(2) )/ param.meshsize;
[gx,gy] = gradient(Z,padding_x,padding_y);
surf(X,Y,Z);
xlabel('x(m)');
ylabel('y(m)');
zlabel('z(m)');
colorbar;

figure
contour(X,Y,Z)
hold on
quiver(X,Y,gx,gy)
hold off


% t(63:68,63:68)
% Zx(63:68,63:68)
INT_MIN = 10^(-10);

ha = 1;
hw = 1;
intensity = zeros(size(Z));
for i=1:length(x)
    for j =1:length(y)
        rpos = [x(i),y(j)];
        [ipos3,err] = find_incident_point(rpos,X,Y,ha,hw,Z,gx,gy);
        if err == inf
            intensity(i,j) = INT_MIN;
            continue;
        end
        rpos3 = [rpos ha+hw];
        intensity(i,j) = cal_intensity(rpos3,ipos3);
    end
end
% rpos = [0,1];
% [ipos3,err] = find_incident_point(rpos,X,Y,ha,hw,Z,gx,gy);
% rpos3 = [rpos ha+hw];
% intensity = cal_intensity(rpos3,ipos3);


figure();
surf(X,Y,intensity);
% hspan  = 1;
% xlim([-hspan,hspan]);
% ylim([-hspan,hspan]);
xlabel('x(m)');
ylabel('y(m)');
zlabel('intensity(w)');
colorbar;
% %% // populate the display panel
% h.fig  = figure('Color','w') ;
% h.ax   = handle(axes) ;                 %// create an empty axes that fills the figure
% h.surf = handle( surf( NaN(2) ) ) ;     %// create an empty "surface" object
% 
% %% // Display the initial wave surface
% set( h.surf , 'XData',X , 'YData',Y , 'ZData',Z )
% set( h.ax   , 'XLim',param.xLim , 'YLim',param.yLim , 'ZLim',param.zLim )
% 
% %% // Change some rendering options
% % axis off                                %// make the axis grid and border invisible
% shading interp                          %// improve shading (remove "faceted" effect) 对曲面或图形对象的颜色着色进行色彩的插值处理，使色彩平滑过渡
% blue = linspace(0.4, 1.0, 25).' ; cmap = [blue*0, blue*0, blue]; %'// create blue colormap
% colormap(cmap)
% %// configure lighting
% h.light_handle = lightangle(-45,30) ;   %// add a light source
% set(h.surf,'FaceLighting','phong','AmbientStrength',.3,'DiffuseStrength',.8,'SpecularStrength',.9,'SpecularExponent',25,'BackFaceLighting','unlit')
% 
% %% // Animate
% view(75,55) %// no need to reset the view inside the loop ;)
% view(az,el)函数的功能是为三维空间图形设置观察点的方向角

% timeStep = 1./25 ;
% nSteps = 2000 ;
% for time = (1:nSteps)*timeStep    
%     %// update wave surface
%     Z = calc_wave( H0,W,time,Grid_Sign ) ;
%     h.surf.ZData = Z ;
%     pause(0.001);
% end


%% // This block of code is only if you want to generate a GIF file
%// be carefull on how many frames you put there, the size of the GIF can
%// quickly grow out of proportion ;)
% 
% nFrame = 55 ;
% gifFileName = 'MyDancingWaves.gif' ;
% 
% view(-70,40)
% clear im
% f = getframe;
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,20) = 0;
% iframe = 0 ;
% for time = (1:nFrame)*.5
%     %// update wave surface
%     Z = calc_wave( H0,W,time,Grid_Sign ) ;
%     h.surf.ZData = Z ;
%     pause(0.001);
% 
%     f = getframe;
%     iframe= iframe+1 ;
%     im(:,:,1,iframe) = rgb2ind(f.cdata,map,'nodither');
% end
% imwrite(im,map,gifFileName,'DelayTime',0,'LoopCount',inf)
% disp([num2str(nFrame) ' frames written in file: ' gifFileName])