classdef WaveEnv
    properties
        meshsize  = 128 ;     %// main grid size
        patchsize = 200 ;     
        windSpeed = 100  ;    %// what unit ? [m/s] ??
        % 注意，初始的 windr 参数现在用一个标量值表示，该值代表风的“方位角”(以度数为单位)(从0到360)。
        winddir   = 90   ;    %// Azimuth
        rng = 1 ;            %// setting seed for random numbers
%         A         = 0.3 ;   %// Scaling factor
%        A         = 1e-10 ;   %// Scaling factor
        A         = 1e-1 ;   %// Scaling factor
        g         = 9.81 ;    %// gravitational constant
        
        xLim = [-2 2] ;     %// domain limits X
        yLim = [-2 2] ;     %// domain limits Y
        zLim = [-1e-4 1e-4]*2 ;
        
        gridSize = 0;
        X = 0;
        Y = 0;
        H0 = 0;
        Grid_Sign = 0;
        W = 0;
        x_padding = 0;
        y_padding = 0;
        
    end
    methods 
       function this = WaveEnv()
            this.gridSize = this.meshsize * [1 1] ;
            this.x_padding = abs(this.xLim(1) - this.xLim(2))/(this.meshsize-1);
            this.y_padding = abs(this.yLim(1) - this.yLim(2))/(this.meshsize-1);
            x = linspace( this.xLim(1) , this.xLim(2) , this.meshsize ) ;
            y = linspace( this.yLim(1) , this.yLim(2) , this.meshsize ) ;
            [X,Y] = meshgrid(x, y);
            this.X = X;
            this.Y = Y;
            [H0, W, Kx,Ky, Grid_Sign] =  initialize_wave( this ) ;
            this.H0 = H0;
            this.W = W;
            this.Grid_Sign = Grid_Sign;

        end
         function Z = calc_wave(this,time)
            wt = exp(1i .* this.W .* time ) ;
            Ht = this.H0 .* wt + conj(rot90(this.H0,2)) .* conj(wt) ;  
            Z = real( ifft2(Ht) .* this.Grid_Sign ) ;
         end
         function intensity = getIntensityWithAlign(this,rpos,time)
            Z = calc_wave(this,time);             
            [gx,gy] = gradient(Z,this.x_padding,this.y_padding);
            INT_MIN = 10^(-10);
            ha = 1;
            hw = 1;
            pos = rpos;
            [interPos3,err] = find_incident_point(pos,this.X,this.Y,ha,hw,Z,gx,gy);
            if err == inf 
                intensity = INT_MIN;
            else
%                     interPos3
                rpos3 = [pos ha+hw];
                intensity =  cal_intensity_aligned(rpos3,interPos3);
            end
         end
         function intensity = get_intensity(this,rpos,n,time)
                Z = calc_wave(this,time);             
                [gx,gy] = gradient(Z,this.x_padding,this.y_padding);
                INT_MIN = 10^(-10);
                ha = 1;
                hw = 1;
                pos = rpos;
                [interPos3,err] = find_incident_point(pos,this.X,this.Y,ha,hw,Z,gx,gy);
                if err == inf 
                    intensity = INT_MIN;
                else
%                     interPos3
                    rpos3 = [pos ha+hw];
                    intensity =  cal_intensity(rpos3,interPos3,n);
                end
         end

         function intensity = getIntensityByYP(this,x,y,yaw,pitch,time)
                R_xyz= rotz(yaw)*rotx(pitch); 
                n = R_xyz*[0 0 -1]';
                rpos = [x, y];
%                 n
                intensity = get_intensity(this,rpos,n,time);
         end
    end
   
end