function intensity = cal_intensity_aligned(rpos,ipos)
        INT_MIN = 10^(-10);
        T_r = rpos - ipos;
        da = norm(T_r);
        dw = norm(ipos);
        a_t = acos(ipos(3)/norm(ipos));
        %对准了a_r_cos = 1；
        a_r_cos = 1;
        % 参数
        
        % 波束发散角
        theta = pi/2;
        % 水中衰减系数
        c = 0.15;
        fov = pi;
        Ar = 0.01 ^ 2;

%         INTENSITY_MIN = 1e-20;
%         if fov/2 < a_r
%             intensity = INTENSITY_MIN;
%             return
%         end


        bw = -log(2) / log(cos(theta / 2));
        tau = 0.9;
        % 接收器的光电效应比
        R = 0.6;
        P_tx = 0.03;
        Is = P_tx * (bw+1)/(2*pi*dw)*cos(a_t)^bw;
        la = 1/(2*theta*da^2);
        lw = exp(-c*dw);

        A_eff = Ar*a_r_cos;
        intensity = R*tau*la*lw*A_eff*Is;
end