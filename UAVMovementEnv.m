classdef UAVMovementEnv < rl.env.MATLABEnvironment
    %LIGHTINTENSITYENV: Template for defining custom environment in MATLAB.    
    
    %% Properties (set properties' attributes accordingly)
    properties
       %X Y 坐标的最大值
        xMax = 10;
        yMax = 10;
        %UAV的初始位置和目标结束位置
        startX = 0;
        startY = 0;
        endX = 0;
        endY = 0;
        %AUV的初始位置
        refX = 0;
        refY = 0;
        refV = 0;
        refHeadingAngle = 0;
        %当前时间
        t = 0;
        Tf = 20;
        % 采样时间
        Ts = 0.1;
        %当前UAV x坐标
        uX = 0;
        %当前UAV y坐标
        uY = 0;
        %当前UAV 速度
        uV = 0;
        %当前UAV 航向角
        uHeadingAngle=0;
        % 控制量的最大值  
        dhMax = 180; %最大转动角度
        dvMax = 3; %最多加减速
        % UAV和AUV之间的相对距离
        preDis = 0;
        count = 0;
        INT_THRESHOLD = 1e-10;
        wave = 0;
        % AUV的轨迹
        trajectory = 0;
        sim = false;
        uIntialV = 3;
        uIntialH = 90;
        HisPos = [];
        HisRef = [];
    end
    
    properties
        % 初始化系统状态 
        State = zeros(8,1)
    end
    
    properties(Access = protected)
        % Initialize internal flag to indicate episode termination
        IsDone = false        
    end
%     1.8474e-7，最大值

    %% Necessary Methods
    methods              
        function this = UAVMovementEnv(trajectory,refV,ts)
            
           % Initialize Observation settings.
%             refPos,refV,refHeadingAngle
%             uPos,uV,uHeadingAngle
            numObs = 8;
            ObservationInfo = rlNumericSpec([numObs 1]);
            ObservationInfo.Name = 'observations';
            %给出造成结果的原因
            ObservationInfo.Description = 'refX,refY,refV,refHeadingAngle,uX,uY,uV,uHeadingAngle';
            
            % Initialize Action settings  
%             numAct = 9;         
%             ActionInfo = rlFiniteSetSpec(num2cell(1:numAct));
            % 每次移动角度大小
            numActions = 2;
            ActionInfo = rlNumericSpec([numActions 1],'LowerLimit',-1,'UpperLimit', 1);
            ActionInfo.Name = 'action';
            % 俯仰角和转向角的增量
            ActionInfo.Description = 'deltaV, deltaHeadingAngle';
            % The following line implements built-in functions of RL env
            this = this@rl.env.MATLABEnvironment(ObservationInfo,ActionInfo);
            this.refX = trajectory(1,1);
            this.refY = trajectory(2,1);
            this.refHeadingAngle = trajectory(3,1);
            this.trajectory = trajectory;
            this.refV = refV;
            this.Ts = ts;
            this.wave = WaveEnv();
            this.startX = trajectory(1,1)+rand()*0.1;
            this.startY = trajectory(2,1)+rand()*0.1;
            this.endX = trajectory(1,end);
            this.endY = trajectory(2,end);
            
        end
        % Apply system dynamics and simulates the environment with the 
        % given action for one step. action is 2x1 matrix
        function [Observation,Reward,IsDone,LoggedSignals] = step(this,action)
            LoggedSignals = [];
%             Observation = [this.refX,this.refY,this.refV,this.refHeadingAngle,this.uX,this.uY,this.uV,this.uHeadingAngle]; %#ok<*PROPLC> 
%             Reward = -100;
%             IsDone = true;
%             return

            act = satufun(this,action);

            dv = act(1);
            dh = act(2);
            
            
            this.uV = this.uV + dv;
            if this.uV < 0
                this.uV = 0;
            end

%             this.uX = this.uX + this.Ts* ( cosd(this.uHeadingAngle) - this.uV*sind(this.uHeadingAngle) );
%             this.uY = this.uY + this.Ts*( sind(this.uHeadingAngle) + this.uV * cosd(this.uHeadingAngle) );
            
            this.uHeadingAngle = this.uHeadingAngle + dh;
            this.uHeadingAngle = mod(this.uHeadingAngle,360);
            this.uX = this.uX + this.Ts* ( cosd(this.uHeadingAngle) );
            this.uY = this.uY + this.Ts*( sind(this.uHeadingAngle) );
            

           
           

            this.t = this.t + this.Ts;
            this.count  = this.count +1;
            idx = ceil(this.t/this.Tf*size(this.trajectory,2));



            % Check terminal condition
            
            IsDone = true;
%             out = outOfRange(this);
%             if out
%                 Reward = -300;
            if norm([this.uX this.uY] - [this.endX this.endY]) < 0.1
                Reward = this.count;
            elseif idx >= size(this.trajectory,2)                
                Reward = 0;
            else
                IsDone = false;
            end

            if ~IsDone                
                refX = this.trajectory(1,idx);
                refY = this.trajectory(2,idx);
                refHeadingAngle = this.trajectory(3,idx);
                refV = this.refV;
                this.refX = refX;
                this.refY = refY;
                this.refV = refV;
                this.refHeadingAngle = refHeadingAngle;
            end
            
            Observation = [this.refX,this.refY ,this.refV,this.refHeadingAngle,this.uX,this.uY,this.uV,this.uHeadingAngle]; %#ok<*PROPLC> 
            this.State = Observation;
           
            if ~IsDone                
                Reward = getReward(this,dv,dh);                
            end

            this.HisRef = [this.HisRef;this.refX, this.refY];
            this.HisPos = [this.HisPos;this.uX,this.uY]; 
            this.IsDone = IsDone;
            notifyEnvUpdated(this);
        end
        function done = outOfRange(this)
            done = abs(this.uX)> this.xMax || abs( this.uY) > this.yMax;
        end
        % Reset environment to initial state and output initial observation
            function Observation = reset(this)
%             "reset env"   ，测试每次是否调用reset函数
            this.IsDone = false;
            this.t = 0;
            this.count = 1;
            this.uX = this.startX;
            this.uY = this.startY;
            this.uV = this.uIntialV;
            this.uHeadingAngle = this.uIntialH;
            this.refX = this.trajectory(1,1);
            this.refY = this.trajectory(2,1);
            this.refHeadingAngle = this.trajectory(3,1);
            this.preDis = norm([this.uX this.uY] - [this.refX this.refY]);
            Observation = [this.refX,this.refY,this.refV,this.refHeadingAngle, ...
                this.uX,this.uY,this.uV,this.uHeadingAngle];
            this.State = Observation;
            this.HisRef = [this.refX, this.refY];
            this.HisPos = [this.uX,this.uY];
            notifyEnvUpdated(this);
        end
    end
    %% Optional Methods (set methods' attributes accordingly)
    methods               
        % Helper methods to create the environment
        % Discrete force 1 or 2
        function act = satufun(this,action)
            dv = action(1)*this.dvMax;
            %将 -1~1 变为 0 -1
            dh = (action(2)+1.0)/2*this.dhMax;
            act = [dv dh];
            v = 0.5*exp(-this.count*1e-3);
            covar = diag([v v]);
            mu = [dv dh];
            act = mvnrnd(mu,covar,1);
        end
        
        % Reward function
        function Reward = getReward(this,dv,dh)
            x = this.uX - this.refX;
            y = this.uY - this.refY;
            intensity = this.wave.getIntensityWithAlign([x y],this.t);
            if intensity < this.INT_THRESHOLD
                rT = -10;
            else            
                dis = norm([this.uX this.uY] - [this.refX this.refY]);
                preDis = this.preDis;
                rT = log(dis/this.preDis);
                this.preDis = dis;

                rT = (preDis - dis)*1/dis+this.t/this.Tf;
%                 rT = - ( dis + 0.1*deg2rad(abs(this.uHeadingAngle - this.refHeadingAngle)) ) ;
                rT = -dis;

            end
            rho = 0.01;
            pE = 1/2*rho*this.uV^3/this.dvMax;
            pE = 0;
            pJ = abs(dh)+abs(dv);
            pJ = 0;
            Reward = 1*rT+ 0.01*pE + 0.01*pJ;            
        end
        % (optional) Visualization method
        function plot(this)
            % Initiate the visualization
            
%             quiver(this.uX,this.uY,cos(this.uHeadingAngle),sin(this.uHeadingAngle),0,'r')
%             plot(this.uX,this.uY)
            scatter(this.refX,this.refY,'red');
            scatter(this.uX,this.uY,'green');
            hold on
%             [this.refX, this.refY]
            % Update the visualization
            envUpdatedCallback(this)
        end
    end
end
