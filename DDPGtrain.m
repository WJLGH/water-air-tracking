%% Initialize environment
clc;
close all;
rng(0);

% generate_circle_traj
generate_trajectory
trajectory = trajectory';
Ts = 0.5;
Tf = 20;
maxEpisodes = 10000;
maxSteps = floor(Tf/Ts);  
env = UAVMovementEnv(trajectory,lv,Ts);

%% create DDPG Agent
obsInfo = env.getObservationInfo;
numObs = obsInfo.Dimension(1);
actInfo = env.getActionInfo;
numAct =  actInfo.Dimension(1);
% load('rlUAVMovementAgent.mat','agent')
createDDPGNetworks
% plot(criticNetwork)
% plot(actorNetwork)
createDDPGAgent
trainOpts = rlTrainingOptions(...
    'MaxEpisodes',maxEpisodes,...
    'MaxStepsPerEpisode',maxSteps,...
    'ScoreAveragingWindowLength',250,...
    'Verbose',true,...
    'Plots','training-progress',...
    'StopTrainingCriteria','AverageReward',...
    'StopTrainingValue',1e4,...                   
    'SaveAgentCriteria','EpisodeReward',... 
    'SaveAgentValue',10);   
%% train 
doTraining = true;
if doTraining    
    % Train the agent.
    trainingStats = train(agent,env,trainOpts);
else
    % Load a pretrained agent for the example.
    load('rlUAVMovementAgent.mat','agent')
end
%% validate
load('ddpg_linear_8v.mat','agent')

generate_trajectory
trajectory = trajectory';
Ts = 0.1;
Tf = 20;
maxEpisodes = 10000;
maxSteps = floor(Tf/Ts);  
env = UAVMovementEnv(trajectory,lv,Ts);
simOptions = rlSimulationOptions('MaxSteps',maxSteps);
experience = sim(env,agent,simOptions);
figure
grid on
for i = 1:length(env.HisPos)
    scatter(env.HisPos(i,1),env.HisPos(i,2),'g');
    hold on
    scatter(env.HisRef(i,1),env.HisRef(i,2),'r');
    pause(0.005);
end
legend('UAV','AUV')

% save('rlUAVMovementAgent.mat','agent')