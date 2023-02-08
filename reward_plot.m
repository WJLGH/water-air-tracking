load('training_reward.mat')
len = length(ave_episode_reward);
idx = 1:4:1400;
plot(idx,episode_reward(idx),'-o')
hold on
idx = 1:10:1400;
plot(idx,ave_episode_reward(idx),'-*','linewidth',0.6)
xlim([0,1400])
xlabel('Episode')
ylabel('Reward')
grid on
legend('Epsiode Reward','Average Reward')