%% Plotting the variance and speed of evolution

% --------------------------------------
% REQUIRES THE FOLLOWING DATA FROM
% "DETERMINISTIC_WITH_WITHOUT_BENEFITS/NO BENEFITS":
%   1. "K_02_1000.MAT",
%   2. "K_1_1000.MAT"
%   3. "K_10_1000.MAT"
% --------------------------------------

close all
clearvars

%% Load Files In
load('k_02_1000');        k02 = I_vec;
load('k_1_1000');         k1 = I_vec;
load('k_10_1000');        k10 = I_vec;


%% Produce Figure 2a

x= 0:0.1:10;
y= 1:50;

fig2 = figure;
subplot(2,3,1)
plot(x*(mean(k02(:,:,:),3))./(sum(mean(k02(:,:,:),3))), y, 'k-', 'LineWidth', 1.5)
yticklabels({'0', '25', '50'})
yticks([0 25 50])
ylim([0 50.25])
ylabel('Evolutionary time steps')
xlim([2 8])
xticklabels('')
xticks([0 2 4 6 8 10])
%xlabel('Pathogen virulence, \alpha')
title('','k = 0.2')
set(gca,'box','off')
set(gca, 'position', [0.15 0.47 0.25 0.35])
grid on
grid minor
ax = gca;
ax.FontSize = 11;
ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;

subplot(2,3,2)
plot(x*(mean(k1(:,:,:),3))./(sum(mean(k1(:,:,:),3))), y, 'k-', 'LineWidth', 1.5)
yticks([0 25 50])
ylim([0 50.25])
yticklabels('')
xticklabels('')
xlim([2 8])
xticks([0 2 4 6 8 10])
title('','k = 1')
set(gca,'box','off')
set(gca, 'position', [0.42 0.47 0.25 0.35])
grid on
grid minor
ax = gca;
ax.FontSize = 11;
ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;

subplot(2,3,3)
plot(x*(mean(k10(:,:,:),3))./(sum(mean(k10(:,:,:),3))), y, 'k-', 'LineWidth', 1.5)
yticks([0 25 50])
ylim([0 50.25])
yticklabels('')
xticklabels('')
xlim([2 8])
xticks([0 2 4 6 8 10])
title('','k = 10')
set(gca,'box','off')
set(gca, 'position', [0.69 0.47 0.25 0.35])
grid on
grid minor
ax = gca;
ax.FontSize = 11;
ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;