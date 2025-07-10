%% Plotting the variance and speed of evolution

% --------------------------------------
% REQUIRES THE FOLLOWING DATA FROM
% "DETERMINISTIC_WITH_VIURLENCE_COSTS":
%   1. "K_02_1000_MARY.MAT",
%   2. "K_1_1000_MARY.MAT"
%   3. "K_10_1000_MARY.MAT"
% --------------------------------------

% This script is for simulations with long evolutionary time steps (t_e =
% 1000)

close all
clearvars

%% Load Files In
load('k_02_1000_tol');        k02 = I_vec;
load('k_1_1000_tol');         k1 = I_vec;
load('k_10_1000_tol');        k10 = I_vec;

%% Produce Figure S4a

x= 0:0.1:10;
y= 1:75;

fig = figure;
subplot(2,3,1)
plot(x*(mean(k02(:,:,:),3))./(sum(mean(k02(:,:,:),3))), y, 'k-', 'LineWidth', 1.5)
yticklabels({'0', '25', '50', '75'})
yticks([0 25 50 75])
ylim([0 75])
ylabel('Evolutionary time steps')
xlim([0 10])
xticks([0 2 4 6 8 10])
xticklabels('')
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
yticks([0 25 50 75])
ylim([0 75])
yticklabels('')
xticklabels('')
xlim([0 10])
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
yticks([0 25 50 75])
ylim([0 75])
yticklabels('')
xticklabels('')
xlim([0 10])
xticks([0 2 4 6 8 10])
title('','k = 10')
set(gca,'box','off')
set(gca, 'position', [0.69 0.47 0.25 0.35])
grid on
grid minor
ax = gca;
ax.FontSize = 11;
ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;
