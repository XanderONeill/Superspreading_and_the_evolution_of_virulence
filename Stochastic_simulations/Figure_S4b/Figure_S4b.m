%% Plotting the variance and speed of evolution for the stochastic model

% --------------------------------------
% REQUIRES THE FOLLOWING DATA FROM
% "STOCHASTIC_WITH_VIRULENCE_COSTS/RUN_STOCHASTIC_SIMS":
%   1. "K_02_MARY.MAT"
%   2. "K_1_MARY.MAT"
%   3. "K_10_MARY.MAT"
% --------------------------------------
clearvars
close all

%% Load Data In

load('k_02');       I_02 = Ivec_3D;
load('k_1');        I_1 = Ivec_3D;
load('k_10');       I_10 = Ivec_3D;

%% Finding means and variances

Mean_k02 = zeros(size(I_02, 1), 1);
Mean_k1  = zeros(size(I_1, 1), 1);
Mean_k10 = zeros(size(I_10, 1), 1);

Var_k02 = zeros(size(I_02, 1), 1);
Var_k1  = zeros(size(I_1, 1), 1);
Var_k10 = zeros(size(I_10, 1), 1);

for i = 1:size(I_02,1)
    Mean_k02(i) = alpha*sum(I_02(i, :, :), 3)'/sum(I_02(i, :, :), 'all');
    Mean_k1(i)  = alpha*sum(I_1(i, :, :), 3)'/sum(I_1(i, :, :), 'all');
    Mean_k10(i) = alpha*sum(I_10(i, :, :), 3)'/sum(I_10(i, :, :), 'all');
    Var_k02(i) = sum(I_02(i,:,:), 3)*(alpha - Mean_k02(i)).^2'/sum(I_02(i,:,:), 'all');
    Var_k1(i) = sum(I_1(i,:,:), 3)*(alpha - Mean_k1(i)).^2'/sum(I_1(i,:,:), 'all');
    Var_k10(i) = sum(I_10(i,:,:), 3)*(alpha - Mean_k10(i)).^2'/sum(I_10(i,:,:), 'all');
end

%% Produce Figures

x = [0 10];
y = [10001 0];

fig = figure;

subplot(1,3,1)
clims = ([0.1 max(I_02(end,:,1))]);
imagesc(x, y, sum(I_02(:,:,1:10),3), clims)
hold on
plot(Mean_k02, 10000:-1:0, 'LineWidth', 2)
%plot(Mean_k02 + 3*sqrt(Var_k02), 10000:-1:0, 'Color', [0.8500 0.3250 0.0980])
%plot(Mean_k02 - 3*sqrt(Var_k02), 10000:-1:0, 'Color', [0.8500 0.3250 0.0980])
ylabel('Time, t')
yticklabels({'100000', '80000', '60000','40000','20000','0'})
ylim([0 10000])
yticks([0 2000 4000 6000 8000 10000])
xlim([0 10])
xticks([0 2 4 6 8 10])
xtickangle(0)
set(gca,'box','off')
set(gca, 'position', [0.15 0.47 0.25 0.35])
grid on
grid minor
ax = gca;
ax.FontSize = 11;
ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;

subplot(1,3,2)
clims = ([0.1 max(I_1(end,:,1))]);
imagesc(x, y, sum(I_1(:,:,1:10),3), clims)
hold on
plot(Mean_k1, 10000:-1:0, 'LineWidth', 2)
%plot(Mean_k1 + 3*sqrt(Var_k1), 10000:-1:0, 'Color', [0.8500 0.3250 0.0980])
%plot(Mean_k1 - 3*sqrt(Var_k1), 10000:-1:0, 'Color', [0.8500 0.3250 0.0980])
yticklabels('')
ylim([0 10000])
yticks([0 2000 4000 6000 8000 10000])
xlim([0 10])
xticks([2 4 6 8 10])
xtickangle(0)
xlabel('Pathogen virulence, \alpha')
set(gca,'box','off')
set(gca, 'position', [0.42 0.47 0.25 0.35])
grid on
grid minor
ax = gca;
ax.FontSize = 11;
ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;

subplot(1,3,3)
clims = ([0.1 max(I_10(end,:,1))]);
imagesc(x, y, sum(I_10(:,:,1:10),3), clims)
hold on
plot(Mean_k10, 10000:-1:0, 'LineWidth', 2)
%plot(Mean_k10 + 3*sqrt(Var_k10), 10000:-1:0,'Color', [0.8500 0.3250 0.0980])
%plot(Mean_k10 - 3*sqrt(Var_k10), 10000:-1:0, 'Color', [0.8500 0.3250 0.0980])
yticklabels('')
ylim([0 10000])
yticks([0 2000 4000 6000 8000 10000])
xlim([0 10])
xticks([2 4 6 8 10])
xtickangle(0)
set(gca,'box','off')
set(gca, 'position', [0.69 0.47 0.25 0.35])
grid on
grid minor
ax = gca;
ax.FontSize = 11;
ax.GridColor = [0 .5 .5]; ax.GridLineStyle = '--'; ax.GridAlpha = 0.5;

colormap(gray*(-1) +1);

%% Mean Variance
mean(Var_k02)
mean(Var_k1)
mean(Var_k10)

Var_k02(end)
Var_k1(end)
Var_k10(end)