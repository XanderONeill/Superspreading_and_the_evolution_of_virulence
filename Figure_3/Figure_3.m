clearvars
close all

% ----------------------------------------------
% REQUIRES NEW DATA FROM:
%       1. "RUN_ESS_VARY_K.M" FROM
%       "DETERMINISTIC_WITH_WITHOUT_BENEFITS/WITH_BENEFITS"
% BEFORE CONTINUING
% ----------------------------------------------

%% Load in Files

load('ESS_k_te_1000')

%% Plot Figures

c=0:1:100;
% Benefit function  d = d(c);
Max_1 = 4.0;    Min_1 = 0.25;   cbar = 10;
A_0 = Max_1;    A_1 = A_0 - Min_1;  A_2 = cbar^2*(A_0 - A_1 - 1)/(1 - A_0);
d = A_0 - A_1*c.^2./(A_2+c.^2);

fig = figure;

subplot(2,1,2)
plot(k_vec, alpha_vec, 'LineWidth', 2)
ylim([2 5])
ylabel('ESS Virulence', 'Color','k')
xlim([0.2 1000])
xlabel('Shape parameter, k')
ax = gca;
ax.FontSize = 11;
%set(gca,'box','off')
set(gca, 'XScale', 'log')
set(gca, 'position', [0.1 0.13 0.86 0.45])

subplot(2,1,1)
plot(c, d, 'LineWidth', 2);
ylabel('Natural death, d(c_i)')
xlabel('Individual transmission, c_i')
ax = gca;
ax.FontSize = 8;
%set(gca,'box','off')
set(gca, 'position', [0.69 0.23 0.25 0.25])

