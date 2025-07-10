%Superspreader Deterministic Simulations (with virulence costs)

% ----------------------------------------------
%   ONCE RUN, MOVED SAVED DATA FILE TO
%       "FIGURES_ESS_TRAJECTORY"
% ----------------------------------------------

% we now split the population into n_v vulnerability classes instead of n_c
% connectivity classes

close all
clearvars
rng("shuffle")

%1. Change k = 0.1, 1, 10 and change the save command (at the end) appropriately
%2. The function f(c) relates to figure 4c in the main manuscript 
%3 i. Check tf = 1000 to represent the time to mutation event
%  ii. Change n_iter to run for sufficient time 

%% (1) determine p(i) and c(i) for the host based on the imposed  truncated gamma distribution

mu = 10;        %mean of the gamma distribution
k = 0.2;        %shape parameter of the gamma distribution
theta = mu/k;   %scale parameter of the gamma distribution
nc = 100;       %number of categories for the connectivity

gam_fun1 = @(x) gampdf(x,k,theta);       %function denoting the gamma ditribution pdf
gam_fun2 =  @(x) x.*gampdf(x,k,theta);  
mu_1 = integral(gam_fun2, 0, nc)/integral(gam_fun1, 0, nc);

% finding theta to get the correct mean of mu = 10;
while abs(mu_1 - mu) > 0.0001
    theta = theta*(mu/mu_1);   %scale parameter of the gamma distribution
    gam_fun1 = @(x) gampdf(x,k,theta);       %function denoting the gamma ditribution pdf
    gam_fun2 =  @(x) x.*gampdf(x,k,theta);
    mu_1 = integral(gam_fun2, 0, nc)/integral(gam_fun1, 0, nc);
end

p = zeros(1, nc);   %vector denoting the probabilities of each connectivity
c = zeros(1, nc);   %vector denoting the connectivities

for i = 1:nc
    p(i) = integral(gam_fun1, i-1, i)/integral(gam_fun1, 0, nc);

    if (p(i) > 0) 
        c(i) = integral(gam_fun2, i-1, i)/integral(gam_fun1, i-1, i);
    else
        c(i) = (2*i-1)/2;
    end
end

%% (2) set parameters for the ODE framework. Including the alpha and beta values along a trade-off

b = 10;         K = 1000;       d = 1;      q=(b - d)/K;
al_min = 0;     al_max = 10; 
alpha = (0:100)/10;
n_alpha = length(alpha);      
beta = trade_off_function(alpha, al_min, al_max);

%%% Tolerant function (see Figure 4c)
Mx = 2.0;       Mn = 0.25;      cbar = mu;
A0 = Mx;        A1 = A0-Mn;     A2 = cbar^2*(A0-A1-1)/(1-A0);
func = A0-A1*c.^2./(A2+c.^2);

%% (3) Set up the routine to call the differential equation solver

tf = 1000;
tspan = 0:tf/2:tf;
n_iter = 150;

x0 = p*151;             %inital susceptible population distributed
                        %according to gamma distribution but scaled to have
                        %a total population equal to endemic SS at
                        %k=infinity


alpha_0 = 7;                        %initial starting virulence
y0 = zeros(n_alpha, nc);            %setting up all infected classes    
y0(alpha==alpha_0, :) = p*584;   
        %setting the number of infected individuals of particular virulence
        %alpha_0 to follow the gamma distribution but scaled to have a 
        %total population equal to endemic SS at k=infinity

% z0 = [x0, reshape(y0,1,[])]';
I = y0;

%% (4) Running with mutation, keeping track of total I density of each virulence class

% Want to save:
% 1. Our final vector of Susceptibles, S, for each simulation
% 2. Our final matrix of Infecteds, I, for each simulation
% 3. The vector for alpha
% 4. The density of infecteds in each alpha type (sum over connectivities)
% for each time period, and over each simulation

n_simulations = 10;
% 1.
S_vec = zeros(1, nc, n_simulations);

% 2.
I_mat = zeros(n_alpha, nc, n_simulations);

% 4.
I_vec = zeros(n_alpha, n_iter, n_simulations);
I_vec(:, 1, 1:n_simulations) = repmat(sum(I,2), [1, n_simulations]);

options = odeset('NonNegative',1:nc*(1+n_alpha));

parfor j = 1:n_simulations

    z0 = [x0, reshape(y0,1,[])]';

    for i = 2:n_iter

        [~,z] = ode45(@(t,z) eqn_ode_ss(t, z, b, d, q, alpha, beta, c, p, func), tspan, z0, options);
    
        S = z(end, 1:nc);
        y = z(end, nc+1:end);
        I = reshape(y, n_alpha, nc);
    
        I_vec(:, i, j) = sum(I,2);
    
        % mutation: 1 - choose a random number, 2 - use this random number to
        % determine which class of individuals evolves (density>0), 3 - evolve
        % that particular virulence to either a lower, or higher virulence
        % level (with equal prob.) at 10% density of the previous individual.
    
        u1 = rand;      % random number
        
        mut_vec = [0; cumsum(I_vec(:, i, j))/(sum(I_vec(:, i, j)))];    % cumul. sum of all possible individuals that could mutate
        f1 = find(mut_vec < u1);                        % find which individual mutates
        ind = f1(end);                                  % "
        mid = (mut_vec(ind) + mut_vec(ind+1))/2;        % find the mid point in the prob. space for that ind.
    
        % mutate above or below. Stipulations: 1. if alpha chosen is at 
        %       maximum it can only evolve to lower virulence. 2. if alpha 
        %       chosen is at minimum it can only evolve to higher virulence
        if u1 <= mid
            if rem(ind, n_alpha) ~= 1     % only choose lower if alpha is not minimum
                I(ind-1,:) = I(ind-1,:) + 0.1*I(ind,:);
    
                I(ind,:) = (1 - 0.1)*I(ind,:);
            end 
        end
        if u1 > mid
            if rem(ind, n_alpha) ~= 0     % only choose higher if alpha is not maximum 
                I(ind+1,:) = I(ind+1,:) + 0.1*I(ind,:);
                I(ind,:) = (1 - 0.1)*I(ind,:);
            end
        end
        % mutation complete
        
        z0 = [S, reshape(I,1,[])]';
        t = [];
        z = [];
    
        if i == n_iter
            S_vec(1,:, j) = S;
            I_mat(:,:, j) = I;
        end
    end
end    
clearvars -except S_vec I_mat alpha I_vec

save('k_02_1000_tol.mat')
