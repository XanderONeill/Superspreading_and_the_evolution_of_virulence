%   Superspreader Deterministic Simulations (with benefits)


% Produces data on how the ESS varies with shape parameter k.
close all
clearvars
rng("shuffle")

% Finding the ESS for k = [0.1:0.3:1, 4:3:10, 40:30:100,400:300:1000];
% Finding the ESS for k = [0.1:0.1:1, 2:1:10,20:10:100,200:100:1000];

%% Fixed parameters and initialising
mu = 10;   %mean of the gamma distribution
nc = 100;   %number of categories for the connectivity

%set parameters for the ODE framework. Including the alpha and beta values along a trade-off
b = 10;         K = 1000;

al_min = 0;     al_max = 10;    alpha = (0:100)/10;
n_alpha = length(alpha); 
beta = trade_off_function(alpha, al_min, al_max);

%initialising:

alpha_0 = 3.5; %initial starting virulence
k_vec = [0.2:0.1:1, 2, 5, 10, 20, 50, 100, 1000];

alpha_vec = zeros(length(k_vec),1);
I_MORT = zeros(length(k_vec), nc);
S_MORT = zeros(length(k_vec), nc);

% time to extinction events
tf = 1000;
tspan = 0:tf/2:tf;

options = odeset('NonNegative',1:nc*(1+n_alpha));

for k = k_vec %shape parameter of the gamma distribution
       
    %% (1) determine p(i) and c(i) for the host based on the imposed truncated gamma distribution.

    theta = mu/k;   %scale parameter of the gamma distribution

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

%% (2) Parameters that will vary with k (our benefits function)
    
    % Benefit function d = d(c)
    Max_1 = 4.0;    Min_1 = 0.25;   cbar = mu;
    A_0 = Max_1;    A_1 = A_0 - Min_1;  A_2 = cbar^2*(A_0 - A_1 - 1)/(1 - A_0);
    d = A_0 - A_1*c.^2./(A_2+c.^2);

    q=(b - d*p')/K;

%% (3) Set up the routine to call the differential equation solver
   
    x0 = p*151;             %inital susceptible population distributed
                            %according to gamma distribution but scaled to have
                            %a total population equal to endemic SS at
                            %k=infinity                      
    y0 = zeros(n_alpha, nc);            %setting up all infected classes    
    y0(alpha==alpha_0, :) = p*584;   
            %setting the number of infected individuals of particular virulence
            %alpha_0 to follow the gamma distribution but scaled to have a 
            %total population equal to endemic SS at k=infinity
    
    z0 = [x0, reshape(y0,1,[])]';
    I = y0;

%% (4) Running with mutation, keeping track of total I density of each virulence class

    I_vec = zeros(n_alpha,1);
    I_vec(:,1) = sum(I,2);
    
    for i = 2:50
        [~,z] = ode45(@(t,z) eqn_ode_ss(t, z, b, d, q, alpha, beta, c, p), tspan, z0, options);
        
        S = z(end, 1:nc);
        y = z(end, nc+1:end);
        I = reshape(y, n_alpha, nc);
    
        I_vec(:,i) = sum(I,2);
    
        % mutation: 1 - choose a random number, 2 - use this random number to
        % determine which class of individuals evolves (density>0), 3 - evolve
        % that particular virulence to either a lower, or higher virulence
        % level (with equal prob.) at 10% density of the previous individual.
    
        u1 = rand;      % random number
        
        mut_vec = [0; cumsum(I_vec(:, i))/(sum(I_vec(:, i)))];    % cumul. sum of all possible individuals that could mutate
        f1 = find(mut_vec < u1);                        % find which individual mutates
        ind = f1(end);                                  % "
        mid = (mut_vec(ind) + mut_vec(ind+1))/2;        % find the mid point in the prob. space for that ind.
    
        % mutate above or below. Stipulations: 1. if alpha chosen is at 
        %       maximum it can only evolve to lower virulence. 2. if alpha 
        %       chosen is at minimum it can only evolve to higher virulence
        if u1 < mid
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
        z = [];
    end

    [~,z] = ode45(@(t,z) eqn_ode_ss(t, z, b, d, q, alpha, beta, c, p), tspan, z0, options);
    S = z(end, 1:nc);
    y = z(end, nc+1:end);
    I = reshape(y, n_alpha, nc);
    
    i = i+1;
    I_vec(:,i) = sum(I,2);

    alpha_mean = alpha*I_vec(:,end)/sum(I_vec(:,end));
    alpha_0 = min(alpha(alpha>alpha_mean)); %initial starting virulence for next round of k
    f = find(k == k_vec);
    alpha_vec(f) = alpha*I_vec(:,end)/sum(I_vec(:,end)); 
    I_MORT(f, :) = sum(I, 1);
    S_MORT(f, :) = S;
end

clearvars -except alpha alpha_vec k_vec I_MORT S_MORT d

save('ESS_k_te_1000')

