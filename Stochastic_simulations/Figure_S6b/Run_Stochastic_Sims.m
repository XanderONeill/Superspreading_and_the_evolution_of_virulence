%%% Semi-stochastic model simulations (with virulence costs)

% The selection of birth, death and transmission (and mutation) events
% are completely random processes. However, we still split the population
% into n_v vulnerability classes (unlike a full stochastic model where the
% vulnerability of an individual is taken from a gamma distribution at
% birth - hence semi-stochastic).

%1. Change k = 0.1, 1, 10 and change the save command (at the end) appropriately
%2. The function f(c) can be amended accordingly, see section (2).
%3. You may need to change tf if the model simulations have not reached
%the ESS by tf. 

%4. n_sims denotes the number of simulations you would like to try, change
% as you please.

clearvars;
close all;
rng("shuffle")

%% (1) determine p(i) and c(i) for the host based on the imposed truncated gamma distribution.

mu = 10;        %mean of the gamma distribution
k = 1;        %shape parameter of the gamma distribution
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

b = 10;         d = 1;          K = 1000;       q=(b - d)/K;

al_min = 0;     al_max = 10;    alpha = (0:100)/10;
n_alpha = length(alpha);
beta = trade_off_function(alpha, al_min, al_max);

%%% Vulnerability function (see Figure 4b)
Mx = 4;     Mn = 0.5;       cbar = 10;
A0 = Mn;    A1 = Mx-A0;     A2 = cbar^2*A1/(1-A0) - cbar^2;
func = A0 + A1*c.^2./(A2+c.^2);

%% (3) Set up the routine to call the differential equation solver

t_f = 100000;
n_sims = 10;
Ivec_3D = zeros(t_f/10+1, n_alpha, n_sims);

parfor i = 1:n_sims
    alpha_0 = 6.5;                        %initial starting virulence

    S_0 = round(151*p);
    I_0 = zeros(n_alpha, nc);
    I_0(alpha==alpha_0, :) = round(584*p);

    S = S_0;
    I = I_0;

    I_vec = zeros(1,  n_alpha);
    I_vec(1, :) = sum(I,2);

    %%%  Gillespie Runs

    % Initial Event Probabilities
    N = sum(S)+sum(I,"all");
    birth = max(0, N*(b-q*N));
    sus_death = sum(S)*d;
    inf_death = sum(d*I, 'all') + alpha*(func*I')';
    trans = (c*S')*(mu*sum(beta*I));
    R_T = birth + sus_death + inf_death + trans;

    j = 1;
    t_int = 0;
    t = 0;

    while t < t_f
        dt = exprnd(1/R_T);
        t = t + dt;
        u1 = rand;
    
        if u1 < birth/R_T
            % Determine which specific event within birth
            prob_vec = [0, cumsum(p)]*birth/R_T;
            f1 = find(prob_vec < u1);
    
            % Event
            S(f1(end)) = S(f1(end)) + 1;
            N = N + 1;
    
            % update event probabilities
            birth = max(0, N*(b - q*N));
            sus_death = sus_death + d;
            trans = (c*S')*(mu*sum(beta*I));
    
        elseif u1 < (birth + sus_death)/R_T
            % Determine specific susceptible to die
            prob_vec = [0, cumsum(d.*S)]/R_T + birth/R_T;
            f1 = find(prob_vec < u1);
    
            % Event
            S(f1(end)) = S(f1(end)) - 1;        
            N = N - 1;
    
            % update event probabilities
            birth = max(0, N*(b - q*N));
            sus_death = sus_death - d;
            trans = (c*S')*(mu*sum(beta*I));
    
        elseif u1 < (birth +sus_death + inf_death)/R_T
            Mat1 = func.*(alpha'.*I) + d*I;
    
            % Determine which vulnerability class
            prob_vec =  [0, cumsum(sum(Mat1,1))]/R_T + (birth+sus_death)/R_T;
            f1 = find(prob_vec < u1);
            vuln_choice = f1(end);
    
            % Determine which virulence class
            prob_vec_2 = [0; cumsum(Mat1(:,vuln_choice)/sum(Mat1(:,vuln_choice)))]*(prob_vec(f1(end)+1) - prob_vec(f1(end)))+ prob_vec(f1(end));
            f2 = find(prob_vec_2 < u1);
            vir_choice = f2(end);
    
            % Event
            I(vir_choice, vuln_choice) = I(vir_choice, vuln_choice) - 1;        
            N = N - 1;
    
            % update event probabilities
            birth = max(0, N*(b - q*N));
            inf_death = inf_death - (d + alpha(vir_choice)*func(vuln_choice));
            trans = (c*S')*(mu*sum(beta*I));
    
        else
            % Determine which specific susceptible becomes infected
            prob_vec = [0, cumsum(c.*S)/(c*S')]*trans/R_T + (birth + sus_death + inf_death)/R_T;
            f1 = find(prob_vec < u1);
            vuln_choice = f1(end);
    
            % Determine which specific infected transmits
            prob_vec_2 = [0, cumsum(mu*beta.*sum(I,2)'/(mu*sum(beta*I)))]*(prob_vec(f1(end)+1) - prob_vec(f1(end))) + prob_vec(f1(end));
            f2 = find(prob_vec_2 < u1);
            vir_choice = f2(end);
    
            % Determine if evolution event occurs
            evol = rand;
            if (evol <= 0.001) && (vir_choice ~= 1)
                vir_choice = vir_choice - 1;
            end
            if (evol > 0.999) && (vir_choice ~= n_alpha)
                vir_choice = vir_choice + 1;
            end
    
            % Event
            S(1, vuln_choice) = S(1, vuln_choice) - 1;
            I(vir_choice, vuln_choice) = I(vir_choice, vuln_choice) + 1;
    
            % update event probabilities
            sus_death = sus_death - d;
            inf_death = inf_death + (d + alpha(vir_choice)*func(vuln_choice)); 
            trans = (c*S')*(mu*sum(beta*I));
        end
        R_T = birth + sus_death + inf_death + trans;
        
        if (floor(t/10)>t_int)
            t_int = floor(t/10);
            j = j + 1;
            I_vec(j,:) = sum(I,2);
        end
    end

    Ivec_3D(:,:,i) = I_vec;
end

clearvars -except Ivec_3D alpha
save('k_1')