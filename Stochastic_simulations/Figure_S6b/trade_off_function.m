function f = trade_off_function(alpha, alpha_min, alpha_max)
a = -0.25;
beta_min = 2e-05;
beta_max = 1e-03;
f = beta_min + ((beta_max - beta_min)*(1 - (alpha_max - alpha)/(alpha_max - alpha_min)))./(1 + a*(alpha_max - alpha)/(alpha_max - alpha_min));
end
