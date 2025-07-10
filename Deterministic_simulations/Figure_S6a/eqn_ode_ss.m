function dzdt = eqn_ode_ss(t, z, b, d, q, al, bt, v, p, f)

% setting up the correct number of ode equations
dxdt = zeros(1, length(v));
dydt = zeros(length(al), length(v));

% reshaping z into a vector (x) containing the susceptibles and a matrix (y)
%containing the infecteds (n_alpha x n_c)
x = z(1:length(v));
y = reshape(z(length(v)+1:end), [length(al), length(v)]);

% total population
N = sum(z);

% force of infection
trans1 = sum(bt'.*sum(v.*y, 2));

dxdt(1,:) = p*(b - q*N)*N - d.*x' - v.*x'*(trans1);

for i = 1:length(al)
    dydt(i,:) = bt(i)*(v.*x')*(sum(v.*y(i,:))) - (d + f*al(i)).*y(i,:);
end

dzdt = [dxdt, reshape(dydt, 1, [])]';
%norm(dzdt,1)
end