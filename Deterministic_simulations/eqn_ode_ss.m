function dzdt = eqn_ode_ss(t, z, b, d, q, al, bt, c, p)

dxdt = zeros(1, length(c));
dydt = zeros(length(al), length(c));

x = z(1:length(c));
y = reshape(z(length(c)+1:end), [length(al), length(c)]);
N = sum(z);

trans1 = sum(bt'.*sum(c.*y, 2));

dxdt(1,:) = p*(b - q*N)*N - d.*x' - c.*x'*(trans1);

for i = 1:length(al)
    dydt(i,:) = bt(i)*(c.*x')*(sum(c.*y(i,:))) - (d + al(i)).*y(i,:);
end

dzdt = [dxdt, reshape(dydt, 1, [])]';
%norm(dzdt,1)
end