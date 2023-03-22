Ix = 1000;
Iz = Ix;
Iy = 500;   % assumption, Iy < Iz from the satellite configuration

omega_0 = 0.727*(10e-4);
%disp(nutation_f);
%disp(omega_0);
%disp(['     Nutation Frequency / Omega_0 =', num2str(nutation_f/omega_0)]);
% h = k (pitch angular momentum) = k(Iy*omega_0);
% for geostationarity of the satellite oemga_0 is nearly 7.2711e-5

tau = 5;
alpha_d = 10;

a = 4 * omega_0^2 * (Iy - Iz);
b = -1 * omega_0 * (Ix - Iy +Iz);
c = omega_0^2 * (Iy - Iz);
d = 3 * omega_0^2 * (Ix - Iz);


K = 0:1:5000;
n = length(K);
p = zeros([4 n]);

for i = 1:n
h = K(i)*omega_0;
nutation_f = h/(Ix*Iz)^0.5;

K1 = a + omega_0*h;
K2 = b + h;
K3 = c + omega_0*h;

% poles :: roots of characteristic equation
ch_eq = [Ix*Iz 0 (Ix*K3 + Iz*K1 + K2^2) 0 K1*K3];
p(:,i) = roots(ch_eq);
end
hold on
plot(K,abs(p))
plot(K,(K*omega_0)/((Ix*Iz)^0.5))
yline(omega_0)
xlabel('h/omega_0');
ylabel('|poles|');
legend('|pole1|','|pole2|','|pole3|','|pole4|')
hold off