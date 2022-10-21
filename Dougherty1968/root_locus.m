Ix = 1000;
Iz = Ix;
h = 150;

nutation_f = h/(Ix*Iz)^0.5;
omega_0 = 0.727*(10e-4);
%disp(nutation_f);
%disp(omega_0);
disp(['     Nutation Frequency / Omega_0 =', num2str(nutation_f/omega_0)]);
% h = k (pitch angular momentum) = k(Iy*omega_0);
% for geostationarity of the satellite oemga_0 is nearly 7.2711e-5

tau = 25;
alpha_d = 10;

nr3 = tau*Ix;
nr2 = (tau*h*tand(alpha_d) + Ix);
nr1 = (tau*omega_0*h + h*tand(alpha_d));
nr0 = omega_0*h;
nr = [nr3 nr2 nr1 nr0];

dr4 = Ix*Iz;
dr3 = 0;
dr2 = Ix*Iz*(omega_0^2) + h*h;
dr1 = 0;
dr0 = h^2 * omega_0^2;
dr = [dr4 dr3 dr2 dr1 dr0];


% realistic K
% 1 + K.H(s) where H(s) = a(s)/b(s) --> b(s) + K.a(s) = 0
K = 0:1:1000;
n = length(K);
p = zeros([4 n]);
for i = 1:n
    ch_eq = [dr4 (dr3+K(i)*nr3) (dr2+K(i)*nr2) (dr1+K(i)*nr1) (dr0+K(i)*nr0)];
    p(:,i) = roots(ch_eq);
end
plot(real(p), imag(p), '.'),grid
xlabel('real');
ylabel('imaginery');


%rlocus(nr,dr),grid