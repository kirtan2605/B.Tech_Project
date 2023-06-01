Ix = 1475.8425;
Iy = 886.242;   % assumption, Iy < Iz from the satellite configuration
Iz = 1532.96;

h = 20;
omega_0 = (2*pi)/86164;

a = 4 * omega_0^2 * (Iy - Iz);
b = -1 * omega_0 * (Ix - Iy +Iz);
c = omega_0^2 * (Iy - Iz);
d = 3 * omega_0^2 * (Ix - Iz);

K = 0:-0.0001:-5;
n = length(K);
p = zeros([4 n]);

for i = 1:n
h = K(i);

A = a + omega_0*h;
B = b + h;
C = c + omega_0*h;

% poles :: roots of characteristic equation
ch_eq = [Ix*Iz 0 (Ix*C + Iz*A + B*B) 0 A*C];
p(:,i) = roots(ch_eq);
end

hold on
plot(-K,log(abs(p(1,:))/omega_0),'LineWidth',1.1)
plot(-K,log(abs(p(3,:))/omega_0),'LineWidth',1.1)
plot(-K,log((K/((Ix*Iz)^0.5))/omega_0),'--','LineWidth',1.1)
plot(-K,log((sqrt(((K-(omega_0*(Ix-Iy))).*(K-(omega_0*(Iz-Iy))))./(Ix*Iz))./omega_0)),'--','LineWidth',1.1)
yline(log(omega_0/omega_0),'--','LineWidth',1.1)
ylim([-1 5])
xlim([0,4])
xlabel('h (Nms)');
ylabel('log(w/w0)');
title('Splitting of Poles due to Momentum Bias')
legend('log(w1/w0) : Nutation Frequency','log(w3/w0) : Precession Frequency','Nutation Frequency Binomial Approximation','Nutation Frequency Approximation by Iwens','Orbit Frequency','Location','northwest')
saveas(gcf,'PoleSplitting.png')
hold off
