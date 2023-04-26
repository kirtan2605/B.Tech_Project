Ix = 1475.8425 ;
Iz = 1532.9 ;
h = 20;

omega_0 = 0.72722*(10^-4) ;

% nutation frequency for large h
nutation_f_approx = h /(( Ix * Iz ) ^0.5) ;

K = 15;
tau = 25;
alpha_d = 10; % alpha in degrees

s = tf('s');

% H is used to denote the transfer function of the plant .
H_roll = cosd ( alpha_d ) *( Ix * s ^2 + h * tand ( alpha_d ) * s + omega_0 *h ) /( Ix * Iz *( s ^2 + nutation_f_approx ^2) *( s ^2 + omega_0 ^2 ) ) ;

figure(1)
rlocus ( H_roll ), grid
xlim([-0.002 0.0001])
ylim([-0.0005 0.0005])
%set(gca,'XTickLabel',[0.0002 0.0001 0]);
%set(gca,'YTickLabel',[1.5 1.45 1.4 1.35 1.3 1.25]);
zeta = [0.99 0.975  0.95 0.9 0.85 0.75 0.5 0.25 0];
wn = [0.0025 0.0015 0.001 0.0005 0];
sgrid(zeta,wn);
title('Orbit Rate Poles Root Locus')
saveas(gcf,'RL_01_OrbitPoles.png')

figure(2)
rlocus ( H_roll ), grid
xlim([-0.0005 0.0025])
ylim([-0.04 0.04])
%set(gca,'XTickLabel',[0.0002 0.0001 0]);
%set(gca,'YTickLabel',[1.5 1.45 1.4 1.35 1.3 1.25]);
zeta = [0.05 0.005 0];
wn = [0.1 0.05 0.03 0];
sgrid(zeta,wn);
title('Nutation Poles Root Locus')
saveas(gcf,'RL_01_NutationPoles.png')

figure(3)
G = ( tau * s - 1) ;
rlocus ( G * H_roll ) , grid
xlim([-0.001 0.0005])
ylim([-0.0002 0.0002])
%set(gca,'XTickLabel',[0.0002 0.0001 0]);
%set(gca,'YTickLabel',[1.5 1.45 1.4 1.35 1.3 1.25]);
zeta = [0.99 0.975  0.95 0.9 0.5 0];
wn = [0.0025 0.0015 0.001 0.0005 0];
sgrid(zeta,wn);
title('Orbit Rate Poles Root Locus')
saveas(gcf,'RL_02_OrbitPoles.png')

figure(4)
G = ( tau * s - 1) ;
rlocus ( G * H_roll ) , grid
xlim([-0.025 0.00001])
ylim([-0.02 0.02])
%set(gca,'XTickLabel',[0.0002 0.0001 0]);
%set(gca,'YTickLabel',[1.5 1.45 1.4 1.35 1.3 1.25]);
zeta = [0.95 0.5 0.1 0];
wn = [0.05 0.03 0.02 0.01 0.005 0];
sgrid(zeta,wn);
title('Nutation Poles Root Locus')
saveas(gcf,'RL_02_NutationPoles.png')

figure(5)
G = ( tau * s + 1) ;
rlocus ( G * H_roll ) , grid
xlim([-0.002 0.00005])
ylim([-0.002 0.002])
%set(gca,'XTickLabel',[0.0002 0.0001 0]);
%set(gca,'YTickLabel',[1.5 1.45 1.4 1.35 1.3 1.25]);
zeta = [0.95 0.75 0.5 0.1 0];
wn = [0.0025 0.0015 0.001 0.0005 0];
sgrid(zeta,wn);
title('Orbit Rate Poles Root Locus')
saveas(gcf,'RL_03_OrbitPoles.png')

figure(6)
G = ( tau * s + 1) ;
rlocus ( G * H_roll ) , grid
xlim([-0.2 0.0001])
ylim([-0.1 0.1])
%set(gca,'XTickLabel',[0.0002 0.0001 0]);
%set(gca,'YTickLabel',[1.5 1.45 1.4 1.35 1.3 1.25]);
zeta = [0.95 0.5 0.1 0];
wn = [0.25 0.2 0.15 0.1 0.05 0];
sgrid(zeta,wn);
title('Nutation Poles Root Locus')
saveas(gcf,'RL_03_NutationPoles.png')
