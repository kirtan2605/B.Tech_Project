Ix = 1000;
Iz = Ix ;
h = 150;

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
rlocus ( H_roll ) , grid

figure(2)
G = ( tau * s - 1) ;
rlocus ( G * H_roll ) , grid

figure(3)
G = ( tau * s + 1) ;
rlocus ( G * H_roll ) , grid
