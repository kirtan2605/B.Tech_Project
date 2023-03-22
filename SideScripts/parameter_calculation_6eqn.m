clear;

Tdx = 5e-6;
Tdz = 5e-6;

%wo = 2*pi/86400;
wo = 7.236e-5;

phi_ss = deg2rad(0.05);
psi_ss = deg2rad(0.4);

Ix = 800;
Iz = 1000;

xi1 = 0.7;
xi2 = 0.7;

J = Ix;
K = Iz;
L = xi1;
M = xi2;

syms w x y z a h

eqn1 = psi_ss == (Tdz/(wo*h)) - (a*Tdx*z)/(wo*h*(wo*h + z));
eqn2 = phi_ss*(wo*h + z) == Tdx;
eqn3 = 2*J*(L*w + M*x) == y;
eqn4 = (w*w + x*x + 4*L*M*x*w)*J*K == wo*h*(J+K) + K*z + h*h + a*h*y;
eqn5 = 2*w*x*(L*x + M*w)*J*K == a*h*z + wo*h*y;
eqn6 = w*w*x*x*J*K == wo*wo*h*h*wo*h*z;

equations =[eqn1, eqn2, eqn3, eqn4, eqn5, eqn6];
vars = [w,x,y,z,a,h];

S = solve(equations,vars,"Real",true, 'IgnoreAnalyticConstraints',true);

disp(S)