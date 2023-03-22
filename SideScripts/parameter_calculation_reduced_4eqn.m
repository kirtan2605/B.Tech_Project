Tdx = 5e-6;
Tdz = 5e-6;

wo = 2*pi/86400;
wo = 7.236e-5;

phi_ss = deg2rad(0.05);
psi_ss = deg2rad(0.4);

Ix = 800;
Iz = 1000;

xi1 = 0.7;
xi2 = 0.7;

A = (Tdx*psi_ss - Tdx*(psi_ss/phi_ss)) + Tdz;
B = -1*Tdx*phi_ss;
C = wo*psi_ss*phi_ss; 
D = -1*wo*phi_ss*phi_ss;
E = Tdx*(psi_ss/phi_ss);
F = psi_ss;
G = -1*phi_ss;
J = Ix;
K = Iz;
L = xi1;
M = xi2;

syms w x y z

eqn1 = 2*(L*w + M*x)*J == y;
eqn2 = (w*w + x*x + 4*L*M*w*x)*J*K == wo*((A + B*z)/(C + D*z))*(J+K) + K*((E)/(F + G*z)) + ((A + B*z)/(C + D*z))*((A + B*z)/(C + D*z)) + z*y*((A + B*z)/(C + D*z));
eqn3 = 2*w*x*(L*x + M*w)*J*K == z*((A + B*z)/(C + D*z))*(E/(F + G*z)) + wo*((A + B*z)/(C + D*z))*y;
eqn4 = w*w*x*x*J*K == wo*wo*((A + B*z)/(C + D*z))*((A + B*z)/(C + D*z)) + wo*((A + B*z)/(C + D*z))*(E/(F + G*z));
equations =[eqn1, eqn2, eqn3, eqn4];
vars = [w,x,y,z];

S = solve(equations,vars,"Real",true, 'IgnoreAnalyticConstraints',true);

disp(S)