clear;

format long

Tdx = 5e-6;
Tdz = 5e-6;

wo = 2*pi/86400;

phi_ss = deg2rad(0.05);
psi_ss = deg2rad(0.4);

Ix = 800;
Iz = 1000;

xi1 = 0.7;
xi2 = xi1;
h = 20;

kx = Tdx/phi_ss - wo*h;

syms wn1 wn2 a_psi kxd

eqn1 = 2*(xi1*wn1 + xi2*wn2)*Ix*Iz == Iz*kxd ;
eqn2 = (wn1^2 + wn2^2 + 4*xi1*xi2*wn1*wn2)*Ix*Iz == (wo*h*(Ix+Iz) + Iz*kx + h^2 + a_psi*h*kxd) ;
eqn3 = (2*wn1*wn2*(xi1*wn2 + xi2*wn1))*Ix*Iz == a_psi*h*kx + wo*h*kxd ;
eqn4 = ((wn1*wn2)^2)*Ix*Iz == (wo*h)^2 + wo*h*kx ;

equations =[eqn1, eqn2, eqn3, eqn4];
vars = [wn1, wn2, a_psi, kxd];

S = solve(equations,vars,"Real",true,'IgnoreAnalyticConstraints',true);

Wn1_combinations = double(S.wn1);
Wn2_combinations = double(S.wn2);
a_psi_combinations = double(S.a_psi);
Kxd_combinations = double(S.kxd);

dimension = length(Wn1_combinations);

%{
disp(Wn1_combinations);
disp(Wn2_combinations);
disp(a_psi_combinations);
disp(Kxd_combinations);
%}

%% eliminating solutions which are not possible

% eliminating solutions with negative frequencies
for i = 1:dimension
    if Wn1_combinations(i) < 0
        Wn1_combinations(i) = 0;
        Wn2_combinations(i) = 0;
        a_psi_combinations(i) = 0;
        Kxd_combinations(i) = 0;
    end

    if Wn2_combinations(i) < 0
        Wn1_combinations(i) = 0;
        Wn2_combinations(i) = 0;
        a_psi_combinations(i) = 0;
        Kxd_combinations(i) = 0;
    end
end

%{
disp("Cleaned Parameters")

disp(Wn1_combinations);
disp(Wn2_combinations);
disp(a_psi_combinations);
disp(Kxd_combinations);
%}

% reducing the cleaned parameters

Wn1_combinations = nonzeros(Wn1_combinations);
Wn2_combinations = nonzeros(Wn2_combinations);
a_psi_combinations = nonzeros(a_psi_combinations);
Kxd_combinations = nonzeros(Kxd_combinations);

%{
disp("reduced Parameters")

disp(Wn1_combinations);
disp(Wn2_combinations);
disp(a_psi_combinations);
disp(Kxd_combinations);
%}

% this will reduce it to 2 possible solutions : Wn1 > Wn2, Wn2 > Wn1
% assuming Wn1 as the orbit rate, hence Wn1 < Wn2

if Wn2_combinations(1) > Wn1_combinations(1)
    Solution_Wn1 = Wn1_combinations(1);
    Solution_Wn2 = Wn2_combinations(1);
    Solution_a_psi = a_psi_combinations(1);
    Solution_Kxd = Kxd_combinations(1);
else
    Solution_Wn1 = Wn1_combinations(2);
    Solution_Wn2 = Wn2_combinations(2);
    Solution_a_psi = a_psi_combinations(2);
    Solution_Kxd = Kxd_combinations(2);
end


% Display Final Parameters
psi_ss_calculated = double( Tdz/(wo*h) + (Solution_a_psi*Tdx*kx)/(wo*h*(wo*h + kx)));
disp(" ");
disp("Final Parameters");
disp(['h : ', num2str(h)]);
disp(['a_psi : ', num2str(Solution_a_psi)]);
disp(['Kx : ', num2str(kx)]);
disp(['Kxd : ', num2str(Solution_Kxd)]);
disp(['Wn1 (orbit rate): ', num2str(Solution_Wn1)]);
disp(['Wn2 (nutation frequency): ', num2str(Solution_Wn2)]);
disp(" ");
disp("Steady State Error Values");
disp(['phi_ss : ', num2str(rad2deg(phi_ss))]);
disp(['psi_ss calculated : ', num2str(rad2deg(psi_ss_calculated))]);
disp(['psi_ss desired : ', num2str(rad2deg(psi_ss))]);


