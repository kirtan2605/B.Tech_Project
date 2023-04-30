% Parameters
Req  = 6378.137;		% km, equatorial radius of the earth   
sma = 42171.798;  		% km, semi-major axis of orbit

% Earth disc and central angle to the horizon: Sec. 5.2 Wertz and Larson

h_GEO  = sma - Req;		% GEO altitude

% Earth central angle and disc radius variation with height

Lamda0 = NaN(1,795);
Rho    = NaN(1,795);

Height = 300:50:40000;		% height km
J      = 795;			% no. of points in Height

SMA    = Req + Height;

for j = 1:J
    Lamda0(1,j) = acos( Req / SMA(1,j)) * 180/pi;
    Rho(1,j)    = 90 - Lamda0(1,j);
end

%%%%%%%%%%%%%%%%%%%% Plots %%%%%%%%%%%%%%%%%%%%

figure (1)
plot (Height, Lamda0, Height, Rho, [h_GEO,  h_GEO], [0, Lamda0(1,J)])
xlabel ('Satellite altitude (km)')
ylabel ('Earth central Angle \lambda_0, disc radius \rho (deg) ')
legend ('\lambda_0','\rho','GEO')
grid on

figure(2)
semilogx(Height, Lamda0, Height, Rho, [h_GEO,  h_GEO], [0, Lamda0(1,J)])
xlabel ('Satellite altitude (km)')
ylabel ('Earth central Angle \lambda_0, disc radius \rho (deg) ')
legend ('\lambda_0','\rho','GEO')
grid on