clear all;
clc;
% Earth parameters
Req  = 6378.137;		% km, equatorial radius of the earth
muE  = 398600.4415;		% mu_E = GM : gravitational parameter of the earth km^3/s^2, 
				        % taking M_E = 5.976 * 10^24 kg, G = 6.67 * 10^-11
omegaE = 0.729212e-4;	% rad/s, calculated from sidereal day

sidereal_day =  86164.09;    % number of seconds in sidereal day
sma_calculated = ((muE*sidereal_day*sidereal_day)/(4*pi*pi))^(1/3);

%sma_c = 42164.16;   	% km, Brij Agrawal, p. 67; ideal geosynchronous
sma_c = 42164.378;      % km, ChatGPT; ideal geosynchronous
%sma_c = sma_calculated;

%sma = 42166.9;         % km, GEO semimajor axis in the paper; it is greater than the ideal.   
sma = 42171.798;  		% GOES 6
%sma   = sma_c;          % ideal
del_sma = sma - sma_c;

e = 0;                 % The offset and asymmetry in lat and other variables disappears with e = 0
%e = 0.00007;
%e = 0.0005;			% Vinod
%e = 0.000316;			% GOES 6

%i = 0.1 * pi/180;  	    % Vinod
i = 2 * pi/180;		% inclination, rad
%i = 0;                 % ideal
%i = 1.97310 * pi/180;	% GOES 6

Omega = (235.28-360) * pi/180;	% rad, ascending node angle [-pi,pi]
AoP   = 23.05 * pi/180;		    % rad, argument of perigee
disp(Omega)
disp(AoP)

omega0  = sqrt(muE/sma^3);      % rad/s; 
n_synch = sqrt(muE/sma_c^3);
del_n   = - 1.5 * n_synch * del_sma / sma_c;

orbit_period = 2*pi/omega0;
OP_hrs  = floor(orbit_period/3600);
OP_mins = floor((orbit_period - OP_hrs*3600)/60);
OP_sec  = orbit_period - OP_hrs * 3600 - OP_mins * 60;

% step size, seconds, traversing the orbit

dT = 10;                
K  = (23*3600 + 56*60 + 4)/dT;     % plotting for sidereal time period 

% Memory for saving the results
% Nadir angle from a deviated lat-lon caused by i (Brij Agrawal, p. 69) to
% the target P -- the ideal longitude of the geostationary satellite with latitude = 0

%Time    = NaN(K+1,1);
%Del_lat = NaN(K+1,1);
%Del_lon = NaN(K+1,1);

%Time(1,1) = 0;
%Del_lat(1,1) = 0;
%Del_lon(1,1) = 0;

for k = 1:K
    Time(k+1,1) = k * dT / 3600;      % hours
    nu          = omega0 * k * dT;    % orbit angle, rad, circular orbit

    % SSP of the geosynchronous satellite: Brij Agrawal, p. 69
    del_lat =              i * sin(nu);
    del_lon = - 0.25 * i * i * sin(2*nu) + del_n * k * dT + 2*e*sin(nu);  % negligible due to i compared to del_lat

    Del_lat(k+1,1) = del_lat * 180/pi;
    Del_lon(k+1,1) = del_lon * 180/pi;
end

disp('Start points')
disp(Time(1,1))
disp(Del_lat(1,1))
disp(Del_lon(1,1))

disp('End points')
disp(Time(end,1))
disp(Del_lat(end,1))
disp(Del_lon(end,1))


figure(1)
hold on
plot (Del_lon, Del_lat,'Color','b','LineWidth',1.1)
plot(Del_lon(1,1), Del_lat(1,1),'Color',[0.9290, 0.6940, 0.1250],'marker','o','markersize',10,'LineWidth',2)
plot(Del_lon(end,1), Del_lat(end,1),'Color','r','marker','x','markersize',10,'LineWidth',2)
xlabel ('\Delta longitude (deg)')
ylabel ('\Delta latitude (deg)')
legend('Variation','Initial Coordinate','Final Coordinate')
grid on
saveas(gcf,'Lat_Lon.png')
hold off

%figure(2)
%hold on
%plot (Time, Del_lon,'Color','b','LineWidth',1.1),
%plot(Time(1,1), Del_lon(1,1),'Color',[0.9290, 0.6940, 0.1250],'marker','o','markersize',10,'LineWidth',2)
%plot(Time(end,1), Del_lon(end,1),'Color','r','marker','x','markersize',10,'LineWidth',2)
%xlabel ('Time (hours)')
%ylabel ('\Delta longitude (deg)')
%set (gca,'XTick',0:3:24)
%legend('Variation','Initial Coordinate : Ascending Node','Final Coordinate : Ascending Node + \Deltan')
%grid on
%saveas(gcf,'Lon_vs_t.png')
%hold off

%figure(3)
%hold on
%plot(Time, Del_lat, 'Color','b','LineWidth',1.1)
%plot(Time(1,1), Del_lat(1,1),'Color',[0.9290, 0.6940, 0.1250],'marker','o','markersize',10,'LineWidth',2)
%plot(Time(end,1), Del_lat(end,1),'Color','r','marker','x','markersize',10,'LineWidth',2)
%xlabel ('Time (hours)')
%ylabel ('\Delta latitude (deg)')
%set (gca,'XTick',0:3:24)
%legend('Variation','Initial Coordinate : Ascending Node','Final Coordinate : Ascending Node + \Deltan')
%grid on
%saveas(gcf,'Lat_vs_t.png')
%hold off

figure(4)
hold on
yyaxis left
plot (Time, Del_lon,'LineWidth',1.1),
plot(Time(1,1), Del_lon(1,1),'Color',[0.9290, 0.6940, 0.1250],'marker','o','markersize',10,'LineWidth',2)
plot(Time(end,1), Del_lon(end,1),'Color','r','marker','x','markersize',10,'LineWidth',2)
ylabel ('\Delta longitude (deg)')

yyaxis right
plot (Time, Del_lat,'LineWidth',1.1),
plot(Time(1,1), Del_lat(1,1),'Color',[0.9290, 0.6940, 0.1250],'marker','o','markersize',10,'LineWidth',2)
plot(Time(end,1), Del_lat(end,1),'Color','r','marker','x','markersize',10,'LineWidth',2)
ylabel ('\Delta latitude (deg)')

legend('Longitude Variation','Initial Coordinate','Final Coordinate','Latitude Variation')
xlabel ('Time (hours)')
set (gca,'XTick',0:3:24)
grid on
saveas(gcf,'LatLon_vs_t.png')
hold off


