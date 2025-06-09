function Model = load_leverhume_urals(DayNumber)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create a simple  idealised field over the Urals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%(it's really pretty basic and just maths)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create the idealised wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Domain (in meters).
x = linspace( -1, 1,      500).*800e3;
y = linspace( -1, 1,      500).*800e3;
z = linspace(  -1, 70e3,  200);
[xi, yi, zi] = ndgrid(x, y, z);

%wave centre location on this grid
x0 = 0; 
y0 = 0; 
z0 = 40e3;

%wave centre location in the real world
CentreLat = 62;
CentreLon = 60;

%wave envelope - stdev of exponential taper
Taper = [3e5,2e5,30e3]; %x,y,z

%wavelengths and amplitude
A  = 10;                        %amplitude,             Kelvin
kx = 2.*pi./1500e3;             %zonal      wavenumber, m^-1
ky = 2.*pi./120e3;              %meridional wavenumber, m^-1
kz = 2.*pi./6e3;                %vertical   wavenumber, m^-1

%generate envelope function
envelope = exp(-(((xi - x0)/Taper(1)).^2 + ((yi - y0)/Taper(2)).^2 + ((zi - z0)/Taper(3)).^2));

%temperature perturbation T
T = 5.*  envelope .* (-cos(kx*xi + ky*yi + kz*zi ));  % Normalized units
clear envelope amplitude kx ky kz  x0 y0 z0 x y A Taper

%create a gridded interpolant for the wave
I = griddedInterpolant(xi,yi,zi,T);
clear xi yi zi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%centre the wave over the Urals, and interpolate to a regional lat/lon grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create a regular lat/lon grid in the chosen region

Lat = (-15:0.05:15) + CentreLat;
Lon = (-15:0.05:15) + CentreLon;
[loni,lati,zo] = ndgrid(Lon,Lat,z);


%convert it to x,y,z coordinates
[r,th] = distance(CentreLat,CentreLon,lati,loni);
r = deg2km(r).*1000;
xo = r.*cosd(90-th); yo = r.*sind(90-th);

%interpolate onto this grid
Tg = I(xo,yo,zo);

clear CentreLat CentreLon loni lati zo xo yo r th I T 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%reformat for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output format:
%struct called Model
%containing fields:
%Lon   - 1d. Runs from -180 to +180
%Lat   - 1d
%Time  - 1d
%Prs   - 1d
%T     - 4d, time x lon x lat x pressure

%first, compute an average pressure scale
%we're in the stratosphere, and have much bigger errors than thuis elsewhere
% stick stuff in a struct
Model.Lon  = Lon;
Model.Lat  = Lat;
Model.Time = DayNumber;
Model.T    = permute(Tg,[4,1,2,3]);
Model.Prs  = h2p(z./1000);

%success!
Model.Error = 0;
return



