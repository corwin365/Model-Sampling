function [Data,PlanetaryWaves] = remove_pws(Data,Settings)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for each height level individually, remove PWs 
%assumes they remian constant over the day
%
%Corwin Wright, c.wright@bath.ac.uk
%10/NOV/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find planetary waves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define a grid to interpolate the data onto inside the loop
LonGrid = -180:20:180; %only used for computation.
LatGrid =  -90:Settings.LatBin:90;

%create an array to store the PW data
PlanetaryWaves = NaN(Settings.NPWs+1, ... %pws and mean
                    numel(LatGrid),   ... %latitude bins
                    size(Data.T,2),   ... %altitude levels
                    numel(LonGrid));      %longitude bins

%meshgrid the lat,lon for easier use
[xi,yi] = meshgrid(LonGrid,LatGrid);



%loop over levels
for iLevel=1:1:size(Data.T,2);
  
  %get data at this level
  T   = Data.T(  :,iLevel);
  Lat = Data.Lat(:,iLevel);
  Lon = Data.Lon(:,iLevel);
  
  %interpolate onto a regular lat-lon grid
  %this will give nonsense outside the geographic range of the instriment
  %but since we don't need background here anyway that's fine
  Good = find(~isnan(T+Lat+Lon)); %NaNs are bad
  F = scatteredInterpolant(Lon(Good),Lat(Good),T(Good));
  if numel(Good) == 0; continue; end
  Tg = F(xi,yi);

  %find mean at each latitude
  PlanetaryWaves(1,:,iLevel,:) = repmat(nanmean(Tg,2),1,numel(LonGrid));
  Tg = Tg-squeeze(PlanetaryWaves(1,:,iLevel,:)); %residual T
  
  %now find all the modes at each latitude
  for iMode=1:1:Settings.NPWs
    for iLat=1:1:numel(LatGrid);
      
      WaveFreq = 1./(360./iMode); %per degrees
      
      [~,Yest] = sinefit(squeeze(Tg(iLat,:)),LonGrid,WaveFreq,0,0,0);
      PlanetaryWaves(1+iMode,iLat,iLevel,:) = Yest; %PW at this lat
      Tg(iLat,:) = Tg(iLat,:)-Yest;
      
      
    end
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%remove planetary waves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data.Tp = Data.T;

for iLevel=1:1:size(Data.T,2);
  
  for iLat=1:1:numel(LatGrid);
    InLatBin = find(Data.Lat(:,iLevel) >= LatGrid(iLat) ...
                  & Data.Lat(:,iLevel) <  LatGrid(iLat)+Settings.LatBin);
    if numel(InLatBin) == 0; continue; end
    
    for iMode = 1:1:Settings.NPWs+1; %i.e. including mean
      
      ToRemove = interp1(LonGrid,squeeze(PlanetaryWaves(iMode,iLat,iLevel,:)),Data.Lon(InLatBin,iLevel));
      
      Data.Tp(InLatBin,iLevel) = Data.Tp(InLatBin,iLevel)-ToRemove;
      

    end
  end
end

%done!

return
