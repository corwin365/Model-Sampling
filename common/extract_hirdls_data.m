function HirdlsData = extract_hirdls_data(MatlabDay,DataDir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract COSMIC data for a given day
%
%Corwin Wright, corwin.wright@trinity.oxon.org
%08/Jan/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%errors:
%0. no error - successful
%1. no file found
%2. no data on this day (not used)
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find the file and load it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if the day is before or after the data record starts, give up!
if MatlabDay < datenum(2005,1,1); HirdlsData.Error = 1; return; end
if MatlabDay > datenum(2008,4,1); HirdlsData.Error = 1; return; end

%find data file for this day
[y,m,~] = datevec(MatlabDay);
Day = date2doy(MatlabDay);
y   = sprintf('%04d',y);
m   = sprintf('%02d',m);
Day = sprintf('%03d',Day);

DayFile = [DataDir,'/',y,'/',m,'/HIRDLS-Aura_L2_v07-00-20-c01_',y,'d',Day,'.he5'];
clear Day y m

%check if the file exists, and load it if so
if exist(DayFile,'file') == 0;
  HirdlsData.Error = 1; 
  return;
end

%extract data
HirdlsData.T      = get_HIRDLS(DayFile,'Temperature')'; HirdlsData.T(HirdlsData.T == -999) = NaN;
HirdlsData.Lat    = get_HIRDLS(DayFile,'Latitude'   );
HirdlsData.Lon    = get_HIRDLS(DayFile,'Longitude'  );
HirdlsData.Prs    = get_HIRDLS(DayFile,'Pressure'   );
HirdlsData.Up     = get_HIRDLS(DayFile,'ScanUpFlag' ); HirdlsData.Up(HirdlsData.Up == 0) = -1;
HirdlsData.Time   = get_HIRDLS(DayFile,'Time'       ); HirdlsData.Time = datenum(1993,1,1,0,0,HirdlsData.Time);
clear DayFile

%duplicate out lat and lon to correct size
%retrieval is on a regular pressure grid, so this is not an approximation
HirdlsData.Prs = repmat(HirdlsData.Prs,1,numel(HirdlsData.Lat))';


%lat and lon are more complicated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HirdlsData.Lon2 = HirdlsData.Prs.*NaN;
HirdlsData.Lat2 = HirdlsData.Prs.*NaN;


TravelAngle = NaN(numel(HirdlsData.Lon),1);
for iProf=2:1:numel(TravelAngle);
  
  %along-track direction (%c/w from E)
  TravelAngle(iProf) = azimuth(HirdlsData.Lat(iProf-1),HirdlsData.Lon(iProf-1), ...
                               HirdlsData.Lat(iProf  ),HirdlsData.Lon(iProf  ));
  
  
end; clear iProf 


%geolocation level
[~,BasisLevel] = min(abs(HirdlsData.Prs(1,:)-50)); %all HIRDLS prs profiles are identical, so this is fine
Drift = 0.6; %km along track per vertical level

Z = p2h(HirdlsData.Prs(1,:)); % an approximation, but a small one relative to the drift term

for iZ=1:1:size(HirdlsData.T,2);
  Shift = (Z(iZ)./1000-BasisLevel).*Drift.*HirdlsData.Up; %km drift along-track
  Shift = distdim(Shift,'km','deg'); %deg of arc
  [HirdlsData.Lat2(:,iZ),HirdlsData.Lon2(:,iZ)] = reckon(HirdlsData.Lat,HirdlsData.Lon,Shift,TravelAngle);
end; clear iZ BasisLevel Drift ScanUp

HirdlsData.Lat = HirdlsData.Lat2;
HirdlsData.Lon = HirdlsData.Lon2;
HirdlsData = rmfield(HirdlsData,{'Lon2','Lat2','Up'});

%time finally
%assume zero time over the profile, not quite right but good enough for our purposes
HirdlsData.Time = repmat(HirdlsData.Time,1,size(HirdlsData.Prs,2));

%done!

HirdlsData.Error = 0; %it worked!

return
end
