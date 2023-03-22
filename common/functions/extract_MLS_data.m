function MLSData = extract_MLS_data(MatlabDay,DataDir)

% % % test
% MatlabDay = 733043;
% DataDir = [LocalDataDir,'MLS/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%extract MLS data for a given day - temperatures (but other variables
%available)
%
%Phoebe Noble, pn399@bath.ac.uk (adapted from Corwin Wright's similar HIRDLS
%script)
%09/Feb/2023
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
if MatlabDay < datenum(2004,1,1); MLSData.Error = 1; return; end
if MatlabDay > datenum(2023,1,1); MLSData.Error = 1; return; end

%find data file for this day
[y,m,~] = datevec(MatlabDay);
Day = date2doy(MatlabDay);
y   = sprintf('%04d',y);
m   = sprintf('%02d',m);
Day = sprintf('%03d',Day);

% % Files have four possible formats, so we loop through them
% DayFile = [DataDir,'T/',y,'/MLS-Aura_L2GP-Temperature_v05-00-c01_',y,'d',Day,'.he5'];
% DayFile2 = [DataDir,'T/',y,'/MLS-Aura_L2GP-Temperature_v05-01-c01_',y,'d',Day,'.he5'];
% DayFile3 = [DataDir,'T/',y,'/MLS-Aura_L2GP-Temperature_v05-01-c02_',y,'d',Day,'.he5'];
% DayFile4 = [DataDir,'T/',y,'/MLS-Aura_L2GP-Temperature_v05-00-c02_',y,'d',Day,'.he5'];

%Search for files of the correct format
Search_for_file = dir([DataDir,'T/',y,'/MLS-Aura_L2GP-Temperature_*_',y,'d',Day,'.he5']);
DayFile = [DataDir,'T/',y,'/',Search_for_file.name];

%clear Day y m

%check if the file exists, and load it if so
if exist(DayFile,'file') == 0
  MLSData.Error = 2; 
  return;
end

%extract data
A=get_MLS(DayFile, 'Temperature');

MLSData.T      = A.L2gpValue'; 
MLSData.T(MLSData.T == -999) = NaN; % setting any -999 values to NaN
MLSData.Lat    = A.Latitude;
MLSData.Lon    = A.Longitude;
MLSData.Prs    = A.Pressure;
MLSData.Time   = A.Time; MLSData.Time = datenum(1993,1,1,0,0,MLSData.Time); % converting to MATLAB time
MLSData.Time   = A.Time; MLSData.Time = datenum(1993,1,1,0,0,MLSData.Time); % converting to MATLAB time
MLSData.LineOfSightAngle = A.LineOfSightAngle;

clear DayFile


% Unlike HIRDLS, we can make the assumption with MLS that the lat/lon positions are valid for the entire vertical profile
% Therefore we don't need to account for the drift as the instrument scans up and down
% Now all we need to do is ensure all T,lat,lon,time variables are the same size.

MLSData.Lat = repmat(MLSData.Lat,1,size(MLSData.Prs,1));
MLSData.Lon = repmat(MLSData.Lon,1,size(MLSData.Prs,1));
MLSData.LineOfSightAngle = repmat(MLSData.LineOfSightAngle,1,size(MLSData.Prs,1));
MLSData.Prs = repmat(MLSData.Prs,1,size(MLSData.Lat,1))';


%time finally
%assume zero time over the profile, not quite right but good enough for our purposes
MLSData.Time = repmat(MLSData.Time,1,size(MLSData.Prs,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I think we don't need the rest here! (to the next %%%%% section) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%duplicate out lat and lon to correct size
%%retrieval is on a regular pressure grid, so this is not an approximation
%MLSData.Prs = repmat(MLSData.Prs,1,numel(MLSData.Lat))';

%%lat and lon are more complicated
%MLSData.Lon2 = MLSData.Prs.*NaN;
%MLSData.Lat2 = MLSData.Prs.*NaN;


%TravelAngle = NaN(numel(MLSData.Lon),1);
%for iProf=2:1:numel(TravelAngle);
  
%  %along-track direction (%c/w from E)
%  TravelAngle(iProf) = azimuth(MLSData.Lat(iProf-1),MLSData.Lon(iProf-1), ...
%                               MLSData.Lat(iProf  ),MLSData.Lon(iProf  ));
  
  
%end; clear iProf 


%%geolocation level
%[~,BasisLevel] = min(abs(MLSData.Prs(1,:)-50)); %all MLS prs profiles are identical, so this is fine
%Drift = 0.6; %km along track per vertical level

%Z = p2h(MLSData.Prs(1,:)); % an approximation, but a small one relative to the drift term

%for iZ=1:1:size(MLSData.T,2)
%  Shift = (Z(iZ)./1000-BasisLevel).*Drift.*MLSData.Up; %km drift along-track
%  Shift = distdim(Shift,'km','deg'); %deg of arc
%  [MLSData.Lat2(:,iZ),MLSData.Lon2(:,iZ)] = reckon(MLSData.Lat,MLSData.Lon,Shift,TravelAngle);
%end; clear iZ BasisLevel Drift ScanUp

%MLSData.Lat = MLSData.Lat2;
%MLSData.Lon = MLSData.Lon2;
%MLSData = rmfield(MLSData,{'Lon2','Lat2','Up'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%done!


MLSData.Error = 0; %it worked!

return

end
