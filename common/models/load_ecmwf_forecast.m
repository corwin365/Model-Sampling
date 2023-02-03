function Model = load_ecmwf_forecast(DayNumber,ForecastHours,MaxPrs)

% % % clearvars
% % % DayNumber = datenum(2018,1,214);
% % % datestr(DayNumber)
% % % ForecastHours = 0;
% % % MaxPrs = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ERA5 data for a particular day, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%1. date not in valid range
%2. file not found

%model data path
ModelPath  = [LocalDataDir,'/ECMWF_fc'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%work out pressure axis before extraction, so we can store only what we need
PrsScale = ecmwf_prs_v2(11.06059,137)'; %close enough for the s'sphere
PrsScale(1) = 0.01; %because it is, but my routine gives a NaN
GoodPrs = find(PrsScale < MaxPrs); %hPa - omits most of t'sphere


%REMEMBER: WE ARE WORKING WITH FORECAST DATA
%THAT MEANS THAT WE DON'T WANT THE RUN WITH THE DAY SPECIFIED IN THE TITLE
%WE WANT THE FILE GENERATED X HOURS BEFORE IT 

%information needed:
%1. Year
%2. DayNumber
%3. time of day (fraction)
%for up to four files (since the data are twelve-hourly, this will ensure we always cover every point)
%the files must be at 00, 12 or 24 hrs UTC

Meta = NaN(4,3); 
for iMeta=1:1:4
  
  ExtraShift = rem(ForecastHours,12);
  
  TimeShift = (iMeta-2).*.5; %half a day is twelve hours
  TimeShift = TimeShift -  (ForecastHours./24) + ExtraShift./24;
  [y,~,~] = datevec(DayNumber+TimeShift);
  dn = floor(date2doy(DayNumber+TimeShift));
  tod = rem(date2doy(DayNumber+TimeShift),1);
  Meta(iMeta,:) = [y,dn,tod];
  
end


%based on this, identify the four files we're actually talking about
Files = {};
Times = NaN(4,1); %these are the times THAT ARE BEING FORECAST
for iFile=1:1:4
  
  %identify filename
  Files{iFile} = ['forecast_',sprintf('%04d',Meta(iFile,1)),'d',sprintf('%02d',Meta(iFile,2)),...
                  'h',sprintf('%02d',24.*Meta(iFile,3)),'f',sprintf('%03d',ForecastHours),'.nc'];
  
  %and time that is being forecast
  Times(iFile) = datenum(Meta(iFile,1),1,Meta(iFile,2))+Meta(iFile,3) + ForecastHours./24;
                
end
clear Meta dn ExtraShift ForecastHours iFile iMeta TimeShift tod y DayNumber


%load the files, and store the data we want
for iFile=1:1:4
  
  
  %load file
  File = [CoreVars.EC_FC.Path,'/',Files{iFile}];
  
  if ~exist(File,'file')
    %no file found - quit
    Model.Error = 2;
    return
  end

  %extract data of interest 
  Data.longitude = full(ncArray(File,'longitude'));
  Data.latitude  = full(ncArray(File,'latitude'));
  Data.T         = full(ncArray(File,'t')); 
  Data.T = Data.T(:,:,GoodPrs);

  %store said data
  if iFile == 1
    AllData.T = Data.T;
    AllData.longitude = Data.longitude;
    AllData.latitude  = Data.latitude;
    AllData.Prs = PrsScale(GoodPrs);
  else
    AllData.T = cat(4,AllData.T,Data.T);
  end

  %and loop
  clear Data
  
end

%store times THAT ARE BEING FORECAST
AllData.Time = Times;
clear Times



    
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
Model.Lon  = AllData.longitude;
Model.Lat  = AllData.latitude;
Model.Time = AllData.Time;
Model.T    = permute(AllData.T,[4,1,2,3]);
Model.Prs  = PrsScale;

%longitude is 0-360, so we need to rejiggle into -180 to 180
Lon = Model.Lon;
Lon(Lon > 180) = Lon(Lon > 180)-360;
[~,idx] = sort(Lon,'ascend');
Lon = Lon(idx);
Model.Lon = Lon;
Model.T   = Model.T(:,idx,:,:);


%latitude is also descending - we want ascending. 
Model.Lat = flip(Model.Lat,1);
Model.T   = flip(Model.T,3);



%success!
Model.Error = 0;
return


