function Model = load_um_forecast(DayNumber,ObsGrid)

% % % clearvars
% % % DayNumber = datenum(2018,1,214);
% % % datestr(DayNumber)
% % % ForecastHours = 0;
% % % MaxPrs = 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load UM forecast data from annelize, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%1. date not in valid range
%2. file not found

%get core variables - needed for model data path
CoreVars = sampling_core_variables;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%only one day of data, so fail if we didn't ask for that
if DayNumber ~= datenum(2018,08,12); 
  Model.Error = 1;
  return
end

%identify the hour before and the hour after the airs granule
[~,~,~,h1,m1,~] = datevec(min(ObsGrid.Track.Time)); h1 = floor(h1+m1./60);
[~,~,~,h2,m2,~] = datevec(max(ObsGrid.Track.Time)); h2 = ceil( h2+m2./60);
clear m1 m2


%load all the timesteps needed, and glue together
clear AllData
for h=h1:1:h2
  Data = cjw_readnetCDF([CoreVars.UM_FC.Path,'/annelize_',num2str(h),'.nc'],1);
  
  
  if ~exist('AllData','var')
    AllData.T = single(Data.STASH_m01s16i004);
    AllData.Prs  = h2p(Data.ML_THETA_zsea_theta./1000);
    AllData.Time = datenum(2018,8,12,h,0,0);
    AllData.Lat  = Data.latitude_t;
    AllData.Lon  = Data.longitude_t;
  else
    AllData.T    = cat(4,AllData.T,single(Data.STASH_m01s16i004));
    AllData.Time = cat(1,AllData.Time,datenum(2018,8,12,h,0,0));
  end
  
  clear Data
end


  
%longitude is 0-360, so we need to rejiggle into -180 to 180
Lon = AllData.Lon;
Lon(Lon > 180) = Lon(Lon > 180)-360;
[~,idx] = sort(Lon,'ascend');
Lon = Lon(idx);
AllData.Lon = Lon;
AllData.T   = AllData.T(idx,:,:,:);


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
Model.Lon  = AllData.Lon;
Model.Lat  = AllData.Lat;
Model.Time = AllData.Time;
Model.T    = double(permute(AllData.T,[4,1,2,3]));
Model.Prs  = AllData.Prs;



%success!
Model.Error = 0;
return


