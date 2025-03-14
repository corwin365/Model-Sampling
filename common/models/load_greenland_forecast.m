function Model = load_greenland_forecast(ObsGrid,ForecastTime,MaxPrs)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ECMEF forecast data, and put into common analysis format
%note that the logic assumes the points are sufficiently close in time to 
%use a single forecast. If points are spread by more than ~6-12 hours, this
%assumption may break dow. If so, subset the day externally and do each bit
%here individually
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%1. date not in valid range
%2. file not found

%model data path
ModelPath  = [LocalDataDir,'/corwin/ecmwf_forecast'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%work out pressure axis before extraction, so we can store only what we need
PrsScale = ecmwf_prs_v3(137,11.06059)'; %close enough for the s'sphere
PrsScale(1) = 0.01; %because it is, but my routine gives a NaN
% GoodPrs = find(PrsScale < MaxPrs); %hPa - omits most of t'sphere
% clear MaxPrs

%load the model file for the specified forecast
%load the file
[~,~,~,h,~,~] = datevec(ForecastTime);
FileName = [ModelPath,'/greenland_',yyyyDDD(ForecastTime),'h',sprintf('%02d',h),'_fullres.nc'];
if ~exist(FileName,'file');  error(['Error: forecast file ',FileName,' does not exist. Have you specified the right time?']); end
Data = rCDF(FileName);
Data.time = datenum(1900,1,1,Data.time,0,0);
Data.Prs  = PrsScale(Data.level);
GoodPrs = find(Data.Prs < MaxPrs);
Data.Prs = Data.Prs(GoodPrs); Data.t = Data.t(:,GoodPrs,:,:);


%now, pull out the timesteps covering the observations
TimeRange = [min(ObsGrid.Track.Time,[],'all','omitnan'),max(ObsGrid.Track.Time,[],'all','omitnan')];

idx = inrange(Data.time,TimeRange);
if     numel(idx) == 0; idx = closest(mean(Data.time,'all'),TimeRange)+[-1:1:1];
elseif numel(idx) == 1; idx = [idx-1:1:idx+1];
else                    idx = [min(idx)-1:1:max(idx)+1];
end

Data.t = Data.t(idx,:,:,:);
Data.time = Data.time(idx);


%put NaNs around the edge, so we don't extrapolate off it
Data.t(:,:,  1,:) = NaN;
Data.t(:,:,end,:) = NaN;
Data.t(:,:,:,  1) = NaN;
Data.t(:,:,:,end) = NaN;

    
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
Model.Lon  = Data.longitude;
Model.Lat  = Data.latitude;
Model.Time = Data.time;
Model.T    = permute(Data.t,[1,4,3,2]);
Model.Prs  = Data.Prs;

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


