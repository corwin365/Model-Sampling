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
%ModelPath  = [LocalDataDir,'/corwin/ecmwf_forecast/2023_greenland2/'];
% ModelPath = '/data2/peter/crossover/forecast/small/';
ModelPath = [LocalDataDir,'/peter/crossover/forecast/small/'];
% ModelPath = '/scratch/b/b382226/ecmwf_forecast';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%work out pressure axis before extraction, so we can store only what we need
PrsScale = ecmwf_prs_v3(137,11.06059)'; %close enough for the s'sphere
PrsScale(1) = 0.01; %because it is, but my routine gives a NaN
% GoodPrs = find(PrsScale < MaxPrs); %hPa - omits most of t'sphere
% clear MaxPrs

%disp(datestr(ObsGrid.Track.Time(1)))
%stop

%load the model file for the specified forecast
%load the file
%[~,~,~,h,~,~] = datevec(ForecastTime);
%FileName = [ModelPath,'/globalgreenland_',yyyyDDD(ForecastTime),'h',sprintf('%02d',h),'_fullres.nc'];

forecastTime = datetime(ForecastTime, 'ConvertFrom', 'datenum');
folderName = fullfile(ModelPath, ['greenland_',num2str(year(forecastTime)),'d',sprintf('%03d',day(forecastTime,'dayofyear')),'h',sprintf('%02d',hour(forecastTime))]);

%This is the forecast from which I want it initialised
%Need to work out which time step is wanted
folderName;

%now, pull out the timesteps covering the observations
TimeRange = [min(ObsGrid.Track.Time,[],'all','omitnan'),max(ObsGrid.Track.Time,[],'all','omitnan')];
if isa(TimeRange, 'double')
    TimeRange = datetime(TimeRange,'ConvertFrom','epochtime','Epoch','2000-01-01');
end
%stop
TimeRange
forecastTime
steps = round(hours(TimeRange - forecastTime))

allTimeSteps = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,93,96,99,102,105,108,111,114,117,120,123,126,129,132,135,138,141,144,150,156,162,168,174,180,186,192,198,204,210,216,222,228,234,240];

[~,idxStart] = min(abs(steps(1) - allTimeSteps));
[~,idxStop] = min(abs(steps(2) - allTimeSteps));

a = dir(fullfile(folderName, 'greenland*'));

for forcSteps = idxStart:idxStop
    a(forcSteps).name
    Data = rCDF(fullfile(folderName, a(forcSteps).name));
    DataN.time(forcSteps-idxStart+1) = datenum(1900,1,1,Data.time,0,0);
    DataN.Prs  = PrsScale(Data.level);
    GoodPrs = find(DataN.Prs < MaxPrs);
    DataN.Prs = DataN.Prs(GoodPrs);
    DataN.t(forcSteps-idxStart+1,:,:,:) = Data.t(GoodPrs,:,:);
    DataN.latitude = Data.latitude;
    DataN.longitude = Data.longitude;

end

DataN.t = permute(DataN.t, [1, 4, 3, 2]);

Data = DataN;


%if ~exist(FileName,'file');  error(['Error: forecast file ',FileName,' does not exist. Have you specified the right time?']); end
%Data = rCDF(FileName);
%Data.time = datenum(1900,1,1,Data.time,0,0);
%Data.Prs  = PrsScale(Data.level);
%GoodPrs = find(Data.Prs < MaxPrs);
%Data.Prs = Data.Prs(GoodPrs); Data.t = Data.t(:,GoodPrs,:,:);



%idx = inrange(Data.time,TimeRange);
%if     numel(idx) == 0; idx = closest(mean(TimeRange),Data.time)+[-1:1:1];
%elseif numel(idx) == 1; idx = [idx-1:1:idx+1];
%else                    idx = [min(idx)-1:1:max(idx)+1];
%end
%idx
%Data.t = Data.t(idx,:,:,:);
%Data.time = Data.time(idx);



    
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
%we're in the stratosphere, and have much bigger errors than this elsewhere

% stick stuff in a struct
Model.Lon  = Data.longitude;
Model.Lat  = Data.latitude;
Model.Time = Data.time;
%Model.T    = permute(Data.t,[1,4,3,2]);
Model.T = Data.t;
size(Model.T)
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

%Get the step that was used
Model.step = allTimeSteps(idxStart);


%success!
Model.Error = 0;
return


