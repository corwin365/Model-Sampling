function Model = load_era5_levante(ObsGrid,MaxPrs,MinPrs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ERA5 data for a particular day, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%1. date not in valid range
%2. file not found

%path to model data
ModelPath   = '/scratch/b/b382226/era5/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MeanTime = nanmean(ObsGrid.Track.Time(:));
DayNumber = floor(MeanTime);

%path for this hour and the next hour
[y,~,~,h,~,~] = datevec(MeanTime); dn = floor(date2doy(DayNumber));
FileName1 = [ModelPath,'/era5_',sprintf('%04d',y),'d',sprintf('%03d',dn),'h',sprintf('%02d',h),'.nc'];
[y,~,~,h,~,~] = datevec(MeanTime+1/24);   dn = floor(date2doy(MeanTime+1/24));
FileName2 = [ModelPath,'/era5_',sprintf('%04d',y),'d',sprintf('%03d',dn),'h',sprintf('%02d',h),'.nc'];


%work out pressure axis and select wanted region
PrsScale = ecmwf_prs_v3(137,11.06059)'; %close enough for the stratosphere. 
PrsScale(1) = 0.01; %because it is, but my routine gives a NaN
PrsScale = PrsScale(1:80); %we trimmed the files, starting at the top
idx = inrange(PrsScale,[MinPrs,MaxPrs]);

%load the two days, cut them down to just the desired region, and merge them
Day1   = rCDF(FileName1);
Day1.t = Day1.t(:,:,idx,:);
Day1.time = DayNumber+Day1.time./24;

Day2   = rCDF(FileName2);
Day2.t = Day2.t(:,:,idx,:);
Day2.time = DayNumber+Day2.time./24+1;

Store      = Day1; clear Day1
Store.t    = cat(4,Store.t,Day2.t); 
Store.time = cat(1,Store.time,Day2.time);
clear Day2


%tidy up
clear idx MaxPrs MinPrs dn y FileName1 FileName2 ModelPath DayNumber


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
Model.Lon  = Store.lon;
Model.Lat  = Store.lat;
Model.Time = Store.time;
Model.T    = permute(Store.t,[4,1,2,3]);
Model.Prs  = PrsScale'; 

%longitude is 0-360, so we need to rejiggle into -180 to 180
Lon = Model.Lon;
Lon(Lon > 180) = Lon(Lon > 180)-360;
[~,idx] = sort(Lon,'ascend');
Lon = Lon(idx);
Model.Lon = Lon;
Model.T   = Model.T(:,idx,:,:);

%duplicate the end points, in case profiles go over the dateline
Model.T = cat(2,Model.T(:,end,:,:),Model.T,Model.T(:,1,:,:));
Model.Lon = [Model.Lon(1)-mean(diff(Model.Lon));Model.Lon;Model.Lon(end)+mean(diff(Model.Lon))];


%latitude is also descending - we want ascending. Oh, ECMWF...
Model.Lat = flip(Model.Lat,1);
Model.T   = flip(Model.T,3);

%success!
Model.Error = 0;
return


