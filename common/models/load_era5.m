function Model = load_era5(DayNumber,MaxPrs,MinPrs)

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
ModelPath   = [LocalDataDir,'/ERA5/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path for this day and the next day
[y,~,~] = datevec(DayNumber);   dn = date2doy(DayNumber);
FileName1 = [ModelPath,sprintf('%04d',y),'/era5_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.nc'];
[y,~,~] = datevec(DayNumber+1); dn = date2doy(DayNumber+1);
FileName2 = [ModelPath,sprintf('%04d',y),'/era5_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.nc'];


%work out pressure axis and select wanted region
PrsScale = ecmwf_prs_v3(137,11.06059)'; %close enough for the stratosphere. 
PrsScale(1) = 0.01; %because it is, but my routine gives a NaN
idx = inrange(PrsScale,[MinPrs,MaxPrs]);

%load the two days, cut them down to just the desired region, and merge them
Day1   = rCDF(FileName1);
try     Day1   = rmfield(Day1,{'u','v','lnsp','level'});%,'MetaData'});
catch;  Day1   = rmfield(Day1,{'u','v',       'level'});%,'MetaData'});
end
Day1.t = Day1.t(:,:,idx,:);

Day2   = rCDF(FileName2);
try     Day2   = rmfield(Day2,{'u','v','lnsp','level'});%,'MetaData'});
catch;  Day2   = rmfield(Day2,{'u','v',       'level'});%,'MetaData'});
end
Day2.t = Day2.t(:,:,idx,:);

Store      = Day1; clear Day1
Store.t    = cat(4,Store.t,Day2.t); 
Store.time = cat(1,Store.time,Day2.time);
clear Day2

%convert time units
Store.time = datenum(1900,1,1,Store.time,0,0);

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
Model.Lon  = Store.longitude;
Model.Lat  = Store.latitude;
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


%latitude is also descending - we want ascending. Oh, ECMWF...
Model.Lat = flip(Model.Lat,1);
Model.T   = flip(Model.T,3);

%success!
Model.Error = 0;
return


