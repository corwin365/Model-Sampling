function Model = load_merra2replay(ObsGrid)

% ObsGrid = load('C:\Data\corwin\sampling_project\tracks\AIRS_3D\track_airs3d_734421_g044.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load MERRA2 replay model data from Laura, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%2. file not found


%model data path
DataDir  = [LocalDataDir,'/corwin/issi/laura/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this input grid
%we're just going to take the closest timestep, as this is for a specific
%gravity wave study and the data are only available every three hours
%so interpolating across time will never be done.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find timestep
Time = mean(ObsGrid.Track.Time(:));

%and pressure grid
MaxPrs = 10.^max(ObsGrid.Track.Prs(:)).*1.2;
MinPrs = 10.^min(ObsGrid.Track.Prs(:))./1.2;


%identify time
[y,M,d,h,m,~] = datevec(Time);

%find nearest three-hour step
ha = 3.*round(h./3);


%and hence the filename
FileName = [DataDir,'/MERRA2_REPLAY_', ...
            sprintf('%04d',y),sprintf('%02d',M),sprintf('%02d',d), ...
            '_',sprintf('%02d',ha),'00z.nc'];

%load file
if ~exist(FileName,'file');
  Model.Error = 2;
  return
end
Data = cjw_readnetCDF(FileName,1);


%make pressure easier to work with
Prs = squeeze(nanmean(Data.P,[1,2])./100);


%drop unwanted pressure levels
Good = find(Prs <= MaxPrs & Prs >= MinPrs);
Data.T  = Data.T(:,:,Good);
Prs     = Prs(Good);


AllData.T    = Data.T;
AllData.Prs  = Prs;
AllData.Time = Time;
AllData.Lat  = Data.lat;
AllData.Lon  = Data.lon;

clear Data Good Prs T FileName y M d h m

%we need 1d lat and lon. 
%So, we need to reinterpolate the data to a regular lat/lon grid
%oversample in both directions, to be safe

%shift the obsgrid lons into 0-360, to make comparison useful
ObsGrid.Track.Lon(ObsGrid.Track.Lon < 0) = ObsGrid.Track.Lon(ObsGrid.Track.Lon < 0)+360;

MinLat = max([min(ObsGrid.Track.Lat(:)),min(AllData.Lat(:))]);
MaxLat = min([max(ObsGrid.Track.Lat(:)),max(AllData.Lat(:))]);
MinLon = max([min(ObsGrid.Track.Lon(:)),min(AllData.Lon(:))]);
MaxLon = min([max(ObsGrid.Track.Lon(:)),max(AllData.Lon(:))]);

Lat = MinLat:0.075:MaxLat;
Lon = MinLon:0.075:MaxLon;

[xi,yi] = meshgrid(Lon,Lat);

%create interpolant object
F = scatteredInterpolant(double(flatten(AllData.Lon)), ...
                         double(flatten(AllData.Lat)), ...
                         double(flatten(AllData.T(:,:,1,1))));

%create storage array
sz = size(AllData.T);
T2 = NaN([size(xi),sz(3)]);

%regrid
for iLevel=1:1:sz(3);
  F.Values = double(flatten(AllData.T(:,:,iLevel)));
  T2(:,:,iLevel) = F(double(xi),double(yi));
end

%and shift the longitudes
Lon(Lon > 180) = Lon(Lon > 180)-360;
[~,idx] = sort(Lon,'asc');
T2 = T2(:,idx,:);
clear idx


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
Model.Time = AllData.Time;
Model.T    = double(permute(T2,[4,2,1,3]));
Model.Prs  = AllData.Prs;



%success!
Model.Error = 0;
return


