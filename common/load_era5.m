function Model = load_era5(DayNumber,MaxPrs)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ERA5 data for a particular day, and put into common analysis format
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

%path for this day
[y,m,d] = datevec(DayNumber);
FileName1 = [CoreVars.ERA5.Path,sprintf('%04d',y),'/',sprintf('%02d',m),'/',sprintf('%02d',d),'/era5_',sprintf('%04d',y),sprintf('%02d',m),sprintf('%02d',d)];
%path for next day
[y2,m2,d2] = datevec(DayNumber+1);
FileName2 = [CoreVars.ERA5.Path,sprintf('%04d',y),'/',sprintf('%02d',m),'/',sprintf('%02d',d),'/era5_',sprintf('%04d',y),sprintf('%02d',m),sprintf('%02d',d)];


%work out pressure axis here - memory demands are huge, so we want to drop
%elements as soon as possible

PrsScale = ecmwf_prs_v2(11.06059,137)'; %close enough for the s'sphere
PrsScale(1) = 0.01; %because it is, but my routine gives a NaN

GoodPrs = find(PrsScale < MaxPrs); %hPa - omits most of t'sphere


PrsScale = PrsScale(GoodPrs);
try
  for iHour=1:1:24;
    
    %load file
    if iHour == 24 %first of next day
      File = [FileName2,sprintf('%02d',0),'.nc'];      
    else %normal hour
      File = [FileName1,sprintf('%02d',iHour-1),'.nc'];
    end
    
    Data.longitude = full(ncArray(File,'longitude'));
    Data.latitude  = full(ncArray(File,'latitude'));
    Data.T         = full(ncArray(File,'t')); 
    Data.T = Data.T(:,:,GoodPrs);

    if iHour == 1
      %create arrays to store all the time points in
      sz = [size(Data.T),24]; sz = sz([4,1,2,3]);
      AllData.T = NaN(sz);
      AllData.Time = NaN(24,1);
    end
    
    %store
    AllData.Time(iHour)  = datenum(y,m,d,iHour,0,0);
    AllData.T(iHour,:,:,:) = Data.T;

  end
catch
  Model.Error = 2;
  return
end

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
Model.Time = AllData.Time;
Model.T    = AllData.T;
Model.Prs  = PrsScale;

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


