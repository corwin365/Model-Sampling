function Model = load_ecmwf_sonja(DayNumber)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ECMWF forecast data supplied by Sonja Gisinger, 
%and put into common analysis format
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

%work out pressure axis
PrsScale = ecmwf_prs_v2(11.06059,137)'; %close enough for the s'sphere
PrsScale(1) = 0.01; %because it is, but my routine gives a NaN


try
    
  [y,~,~] = datevec(DayNumber);
  dn = date2doy(DayNumber);
  File = [CoreVars.SONJA.Path,'/',sprintf('%04d',y),'/sonja_',sprintf('%04d',y),'d',sprintf('%03d',dn),'.nc'];
  clear y dn
  
  Data.longitude = full(ncArray(File,'longitude'));
  Data.latitude  = full(ncArray(File,'latitude'));
  Data.T         = full(ncArray(File,'t'));
  Data.Time      = datenum(1900,1,1,double(full(ncArray(File,'time'))),0,0);
  levs           = full(ncArray(File,'level'));
  
  Data.Prs = PrsScale(levs);
  clear PrsScale

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


% stick stuff in a struct
Model.Lon  = Data.longitude;
Model.Lat  = Data.latitude;
Model.Time = Data.Time;
Model.T    = Data.T;
Model.Prs  = Data.Prs;

%reorder Model.T
Model.T = permute(Model.T,[4,1,2,3]);

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


