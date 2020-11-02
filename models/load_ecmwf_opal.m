function Model = load_ecmwf_opal(DayNumber)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load ECMWF operational analysis data for a particular day, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%1. date not in valid range
%2. file not found

%get core variables - needed for model data path
CoreVars = sampling_core_variables;
CoreVars.EcOpAl.Path   = [LocalDataDir,'/ECMWF/'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this day
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%ERAi file format:
%monthly files, with p based on a standard grid
Settings.DataDir    = CoreVars.EcOpAl.Path;

%load the appropriate files
[y,m,d] = datevec(DayNumber);




%load data
for iChunk=1:1:5;
  
  ChunkStart = (iChunk-1).*6;
  if ChunkStart ~= 24; 
    [y,m,d] = datevec(DayNumber);
  else
    [y,m,d] = datevec(DayNumber+1); 
    ChunkStart = 0;
  end
  
  FileString = [Settings.DataDir,'/',sprintf('%04d',y), ...
              '/ecmwf91lev_',sprintf('%04d',y),sprintf('%02d',m),sprintf('%02d',d),sprintf('%02d',ChunkStart),'.mat'];
  if ~exist(FileString); continue; end            
  Data = load(FileString);
  Data = Data.EcmwfData;
  
  if ~exist('T');
    T = Data.T;
    U = Data.U;
    V = Data.V;
    Lat = Data.latitude;
    Lon = Data.longitude;
    LNSP = Data.LNSP;
    Time = ChunkStart./24;
  else
    T = cat(4,T,Data.T);
    U = cat(4,U,Data.U);
    V = cat(4,V,Data.V);
    Time(end+1) = ChunkStart./24;
  end
  
end

if ~exist('T'); Model.Error = 2; return; end

%work out prs axis
Pressure = ecmwf_prs_v2(LNSP,size(T,1));
Pressure = squeeze(nanmean(nanmean(Pressure,1),2));
Pressure(1) = 0.01;

%and order lons right
Lon(Lon > 180) = Lon(Lon > 180)-360;
[~,Order] = sort(Lon);
Lon = Lon(Order);
T = T(:,:,Order,:);
U = U(:,:,Order,:);
V = V(:,:,Order,:);

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


%then, stick stuff in a struct
Model.Lon  = Lon;
Model.Lat  = Lat;
Model.Time = DayNumber+Time;
Model.T    = permute(T,[4,3,2,1]);
Model.U    = permute(U,[4,3,2,1]);
Model.V    = permute(V,[4,3,2,1]);
Model.Prs  = Pressure;

%success!
Model.Error = 0;
return


