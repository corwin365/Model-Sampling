function Model = load_ecmwf_issi(ObsGrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load UM data for ISSI period from Annelize, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%2. file not found


%path to data
DataDir = [LocalDataDir,'/corwin/issi/inna/'];

%time handling for file names
%this does not correspond to normal ECMWF time handling due to the way Inna 
%set up the runs, so we need to be careful.
BasisTime      = datenum(2010,10,09,12,0,0); %this is file "ifs_t001.nc"
TimeStep       = 15./24./60;                 %this is the time step for each additional file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this input grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define timestep grid
Step = 1./24./60.*15;
First = min(ObsGrid.Track.Time(:));
Last  = max(ObsGrid.Track.Time(:));

%and pressure grid
MaxPrs = 10.^max(ObsGrid.Track.Prs(:)).*1.2;
MinPrs = 10.^min(ObsGrid.Track.Prs(:))./1.2;


%round 'first' and 'last' time to the nearest 15-min step, rounding away from the sample
First = floor(First./Step).*Step;
Last  = ceil(  Last./Step).*Step;

for Time=First:Step:Last;

  
  %identify file
  ThisStep = round((Time-BasisTime)./TimeStep);
  FileName = [DataDir,'ifs_t',sprintf('%03d',ThisStep),'.nc'];
          

  %load file
  if ~exist(FileName,'file');
    Model.Error = 2;
    return
  end
  Data = cjw_readnetCDF(FileName);  
  
  %pull out vars
  T = Data.t;
  Prs = ecmwf_prs_v2([],137);
  Lat = Data.latitude;
  Lon = Data.longitude;
  
  
  %pull out and reformat data
  if Time == First;
    AllData.T    = T;
    AllData.Prs  = Prs;
    AllData.Time = Time;
    AllData.Lat  = Lat;
    AllData.Lon  = Lon;
  else
    AllData.T    = cat(4,AllData.T,single(T));
    AllData.Time = cat(1,AllData.Time,Time);
    
  end
  
  clear Data Good Prs T FileName y M d h m 

end; clear Time First Last Step


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
Model.T    = double(permute(AllData.T,[4,2,1,3]));
Model.Prs  = AllData.Prs;

%shift longitudes into range
Model.Lon(Model.Lon > 180) = Model.Lon(Model.Lon > 180) - 360;

%success!
Model.Error = 0;
return


