function Model = load_um_issi(ObsGrid)

% % % ObsGrid = load('C:\Data\corwin\sampling_project\tracks\AIRS_3D\track_airs3d_734419_g195.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load UM data for ISSI period from Annelize, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%2. file not found


%get core variables - needed for model data path
CoreVars = sampling_core_variables;
CoreVars.UM_ISSI.Path = [LocalDataDir,'/corwin/issi/annelize/'];


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

  
  
%   2pt5km_L242_AP_TP_20101008T1345_regular_grid.nc
  
  %identify file
  [y,M,d,h,m,~] = datevec(Time);
  FileName = [CoreVars.UM_ISSI.Path, ...
              '2pt5km_L242_AP_TP_', ...
              sprintf('%04d',y),sprintf('%02d',M),sprintf('%02d',d), ...
              'T',sprintf('%02d',h),sprintf('%02d',m),'_regular_grid.nc'];
%   FileName = 'C:\Data\CESM\wrfout_d01_2010-10-08_13-00-00.nc' %temporary override for local testing                      
                   

  %load file
  if ~exist(FileName,'file');
    Model.Error = 2;
    return
  end
  Data = cjw_readnetCDF(FileName,1);  
  
  %pull out vars
  T = Data.STASH_m01s30i004;
  Prs = h2p(Data.STASH_m01s15i102./1000);
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



%success!
Model.Error = 0;
return


