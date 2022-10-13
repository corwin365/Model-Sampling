function Model = load_1kmIFS(ObsGrid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load hi-res IFS data for desired period, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%2. file not found

%get core variables - needed for model data path
CoreVars = sampling_core_variables;
CoreVars.OneKmIFS.Path   = [LocalDataDir,'/emily/1km_data'];

%time handling for filenames
BasisTime = datenum(2018,11,1,0,0,0); %all other files are named in MINUTES offset from this
FileStringA = 'ICMSHh3f7+';
FileStringB = '_gg_interpGBA_TenthDeg_NH.nc';
TimeStep    = 1./24./60.*180; %one every three hours
%so, the files are named [CoreVars.OneKmIFS.Path,'/',FileStringA,sprintf('%06d',MINUTES),FileStringB]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this input grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define timestep grid
First = min(ObsGrid.Track.Time(:));
Last  = max(ObsGrid.Track.Time(:));

% %and pressure grid
% MaxPrs = 10.^max(ObsGrid.Track.Prs(:)).*1.2;
% MinPrs = 10.^min(ObsGrid.Track.Prs(:))./1.2;

%generate a list of timesteps. 
FirstStep = (First-BasisTime).*24.*60;
LastStep  = (Last -BasisTime).*24.*60;
t = 0:180:1e8;
[~,idx1] = min(abs(t-FirstStep)); idx1 = idx1.*180;
[~,idx2] = min(abs(t-LastStep )); idx2 = idx2.*180;

Steps = idx1:180:idx2;

clear FirstStep LastStep First Last t idx1 idx2


for Step=Steps

  %identify file
  FileName =  [CoreVars.OneKmIFS.Path,'/',FileStringA,sprintf('%06d',Step),FileStringB];
          
  %load file
  if ~exist(FileName,'file');
    FileName
    disp('file load failed')
    Model.Error = 2;
    return
  end
  Data = cjw_readnetCDF(FileName);  
  
  %replace the outer edge with NaNs, so that data can't be extrapolated beyond the edge
  Data.t(  1,  :,:) = NaN;
  Data.t(end,  :,:) = NaN;
  Data.t(  :,  1,:) = NaN;
  Data.t(  :,end,:) = NaN;
  
  %pull out vars
  T = Data.t;
  Prs = ecmwf_prs_v3(137);
  Lat = Data.latitude;
  Lon = Data.longitude; 
  
  %pull out and reformat data
  if Step == Steps(1);
    AllData.T    = T;
    AllData.Prs  = Prs;
    AllData.Time = Step;
    AllData.Lat  = Lat;
    AllData.Lon  = Lon;
  else
    AllData.T    = cat(4,AllData.T,single(T));
    AllData.Time = cat(1,AllData.Time,Step);
    
  end
  
  clear Data Good Prs T FileName y M d h m 

end; 




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
Model.T    = double(permute(AllData.T,[4,1,2,3]));
Model.Prs  = AllData.Prs;

%shift longitudes into range
Model.Lon(Model.Lon > 180) = Model.Lon(Model.Lon > 180) - 360;

%make lat increase monotonically
Model.Lat = Model.Lat(end:-1:1);
Model.T   = Model.T(:,:,end:-1:1,:);

%success!
Model.Error = 0;
return



