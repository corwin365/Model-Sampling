function Model = load_1kmIFS(ObsGrid,BlobScale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load hi-res IFS data for desired period, and put into common analysis format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%2. file not found

% model data path
ModelPath   = '/perm/paip/AIRS_SAMPLING_CODE/data/Jan20onwards/';

%time handling for filenames
BasisTime = datenum(2018,11,1,0,0,0); %all other files are named in MINUTES offset from this
FileStringA = 't_130_h3f7-ndjf_';
FileStringB = '.ml-1to80.0.075x0.075GBA.global.nc';
TimeStep    = 1./24./60.*180; %one every three hours
%so, the files are named [ModelPath,'/',FileStringA,sprintf('%06d',MINUTES),FileStringB]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load the data for this input grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%put obsgrid onto 0-360, same as the model, so we can subset correctly UNLIKE LAST TIME.
%AGAIN, YOU'VE DONE THIS LOADS OF TIMES WRONG YOU MORON
%relative geometry doesn't matter as it's just a list of points which we treat independently.
ObsGrid.Track.Lon(ObsGrid.Track.Lon <0) = ObsGrid.Track.Lon(ObsGrid.Track.Lon < 0) + 360;

%define timestep grid
First = min(ObsGrid.Track.Time(:));
Last  = max(ObsGrid.Track.Time(:));

%define range covered by the dataset in terms of measurement centre locations
Range = NaN(5,4);
Range(1,:) = [floor(min(ObsGrid.Track.Lon)),ceil(max(ObsGrid.Track.Lon)), ...
              floor(min(ObsGrid.Track.Lat)),ceil(max(ObsGrid.Track.Lat))];

%find the maximum extra distance a blob could add to the edge of the region (with 25% padding)
MaxDistance = max([ObsGrid.Weight.X;ObsGrid.Weight.Y]).*BlobScale.*1.25;

%find where this puts the range, by taking a line diagonally from each corner
Range(2,[4,2]) = reckon2(Range(1,4),Range(1,2),km2deg2(MaxDistance),45);
Range(3,[3,2]) = reckon2(Range(1,3),Range(1,2),km2deg2(MaxDistance),135);
Range(4,[3,1]) = reckon2(Range(1,3),Range(1,1),km2deg2(MaxDistance),225);
Range(5,[4,2]) = reckon2(Range(1,4),Range(1,1),km2deg2(MaxDistance),315);

r(1) = min(Range(:,[1,2]),[],'all');
r(2) = max(Range(:,[1,2]),[],'all');
r(3) = min(Range(:,[3,4]),[],'all');
r(4) = max(Range(:,[3,4]),[],'all');
Range = r; clear r BlobScale MaxDistance

if Range(1) < -180; Range(1) = -180; end
if Range(2) >  180; Range(2) =  180; end
if Range(3) <  -90; Range(3) =  -90; end
if Range(4) >   90; Range(4) =   90; end


% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % %put obsgrid onto 0-360, same as the model, so we can subset correctly UNLIKE LAST TIME.
% % % % % % %relative geometry doesn't matter as it's just a list of points which we treat independently.
% % % % % % ObsGrid.Track.Lon(ObsGrid.Track.Lon <0) = ObsGrid.Track.Lon(ObsGrid.Track.Lon < 0) + 360;
% % % % % % 
% % % % % % %define range covered by the dataset, with a bit of padding (20% of the range on each side)
% % % % % % Range = [floor(min(ObsGrid.Track.Lon(:))),ceil(max(ObsGrid.Track.Lon(:))), ...
% % % % % %          floor(min(ObsGrid.Track.Lat(:))),ceil(max(ObsGrid.Track.Lat(:)))];
% % % % % % Range = Range ...
% % % % % %       + [-1,1,0,0].*0.2.*(max(Range([1,2])) - min(Range([1,2]))) ...
% % % % % %       + [0,0,-1,1].*0.2.*(max(Range([3,4])) - min(Range([3,4])));
% % % % % % if Range(1) <    0; Range(1) =    0; end   %THESE TWO LINES SHOULD BE FROM 0-360
% % % % % % if Range(2) >  360; Range(2) =  360; end   %NOT -180 TO +180, BECAUSE I AM AN IDIOT.
% % % % % % if Range(3) <  -90; Range(3) =  -90; end
% % % % % % if Range(4) >   90; Range(4) =   90; end


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
   FileName =  [ModelPath,'/',FileStringA,sprintf('%06d',Step),FileStringB];
%   FileName =  [ModelPath,'/',FileStringA,FileStringB];
  FileName         
  %load file
  if ~exist(FileName,'file');
    FileName
    disp('file load failed')
    Model.Error = 2;
    return
  end
  disp('before read')                                 
  Data = rCDF(FileName);  
  disp('after read')   

  %for testing
  % % Data = rCDF('C:\Data\ERA5\2010\era5_2010d002.nc');
  % % Data = load('C:\Data\ERA5\2010\era5_2010d002_360.mat');  


  %pull out vars
  T = Data.t;
  Prs = ecmwf_prs_v3(137);
  Prs = Prs(1:1:80);
  Lat = Data.latitude;
  Lon = Data.longitude; 

  %cut out the geographic region we want, to save memory
  InLatRange = find(Lat >= Range(3) & Lat <= Range(4));
  InLonRange = find(Lon >= Range(1) & Lon <= Range(2));
  Lat = Lat(InLatRange); Lon = Lon(InLonRange);
  T = T(InLonRange,:,:,:); T = T(:,InLatRange,:,:);
  clear InLatRange InLonRange


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

clear FileStringA FileStringB ModelPath ObsGrid Range Step Steps BasisStep TimeStep


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
[~,idx] = sort(Model.Lon,'asc');
Model.Lon = Model.Lon(idx);
Model.T = Model.T(:,idx,:,:);


%make lat increase monotonically
Model.Lat = Model.Lat(end:-1:1);
Model.T   = Model.T(:,:,end:-1:1,:);


%success!
Model.Error = 0;
return



