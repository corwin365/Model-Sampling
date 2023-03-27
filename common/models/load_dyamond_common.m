function Model = load_dyamond_common(ObsGrid,WantedModel)

% % ObsGrid = load('C:\Data\corwin\sampling_project\tracks\AIRS_3D\track_airs3d_734419_g195.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load DYAMOND-WINTER data using CDO to trim the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%2. file not found
%3. model does not exist in list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%input data
switch WantedModel
  case 'icon2km';
    Paths.Root = '/work/ka1081/DYAMOND_WINTER/MPIM-DWD-DKRZ/ICON-NWP-2km/DW-ATM/atmos/';
    Paths.Grid = [Paths.Root,'/fx/gn/grid.nc'];
    Paths.T    = [Paths.Root,'/3hr/ta/r1i1p1f1/ml/gn/'];
    Paths.P    = [Paths.Root,'/3hr/pa/r1i1p1f1/ml/gn/'];
  otherwise
    disp(['DYAMOND model ',WantedModel,' not configured, stopping')
    Model.Error = 3;
    return
end

%working scratch space
ScratchPath = '~/scratch/working/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find which files contain the droids^Wtimestamps we're looking for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the time period covered by each file in the directory
Files = wildcardsearch(ModelPath,'*.nc');
for iFile=1:1:numel(Files)

  %call CDO to get the times of the timesteps
  [~,cmdout] = system(['cdo showtimestamp ',Files{iFile}]);

  %split the output into strings, and find all those that are a valid format
  TimeStamps = split(cmdout);
  Times = zeros(size(TimeStamps));ValidCount =0;
  for iTS=1:1:numel(Times);
    if regexp(TimeStamps{iTS},'20[0-9][0-9]-[0-9]+-[0-9]+T') == 1; 
      %valid timestep, parse it
      ValidCount = ValidCount+1;
      Times(ValidCount) = datenum(TimeStamps{iTS},'yyyy-mm-ddTHH:MM:SS');
    end
  end
  Times = Times(1:ValidCount);
  clear cmdout TimeStamps iTS IsT ValidCount

  %if this is the first file, create a time storage array
  if ~exist('TimeStore','var')
    TimeStore = NaN(numel(Files),numel(Times));
  end

  %store
  TimeStore(iFile,1:numel(Times)) = Times;

end; clear iFile


%ok, now find the first and last timestep in each file
FirstTime = min(TimeStore,[],2);
LastTime  = max(TimeStore,[],2);

%and hence the files we need
[delta1,FirstFile] = min(abs(FirstTime-min(ObsGrid.Track.Time(:))));
[delta2,LastFile]  = min(abs( LastTime-max(ObsGrid.Track.Time(:))));
if delta1 > 0; FirstFile = FirstFile-1; end
if delta2 < 0; LastFile  = LastFile+1; end
if FirstFile < 1; FirstFile = 1; end
if LastFile  > size(TimeStore,1); LastFile = size(TimeStore,1); end


Files = Files(FirstFile:1:LastFile);
TimeStore = TimeStore(FirstFile:1:LastFile,:);
clear delta1 delta2 FirstTime LastTime FirstFile LastFile

TimeStore
stop


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, use cdo to cut down the data to take these files and cut them down
%to just the geographic area we need
%this could take a while, but it's much faster than trying to load the whole 
%dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we're going to do this twice - once for pressure and once for temperature - and then merge them at the end
Range = [floor(min(ObsGrid.Track.Lon)),ceil(max(ObsGrid.Track.Lon)),floor(min(ObsGrid.Track.Lat)),ceil(max(ObsGrid.Track.Lat))];
TempFiles = cell(numel(Files),3); %the third dim is for the merged file
for iFile=1:1:numel(Files)

  for iSource=1:2

    switch iSource
      case 1; Source = Paths.P; 
      case 2; Source = Paths.T;
    end

    disp(['Subsetting model file ',num2str(iFile),' of ',num2str(numel(Files)),' using CDO - variable: ',Source])

    %create a temporary filename
    ff = strrep(Files{iFile},ModelPath,ScratchPath);

    %add a unique identifier, just in case multiple processes are working with this file
    %combination of the time to the tenth of a millisecond and a random number between one and a million should hopefully be unique?
    TempFiles{iFile,iSource} = strrep(ff,'.nc',['_',strrep(num2str(datenum(now)),'.',''),'_',num2str(randi(1e6,[1])),'_',Source,'.nc']);
%%%%%%%%%%    TempFiles{iFile,iSource} = strrep(ff,'.nc',['_',Source,'.nc']);
%%%%%%%%%%    if exist(TempFiles{iFile,iSource},'file'); continue; end

    %now, generate a CDO command to subset the data down in space
    Command = 'cdo ';
    Command = [Command,' -sellonlatbox,',num2str(Range(1)),',',num2str(Range(2)),',',num2str(Range(3)),',',num2str(Range(4))];


    %work out which timesteps this *observation* file covers. If none,find the single closest. 
    %Pad by 3 hours each way to make sure we're not *just* missing something.
    InRange = find(TimeStore(iFile,:) <= max(ObsGrid.Track.Time(:)-3/24) ...
                 & TimeStore(iFile,:) >= min(ObsGrid.Track.Time(:)+3/24));
    if numel(InRange) == 0; [~,InRange] = min(abs(TimeStore(iFile,:)-mean(ObsGrid.Track.Time(:)))); end

    %if the number of timesteps is less than the total in the file, subset the file in time

    if numel(InRange) < size(TimeStore,2)
      Command = [Command, ' -seltimestep'];
      for iTime=1:1:numel(InRange); Command = [Command,',',num2str(InRange(iTime))]; end
    end

    %finalise the command with the names of the in and out files
    Command = [Command,' ',strrep(Files{iFile},'ta',Source),' ',TempFiles{iFile,iSource}];

    %and execute it
    status = system(Command);
    if status ~=0;
      %error!
      Model.Error = 2;
      return
    end

  end; clear iSource Source

  %generate a CDO command to merge the pressure data with the temperature data
  TempFiles{iFile,3} = strrep(TempFiles{iFile,2},'ta','merged'); 
  Command = ['cdo merge ',TempFiles{iFile,1},' ',TempFiles{iFile,2},' ',TempFiles{iFile,3}];
%%%%%%%%%%if exist(TempFiles{iFile,3},'file'); continue; end
  %and execute it
  status = system(Command);
  if status ~=0;
    %error!
    Model.Error = 2;
    return
  end
  disp(['File ',num2str(iFile),' of ',num2str(numel(Files)),' merged'])


end; clear iFile
clear Range Files ff InRange Command status PrsFile TimeStore


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% great! Let's load the files and get what we need out of them
%memory can be *really* tight here, so try to use as little as possible
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iFile=1:1:size(TempFiles,1)

  %load the file
  if iFile ==1;
    %if it's the first file, do some prep  
    Model.Lon  = ncread(TempFiles{iFile,3},'longitude');
    Model.Lat  = ncread(TempFiles{iFile,3},'latitude');
    Model.Prs  = ncread(TempFiles{iFile,3},'pfull');
    Model.Time = datenum(1970,1,1,ncread(TempFiles{iFile,3},'time'),0,0); 
    Model.T    = ncread(TempFiles{iFile,3},'ta');
  else
    %otherwise, append
    stop
    Model.Time = [Model.Time,datenum(1970,1,1,Data.time,0,0)];
    Model.T    = cat(1,Model.T,  permute(   Data.ta,[4,1,2,3]));
    Model.Prs  = cat(1,Model.Prs,permute(Data.pfull,[4,1,2,3]));
  end

 
end; clear iFile


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% delete working files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iFile=1:1:numel(TempFiles)
  delete(TempFiles{iFile})
  disp(['Tidiying up: ',TempFiles{iFile},' deleted'])
end; clear iFile



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the data onto a common pressure scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%what dim is the one we want to interpolate?
Dim = 3;

%permute desired dimension to front
sz = size(Model.T);
DimOrder = unique([Dim,1:1:numel(sz)],'stable');

%get the mean pressure scale
PrsScale = squeeze(mean(Model.Prs,DimOrder(2:end),'omitnan'));


%reshape to make all other dimensions lines
T = reshape(permute(  Model.T,DimOrder),[sz(Dim),prod(sz(DimOrder(2:end)))]);
P = reshape(permute(Model.Prs,DimOrder),[sz(Dim),prod(sz(DimOrder(2:end)))]);

%interpolate
T2 = NaN(numel(PrsScale),size(T,2));
for iProfile=1:1:size(T,2); T2(:,iProfile) = interp1(P(:,iProfile),T(:,iProfile),PrsScale); end
  
%reshape back
Ti = reshape(T2,[numel(PrsScale),sz(DimOrder(2:end))]);

%and permute back
NewOrder = 1:1:numel(sz);
NewOrder = [NewOrder(NewOrder < Dim)+1,1,NewOrder(NewOrder > Dim)];
Model.T   = permute(Ti,NewOrder);
Model.Prs = PrsScale./100; %also converting Pa to hPa here

clear Dim PrsScale DimOrder T P Ti Pi



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% finally, permute the results to the desired output order
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output format:
%struct called Model
%containing fields:
%Lon   - 1d. Runs from -180 to +180
%Lat   - 1d
%Time  - 1d
%Prs   - 1d
%T     - 4d, time x lon x lat x pressure

%switch around T dims
Model.T = permute(Model.T,[4,1,2,3]);

%shift longitudes into range
Model.Lon(Model.Lon > 180) = Model.Lon(Model.Lon > 180) - 360;

%success!
Model.Error = 0;
return


