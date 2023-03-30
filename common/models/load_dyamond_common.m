function Model = load_dyamond_common(ObsGrid,WantedModel,FixedPFlag)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load DYAMOND-WINTER data, using CDO to trim it for memory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%2. file not found
%3. model does not exist in list

%if FixedPFlag exists and is nonzero, then a fixed pressure scale specified 
%below will be used, otherwise the true pressure from the model will be 
%computed (which basically doubles the runtime)


FixedPFlag = 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%input data
switch WantedModel
  case 'icon5km';
    PressureVar = 'pfull';
    Paths.Root  = '/work/ka1081/DYAMOND_WINTER/MPIM-DWD-DKRZ/ICON-SAP-5km/DW-ATM/atmos/';
    Paths.Grid  = [Paths.Root,'/3hr/ta/dpp0014/ml/gn/grid.nc'];
    Paths.T     = [Paths.Root,'/3hr/ta/dpp0014/ml/gn/'];
    Paths.P     = [Paths.Root,'/3hr/',PressureVar,'/dpp0014/ml/gn/'];
    GridSize    = [0.03,0.03]; %degrees lon,lat
    FixedP      = [1.57, 1.97, 2.46, 3.05, 3.78, 4.66, 5.74, 7.04, 8.60, 10.5, 12.7, 15.3, 18.3, 21.7, 25.6, 30, 34.9, 40.3, 46.2, 52.6, 59.5, 66.8, 74.5, 82.7, 91.1, 99.8, 109, 118, 127, 136, 146, 156, 166, 177, 189, 201, 214, 228, 242, 257, 273, 289, 307, 325, 344, 364, 385, 407, 429, 453, 478, 504, 530, 557, 583, 610, 637, 663, 689, 715, 740, 765, 789, 813, 835, 857, 878, 898, 916, 934, 950, 964, 977, 987, 996, 1000, 1010];
  otherwise
    disp(['DYAMOND model ',WantedModel,' not configured, stopping'])
    Model.Error = 3;
    return
end

%working scratch space
ScratchPath = '~/scratch/working/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find which files contain the droids^Wtimestamps we're looking for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the time period covered by each file in the directory
Files = wildcardsearch(Paths.T,'*.nc');

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, use cdo to cut down the data to take these files and cut them down
%to just the geographic area we need
%this could take a while, but it's much faster than trying to load the whole 
%dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we're going to do this twice - once for pressure and once for temperature - and then merge them at the end
Range = [floor(min(ObsGrid.Track.Lon)),ceil(max(ObsGrid.Track.Lon)),floor(min(ObsGrid.Track.Lat)),ceil(max(ObsGrid.Track.Lat))];


if FixedPFlag == 1;  InnerLoop = 1; TempFiles = cell(numel(Files),1); 
else                 InnerLoop = 2; TempFiles = cell(numel(Files),3); %the first dim is for the merged file
end

for iFile=1:1:numel(Files)
  for iSource=1:1:InnerLoop %only do second loop if we want to extract true pressure
    switch iSource
      case 2; Source = PressureVar;
      case 1; Source = 'ta';
    end

    disp(['Subsetting model file ',num2str(iFile),' of ',num2str(numel(Files)),' using CDO - variable: ',Source])


    %add a unique identifier, just in case multiple processes are working with this file
    %combination of the time to the tenth of a millisecond and a random number between 
    %one and a million should hopefully be unique?
    % Identifier = [strrep(num2str(datenum(now)),'.',''),'_',num2str(randi(1e6,[1]))]);
    Identifier = '';
   
    %for testing, non-unique filenames:
    TempFiles{iFile,iSource+1} = [ScratchPath,Source,'_',sprintf('%06d',iFile),Identifier,'.nc'];
    if exist(TempFiles{iFile,iSource+1},'file'); continue; end

    %now, generate a CDO command to subset the data down in space
    Command = 'cdo -v -P 16'; %verbose output, 16 cores permitted
    Command = [Command,' -sellonlatbox,',num2str(Range(1)),',',num2str(Range(2)),',',num2str(Range(3)),',',num2str(Range(4))];

    %interpolate onto a lat/lon grid that is representative of the data
    if exist('GridSize','var')
      nx = round(360./GridSize(1));
      ny = round(180./GridSize(2));
      Command = [Command,' -remapnn,r',num2str(nx),'x',num2str(ny),' '];
      clear nx ny
    end

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

    %call the grid
    if isfield(Paths,'Grid')
      Command = [Command,' -setgrid,',Paths.Grid,' '];
    end


    %finalise the command with the names of the in and out files
    Command = [Command,' ',strrep(Files{iFile},'ta',Source),' ',TempFiles{iFile,iSource+1}]

    
    %and execute it
    status = system(Command);
    if status ~=0;
      %error!
      Model.Error = 2;
      return
    end

  end; clear iSource Source

  %generate a CDO command to merge the pressure data with the temperature data
  if FixedPFlag ~= 1; 
    TempFiles{iFile,1} = strrep(TempFiles{iFile,2},'ta','merged'); 

    Command = ['cdo merge ',TempFiles{iFile,2},' ',TempFiles{iFile,3},' ',TempFiles{iFile,1}];
    if exist(TempFiles{iFile,3},'file'); continue; end
    
    %and execute it
    status = system(Command);
    if status ~=0;
      %error!
      Model.Error = 2;
      return
    end
    disp(['File ',num2str(iFile),' of ',num2str(numel(Files)),' merged'])
  end

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
    Model.Lon  = ncread(TempFiles{iFile,1},'lon');
    Model.Lat  = ncread(TempFiles{iFile,1},'lat');
    if FixedPFlag~=1;Model.Prs  = ncread(TempFiles{iFile,1},PressureVar); end
    Model.Time = datenum(1970,1,1,ncread(TempFiles{iFile,1},'time'),0,0); 
    Model.T    = ncread(TempFiles{iFile,1},'ta');
    if FixedPFlag~=1; Model.Prs  = ncread(TempFiles{iFile,1},PressureVar); end
  else
    %otherwise, append
    Model.Time = [Model.Time,datenum(1970,1,1,Data.time,0,0)];
    Model.T    = cat(1,Model.T,  permute(           Data.ta,[4,1,2,3]));
    if FixedPFlag~=1; Model.Prs  = cat(1,Model.Prs,permute(Data.(PressureVar),[4,1,2,3])); end
  end


end; clear iFile


% % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % %% delete working files
% % % % % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % 
% % % % % % % % 
% % % % % % % % for iFile=1:1:numel(TempFiles)
% % % % % % % %   delete(TempFiles{iFile})
% % % % % % % %   disp(['Tidiying up: ',TempFiles{iFile},' deleted'])
% % % % % % % % end; clear iFile
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the data onto a common pressure scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FixedPFlag ~= 1

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

end
  

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

