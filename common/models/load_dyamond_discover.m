function Model = load_dyamond_discover(ObsGrid,WantedModel,FixedPFlag)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load DYAMOND-WINTER data, using CDO to trim it for memory.This is the version
%of the GEOS data on NASA's DISCOVER system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible errors:
%0: success
%2. file not found
%3. model does not exist in list

%if FixedPFlag exists and is nonzero, then a fixed pressure scale specified 
%below will be used, otherwise the true pressure from the model will be 
%computed (which basically doubles the runtime)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input data
PressureVar    = 'pfull';
TemperatureVar = 'ta';
    Paths.Root     = '/work/ka1081/DYAMOND_WINTER/MPIM-DWD-DKRZ/ICON-SAP-5km/DW-ATM/atmos/';
    Paths.Grid     = [Paths.Root,'/3hr/ta/dpp0014/ml/gn/grid.nc'];
    Paths.T        = [Paths.Root,'/3hr/ta/dpp0014/ml/gn/'];
GridSize       = [0.01,0.01]; %degrees lon,lat
FixedP         = [1008.3652       1006.1125       1003.6658       1001.0101       998.13165       995.01276       991.63477        987.98083       984.02991       979.76385       975.15991       970.19617       964.85162       959.10510        952.93085       946.30682       939.21198       931.62286       923.52069       914.88068       905.68677        895.92084       885.56885       874.61884       863.06305       850.89728       838.11896       824.73474        810.75177       796.18555       781.05774       765.39398       749.22852       732.59888       715.55164        698.13708       680.41010       662.42670       644.25299       625.95697       607.60797       589.27515        571.02869       552.93750       535.06641       517.47650       500.26202       483.49414       467.20297        451.38226       436.02444       421.11563       406.64743       392.61063       378.99368       365.78625        352.97913       340.56210       328.52609       316.86224       305.56061       294.61130       284.00595        273.73602       263.79272       254.16740       244.85287       235.83958       227.11940       218.68471        210.52814       202.64204       195.01849       187.64998       180.53036       173.65205       167.00713        160.59064       154.39568       148.41571       142.64375       137.07347       131.70009       126.51762        121.51954       116.70026       112.05444       107.57777       103.26315       99.107353       95.104362        91.248917       87.538086       83.963913       80.518509       77.193420       73.981903       70.876968        67.872124       64.961655       62.141167       59.405766       56.750366       54.173119       51.673130        49.250134       46.903931       44.634087       42.440006       40.321148       38.276867       36.306076        34.407806       32.580643       30.823042       29.133589       27.511152       25.954039       24.460752        23.030113       21.660688       20.350769       19.098953       17.904091       16.764410       15.678376        14.644983       13.662410       12.729040       11.843728       11.004971       10.211002       9.4605923        8.7524252       8.0846062       7.4559150       6.8651080       6.3103070       5.7903557       5.3039913        4.8493214       4.4253492       4.0306673       3.6635652       3.3233337       3.0081916       2.7167423        2.4482894       2.2010057       1.9739232       1.7659978       1.5756230       1.4023839       1.2446179        1.1013640      0.97193444      0.85479808      0.74967796      0.65521342      0.57071102      0.49539959       0.42822772      0.36900148      0.31642285      0.27038151      0.22983353      0.19446401      0.16369103       0.13698566      0.11397021     0.094084047     0.077148452     0.062610246     0.050387368     0.039958704      0.031233277     0.023860855     0.017500000     0.012000000    0.0072500003    0.0032572495] ;


%working scratch space
ScratchPath = '~/data/Corwin/';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find which files contain the droids^Wtimestamps we're looking for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get the time period covered by each file in the directory
Files = wildcardsearch(Paths.T,'*.nc');

for iFile=1:1:numel(Files)
  % % % % % 
  % % % % % %call CDO to get the times of the timesteps
  % % % % % [~,cmdout] = system(['cdo showtimestamp ',Files{iFile}]);
  % % % % % 
  % % % % % %split the output into strings, and find all those that are a valid format
  % % % % % TimeStamps = split(cmdout);
  % % % % % Times = zeros(size(TimeStamps));ValidCount =0;
  % % % % % for iTS=1:1:numel(Times);
  % % % % %   if regexp(TimeStamps{iTS},'20[0-9][0-9]-[0-9]+-[0-9]+T') == 1; 
  % % % % %     %valid timestep, parse it
  % % % % %     ValidCount = ValidCount+1;
  % % % % %     Times(ValidCount) = datenum(TimeStamps{iTS},'yyyy-mm-ddTHH:MM:SS');
  % % % % %   end
  % % % % % end
  % % % % % Times = Times(1:ValidCount);
  % % % % % clear cmdout TimeStamps iTS IsT ValidCount

  f = Files{iFile};
  a = strfind(f,'Mv.');
  b = strfind(f,'z.nc4');
  TimeStamp = f(a+3:b);
  yy = TimeStamp(1:4);
  mm = TimeStamp(5:6);
  dd = TimeStamp(7:8);
  hh = TimeStamp(10:11);
  Times = datenum(str2num(yy),str2num(mm),str2num(dd),str2num(hh),00,00);

  %if this is the first file, create a time storage array
  if ~exist('TimeStore','var')
    TimeStore = NaN(numel(Files),numel(Times));
  end

  %store
  TimeStore(iFile,1:numel(Times)) = Times;

end; clear iFile


%ok, now find the first and last timestep in each file
TimeStore(TimeStore == 0) = NaN;
FirstTime = nanmin(TimeStore,[],2);
LastTime  = nanmax(TimeStore,[],2);

%and hence the files we need
%an edge case is arising in COSMIC regional subsets where these are reserved - if so, flip them
[delta1,FirstFile] = min(abs(FirstTime-min(ObsGrid.Track.Time(:))));
[delta2,LastFile]  = min(abs( LastTime-max(ObsGrid.Track.Time(:))));
if FirstFile > LastFile; 
  a = FirstFile; b = delta1; c = LastFile; d = delta2; 
  LastFile = a; delta2 = b; FirstFile = a; delta1 = d;
end


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
%define range covered by the dataset, with a bit of padding (10% of the range on each side)
Range = [floor(min(ObsGrid.Track.Lon)),ceil(max(ObsGrid.Track.Lon)), ...
         floor(min(ObsGrid.Track.Lat)),ceil(max(ObsGrid.Track.Lat))];
Range = Range ...
      + [-1,1,0,0].*0.1.*range(Range([1,2])) ...
      + [0,0,-1,1].*0.1.*range(Range([3,4]));
if Range(1) < -180; Range(1) = -180; end
if Range(2) >  180; Range(2) =  180; end
if Range(3) <  -90; Range(3) =  -90; end
if Range(4) >   90; Range(4) =   90; end


if FixedPFlag == 1;  InnerLoop = 1; TempFiles = cell(numel(Files),1); 
else;                InnerLoop = 2; TempFiles = cell(numel(Files),3); %the first dim is for the merged file
end

FileTimes = NaN(numel(Files),1);


%reinitialise rngs with the range plus the time, which is hopefully unique (often things run at the exact same time, so jsut time isn't)
rng(round(sum(Range)+datenum(now))) 
%  pause(30+randi(90))

%create a unique identifier, in case multiple processes are working at once (likely)
%combination of the time to the tenth of a millisecond and a random number between
%one and a million should hopefully be unique?
Identifier = [strrep(num2str(datenum(now)),'.',''),'_',num2str(randi(1e6,[1]))];
% Identifier = '';


for iFile=1:1:numel(Files)
  for iSource=1:1:InnerLoop %only do second loop if we want to extract true pressure
    switch iSource
      case 1; Source = TemperatureVar;
      case 2; Source = PressureVar;
    end

    disp(' ');disp(' ');
    disp(['Subsetting model file ',num2str(iFile),' of ',num2str(numel(Files)),' using CDO - variable: ',Source])


    %produce filename
    TempFiles{iFile,iSource+1} = [ScratchPath,Source,'_',sprintf('%06d',iFile),'_',Identifier,'.nc'];

    %now, generate a CDO command to subset the data down in space
    Command = 'cdo -P 16'; %16 cores permitted. Add '-v' to produce verbose output
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

    %store the times in the steps
    FileTimes(iFile,:) = TimeStore(iFile,InRange(iTime));

    %call the grid
    if isfield(Paths,'Grid')
      Command = [Command,' -setgrid,',Paths.Grid,' '];
    end


    %finalise the command with the names of the in and out files
    Command = [Command,' ',strrep(Files{iFile},'ta',Source),' ',TempFiles{iFile,iSource+1}];


    %and execute it
    if exist(TempFiles{iFile,iSource+1},'file'); continue; end %for testing with non-unique identifiers
    disp(['Executing CDO command: ',Command])
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
    if exist(TempFiles{iFile,1},'file'); continue; end %for testing with non-unique identifiers

    %and execute it
    disp(' ');disp(' ');    
    disp(['Executing CDO command: ',Command])    
    status = system(Command);
    if status ~=0;
      %error!
      Model.Error = 2;
      return
    end
    disp(['File ',num2str(iFile),' of ',num2str(numel(Files)),' merged'])
  else
    %put the link to the T file in the first element of TempFiles
    TempFiles{iFile,1} = TempFiles{iFile,2};
    TempFiles = TempFiles(1);
  end

end; clear iFile
clear Range Files ff InRange Command status PrsFile 


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
    if FixedPFlag~=1; Model.Prs = permute(ncread(TempFiles{iFile,1},PressureVar),[4,1,2,3]);
    else;             Model.Prs = FixedP;
    end
    Model.Time = FileTimes(~isnan(FileTimes(1,:)));
    Model.T    = permute(ncread(TempFiles{iFile,1},'ta'),[4,1,2,3]);
  else
    %otherwise, append
    Model.Time = [Model.Time,FileTimes(~isnan(FileTimes(iFile,:)))];
    Model.T    = cat(1,Model.T,permute(ncread(TempFiles{iFile,1},'ta'),[4,1,2,3]));
    if FixedPFlag~=1; Model.Prs  = cat(1,Model.Prs,permute(ncread(TempFiles{iFile,1},PressureVar),[4,1,2,3])); end
  end


end; clear iFile


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% delete working files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');disp(' ');
for iFile=1:1:numel(TempFiles)
 delete(TempFiles{iFile})
 disp(['Tidying up: ',TempFiles{iFile},' deleted'])
end; clear iFile



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% interpolate the data onto a common pressure scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if FixedPFlag ~= 1

  disp(' ');disp(' ');
  disp('Interpolating to single pressure scale')

  %what dim is the one we want to interpolate?
  Dim = 4;

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

%all should be correct...

%shift longitudes into range
% % % % % % % % % % % Model.Lon(Model.Lon > 180) = Model.Lon(Model.Lon > 180) - 360;

%success!
Model.Error = 0;
return

