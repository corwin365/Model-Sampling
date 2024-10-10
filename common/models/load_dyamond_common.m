function Model = load_dyamond_common(ObsGrid,WantedModel,FixedPFlag,BlobScale)


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
    GridSize    = [0.04,0.04]; %degrees lon,lat
    FixedP      = [1.57, 1.97, 2.46, 3.05, 3.78, 4.66, 5.74, 7.04, 8.60, 10.5, 12.7, 15.3, 18.3, 21.7, 25.6, 30, 34.9, 40.3, 46.2, 52.6, 59.5, 66.8, 74.5, 82.7, 91.1, 99.8, 109, 118, 127, 136, 146, 156, 166, 177, 189, 201, 214, 228, 242, 257, 273, 289, 307, 325, 344, 364, 385, 407, 429, 453, 478, 504, 530, 557, 583, 610, 637, 663, 689, 715, 740, 765, 789, 813, 835, 857, 878, 898, 916, 934, 950, 964, 977, 987, 996, 1000, 1010];
  case 'geos1p5km';
    %this doesn't work yet
    Paths.T     = '/work/bk1040/DYAMOND/data/winter_data/DYAMOND_WINTER/NASA/GEOS-1km/DW-ATM/atmos/1hr/wap/r1i1p1f1/3d/gn';
    GridSize    = [0.01,0.01]; %degrees lon,lat
    FixedP      =  [1008.3652,1006.1125,1003.6658,1001.0101,998.13165,995.01276,991.63477,987.98083,984.02991,979.76385,975.15991,970.19617,964.85162,959.10510,952.93085,946.30682,939.21198,931.62286,923.52069,914.88068,905.68677,895.92084,885.56885,874.61884,863.06305,850.89728,838.11896,824.73474,810.75177,796.18555,781.05774,765.39398,749.22852,732.59888,715.55164,698.13708,680.41010,662.42670,644.25299,625.95697,607.60797,589.27515,571.02869,552.93750,535.06641,517.47650,500.26202,483.49414,467.20297,451.38226,436.02444,421.11563,406.64743,392.61063,378.99368,365.78625,352.97913,340.56210,328.52609,316.86224,305.56061,294.61130,284.00595,273.73602,263.79272,254.16740,244.85287,235.83958,227.11940,218.68471,210.52814,202.64204,195.01849,187.64998,180.53036,173.65205,167.00713,160.59064,154.39568,148.41571,142.64375,137.07347,131.70009,126.51762,121.51954,116.70026,112.05444,107.57777,103.26315,99.107353,95.104362,91.248917,87.538086,83.963913,80.518509,77.193420,73.981903,70.876968,67.872124,64.961655,62.141167,59.405766,56.750366,54.173119,51.673130,49.250134,46.903931,44.634087,42.440006,40.321148,38.276867,36.306076,34.407806,32.580643,30.823042,29.133589,27.511152,25.954039,24.460752,23.030113,21.660688,20.350769,19.098953,17.904091,16.764410,15.678376,14.644983,13.662410,12.729040,11.843728,11.004971,10.211002,9.4605923,8.7524252,8.0846062,7.4559150,6.8651080,6.3103070,5.7903557,5.3039913,4.8493214,4.4253492,4.0306673,3.6635652,3.3233337,3.0081916,2.7167423,2.4482894,2.2010057,1.9739232,1.7659978,1.5756230,1.4023839,1.2446179,1.1013640,0.97193444,0.85479808,0.74967796,0.65521342,0.57071102,0.49539959,0.42822772,0.36900148,0.31642285,0.27038151,0.22983353,0.19446401,0.16369103,0.13698566,0.11397021,0.094084047,0.077148452,0.062610246,0.050387368,0.039958704,0.031233277,0.023860855,0.017500000,0.012000000,0.0072500003,0.0032572495];
    FixedP = FixedP(end:-1:1); %above array is in reverse order
  case 'geos3km';
    Paths.T     = '/fastdata/ka1081/DYAMOND_WINTER/NASA/GEOS-3km/DW-ATM/atmos/1hr/ta/r1i1p1f1/ml/gn';
    GridSize    = [0.025,0.025]; %degrees lon,lat
    FixedP      =  [1008.3652,1006.1125,1003.6658,1001.0101,998.13165,995.01276,991.63477,987.98083,984.02991,979.76385,975.15991,970.19617,964.85162,959.10510,952.93085,946.30682,939.21198,931.62286,923.52069,914.88068,905.68677,895.92084,885.56885,874.61884,863.06305,850.89728,838.11896,824.73474,810.75177,796.18555,781.05774,765.39398,749.22852,732.59888,715.55164,698.13708,680.41010,662.42670,644.25299,625.95697,607.60797,589.27515,571.02869,552.93750,535.06641,517.47650,500.26202,483.49414,467.20297,451.38226,436.02444,421.11563,406.64743,392.61063,378.99368,365.78625,352.97913,340.56210,328.52609,316.86224,305.56061,294.61130,284.00595,273.73602,263.79272,254.16740,244.85287,235.83958,227.11940,218.68471,210.52814,202.64204,195.01849,187.64998,180.53036,173.65205,167.00713,160.59064,154.39568,148.41571,142.64375,137.07347,131.70009,126.51762,121.51954,116.70026,112.05444,107.57777,103.26315,99.107353,95.104362,91.248917,87.538086,83.963913,80.518509,77.193420,73.981903,70.876968,67.872124,64.961655,62.141167,59.405766,56.750366,54.173119,51.673130,49.250134,46.903931,44.634087,42.440006,40.321148,38.276867,36.306076,34.407806,32.580643,30.823042,29.133589,27.511152,25.954039,24.460752,23.030113,21.660688,20.350769,19.098953,17.904091,16.764410,15.678376,14.644983,13.662410,12.729040,11.843728,11.004971,10.211002,9.4605923,8.7524252,8.0846062,7.4559150,6.8651080,6.3103070,5.7903557,5.3039913,4.8493214,4.4253492,4.0306673,3.6635652,3.3233337,3.0081916,2.7167423,2.4482894,2.2010057,1.9739232,1.7659978,1.5756230,1.4023839,1.2446179,1.1013640,0.97193444,0.85479808,0.74967796,0.65521342,0.57071102,0.49539959,0.42822772,0.36900148,0.31642285,0.27038151,0.22983353,0.19446401,0.16369103,0.13698566,0.11397021,0.094084047,0.077148452,0.062610246,0.050387368,0.039958704,0.031233277,0.023860855,0.017500000,0.012000000,0.0072500003,0.0032572495];
    FixedP = FixedP(end:-1:1); %above array is in reverse order
  case 'ifs4km';
%      Paths.Root  = '/home/b/b382226/DYAMOND_winter/ECMWF/IFS-4km/DW-CPL/atmos/';
    Paths.T     = '/scratch/b/b382226/ifs_4'; %data needs to be pre-converted from spherical harmonics to a reduced gaussian grid using CDO, and then put here.
    GridSize    = [0.03,0.03]; %degrees lon,lat
    FixedP      = [0.0100018250000000, 0.0255130300000000,0.0388416250000000,0.0574703050000000,0.0828747150000000,0.116761950000000,0.161071775000000,0.217973245000000,0.289857140000000,0.379324760000000,0.489173525000000,0.622380195000000,0.782082290000000,0.971558115000000,1.19420624000000,1.45352455500000,1.75308983000000,2.09653755000000,2.48754264500000,2.92980163500000,3.42701553500000,3.98287430000000,4.60104263500000,5.28514740000000,6.03876678500000,6.86542023000000,7.76855987500000,8.75156372000000,9.81773041000000,10.9702740500000,12.2123199500000,13.5469024700000,14.9769653350000,16.5053595000000,18.1348413100000,19.8680761700000,21.7076379400000,23.6560107450000,25.7155932650000,27.8886975100000,30.1775537100000,32.5843151850000,35.1110595700000,37.7597937050000,40.5321057150000,43.4287402350000,46.4498339850000,49.5952221700000,52.8644287150000,56.2566772500000,59.7720947300000,63.4151049800000,67.1940869100000,71.1186987300000,75.1999062500000,79.4493305650000,83.8808969750000,88.5101289100000,93.3509150400000,98.4128056850000,103.704051800000,109.232469750000,115.005147450000,121.030576150000,127.316604500000,133.871710000000,140.704500050000,147.822853550000,155.235774400000,162.952860350000,170.982860350000,179.334187500000,188.017056650000,197.041511700000,206.416171850000,216.151273400000,226.257761700000,236.745539050000,247.624503900000,258.905773450000,270.600168000000,282.718611350000,295.271978500000,308.271136700000,321.727359350000,335.652433600000,350.057984400000,364.955175800000,380.356140650000,396.272988300000,412.717556650000,429.701488300000,447.236994150000,465.336437500000,484.011839850000,503.251664050000,522.990328100000,543.111136700000,563.497671900000,584.051710950000,604.672742200000,625.260386750000,645.714968750000,665.940136700000,685.844488250000,705.343222650000,724.359402350000,742.824457025000,760.679343750000,777.875140625000,794.372593750000,810.141828125000,825.162265625000,839.422070315000,852.917398440000,865.651859375000,877.635179690000,888.882476565000,899.413000000000,909.250531250000,918.420914065000,926.952710940000,934.876882815000,942.223617190000,949.024484375000,955.310789065000,961.114500005000,966.466273440000,971.394882815000,975.929953130000,980.099804690000,983.929929690000,987.445414065000,990.670085940000,993.626468755000,996.335289065000,998.815000000000];    
    FixedP = FixedP(1:80); %the input files have been trimmed to the first 80 levels to save space and reduce runtime
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


%first, get a list of available files for this model
Files = wildcardsearch(Paths.T,'*.nc');

%identifying the timestamps of all the files can be slow, so let's use a lookup file unless the list of files present has changed
Lookup = [ScratchPath,'/',WantedModel,'.mat'];
if exist(Lookup,'file')
  %check it contains what we need
  A = load(Lookup);
  %generate a checksum of the list of file names
  if sum(getByteStreamFromArray(Files)) == sum(getByteStreamFromArray(A.Files))
    disp('Using timestamps from lookup file')
    Files = A.Files;
    TimeStore = A.TimeStore;
  end
  clear A;
end

if ~exist('TimeStore','var')

  disp(['Generating timestamp lookup file for ',WantedModel]);
  for iFile=1:1:numel(Files)
    disp(Files{iFile})

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
    TimeStore(TimeStore == 0) = NaN;

    save(Lookup,'TimeStore','Files');
    
  end; clear iFile
else; load(Lookup);
end
clear Lookup




%hence, find the files which cover the period in our data
%we just need the closest timestep to the start and end and don't need to pad it
%this is because the sampling won't interpolate, so only the closest points are used
[~,idx] = min(abs(nanmin(ObsGrid.Track.Time,[],'all')-TimeStore(:))); [FirstFile,~] = ind2sub(size(TimeStore),idx);
[~,idx] = min(abs(nanmax(ObsGrid.Track.Time,[],'all')-TimeStore(:))); [LastFile, ~] = ind2sub(size(TimeStore),idx);

Files = Files(FirstFile:1:LastFile);
TimeStore = TimeStore(FirstFile:1:LastFile,:);
clear idx FirstFile LastFile 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, use cdo to cut down the data to take these files and cut them down
%to just the geographic area we need
%this could take a while, but it's much faster than trying to load the whole 
%dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we're going to do this twice - once for pressure and once for temperature - and then merge them at the end


%define range covered by the dataset in terms of measurement centre locations
Range = NaN(5,4);
Range(1,:) = [floor(min(ObsGrid.Track.Lon)),ceil(max(ObsGrid.Track.Lon)), ...
              floor(min(ObsGrid.Track.Lat)),ceil(max(ObsGrid.Track.Lat))];

%find the maximum extra distance a blob could add to the edge of the region (with 25% padding)
MaxDistance = max([ObsGrid.Weight.X;ObsGrid.Weight.Y]).*BlobScale.*1.25;

%find where this puts the range, by taking a line diagonally from each corner
Range(2,[4,2]) = reckon(Range(1,4),Range(1,2),km2deg(MaxDistance),45);
Range(3,[3,2]) = reckon(Range(1,3),Range(1,2),km2deg(MaxDistance),135);
Range(4,[3,1]) = reckon(Range(1,3),Range(1,1),km2deg(MaxDistance),225);
Range(5,[4,2]) = reckon(Range(1,4),Range(1,1),km2deg(MaxDistance),315);

r(1) = min(Range(:,[1,2]),[],'all');
r(2) = max(Range(:,[1,2]),[],'all');
r(3) = min(Range(:,[3,4]),[],'all');
r(4) = max(Range(:,[3,4]),[],'all');
Range = r; clear r BlobScale MaxDistance

if Range(1) < -180; Range(1) = -180; end
if Range(2) >  180; Range(2) =  180; end
if Range(3) <  -90; Range(3) =  -90; end
if Range(4) >   90; Range(4) =   90; end



if FixedPFlag == 1;  InnerLoop = 1; TempFiles = cell(numel(Files),1); 
else;                InnerLoop = 2; TempFiles = cell(numel(Files),3); %the first dim is for the merged file
end

FileTimes = NaN(numel(Files),1);


%reinitialise rngs with the range plus the time, which is hopefully unique (often things run at the exact same time, so just time isn't)
rng(round(sum(Range)+datenum(now))) 
%  pause(30+randi(90))

%create a unique identifier, in case multiple processes are working at once (likely)
%combination of the time to the tenth of a millisecond and a random number between
%one and a million should hopefully be unique?
Identifier = [strrep(num2str(datenum(now)),'.',''),'_',num2str(randi(1e6,[1]))];
%  Identifier = '';


for iFile=1:1:numel(Files)
  for iSource=1:1:InnerLoop %only do second loop if we want to extract true pressure
    switch iSource
      case 1; Source = 'ta';
      case 2; Source = PressureVar;
    end

    disp(' ');disp(' ');
    disp(['Subsetting model file ',num2str(iFile),' of ',num2str(numel(Files)),' using CDO - variable: ',Source])

   
    %produce filename
    TempFiles{iFile,iSource+1} = [ScratchPath,Source,'_',sprintf('%06d',iFile),'_',Identifier,'.nc'];
    
    %if we've disabled identifiers, and the file exists, then just load the old file
    %we still need to generate a time value for it
    if numel(Identifier) == 0 && exist(TempFiles{iFile,iSource+1},'file');
      FileTimes(iFile,:) = TimeStore(iFile);
      continue
    end 

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
    InRange = find(TimeStore(iFile,:) <= nanmax(ObsGrid.Track.Time(:)-3/24) ...
                 & TimeStore(iFile,:) >= nanmin(ObsGrid.Track.Time(:)+3/24));
    if numel(InRange) == 0; [~,InRange] = min(abs(TimeStore(iFile,:)-nanmean(ObsGrid.Track.Time(:)))); end

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
    if iFile == numel(Files);TempFiles = TempFiles(:,1); end
  end

end; clear iFile
clear Range ff InRange Command status PrsFile 

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

if numel(Identifier) ~= 0
  disp(' ');disp(' ');
  for iFile=1:1:size(TempFiles,1)
   delete(TempFiles{iFile,1})
   disp(['Tidying up: ',TempFiles{iFile},' deleted'])
  end; clear iFile

end

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
Model.Lon(Model.Lon > 180) = Model.Lon(Model.Lon > 180) - 360;

%success!
Model.Error = 0;
return

