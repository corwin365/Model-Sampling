function [Error,OldData] = sampling_core_v2(Instrument,ModelName,DayNumber,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sample models and reanalyses with satellite scan tracks
%
%originally written 2019, extensively updated in February 2023 to make more
%user-friendly and include several formerly-manual steps automatically. Not 
%backwards-compatible to before this rewrite.
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/02/03
%
%possible error states:
%-1. unknown error (shouldn't occur)
%0. successful
%1. problem loading daily obs file
%2. problem loading model file
%3. already done, and clobber not set
%4. sensitivity testing, and output path not set
%5. instrument not supported
%6. model not supported
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get functions
addpath(genpath('common/'));

%default error handling
Error = -1; OldData.ModelID = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% input parsing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%create input parser and testing functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
IsPositiveInteger = @(x) validateattributes(x,{'numeric'},{'positive','integer'});
IsPositive        = @(x) validateattributes(x,{'numeric'},{'positive'});
IsNonNegative     = @(x) validateattributes(x,{'numeric'},{'>=',0});

%inputs - required
%%%%%%%%%%%%%%%%%%%

addRequired(p,'Instrument',  @ischar);  %just a type check - will also be checked against the list of valid options below
addRequired(p,'ModelName',   @ischar);  %just a type check - will also be checked against the list of valid options below
addRequired(p,'DayNumber',   @isnumeric);

%inputs - optional
%%%%%%%%%%%%%%%%%%%

%flags
addParameter(p, 'SensTestMode',  false,        @islogical);     %use sensitivity-testing mode? Note that this overrides many of these settings
addParameter(p,      'Clobber',  false,        @islogical);     %overwrite previous output files (assume no)
addParameter(p,     'Parallel',  false,        @islogical);     %parallel or linear mode? 
addParameter(p,   'TextUpdate',  false,        @islogical);     %display progress of inner sampling loop to screen in linear mode? 
addParameter(p,  'IncludeNaNs',  true,         @islogical);     %if a measurement blob includes NaNs or goes off the edge of the field, do we include it in the result? Useful for e.g. limited-area model runs; if set, output field can include partially-sampled blobs

%choices we need to make
addParameter(p,   'ReportEvery', 10000, IsPositiveInteger);     %update to screen how often (in samples).
addParameter(p,       'MinPrs' ,     0,        IsPositive);     %minimum pressure level, i.e. top height. If set to zero, uses model top
addParameter(p,       'MaxPrs' ,  1100,        IsPositive);     %maximum pressure level, i.e. bottom height.

%data selection
addParameter(p,       'SubSet',      0, IsPositiveInteger);     %which subfile within the specified day to work on - used e.g. for AIRS granule numbers
addParameter(p,   'HoursAhead',      0,     IsNonNegative);     %if using forecast data, how many hours in advance?

%passed-through data
addParameter(p,     'OldData','NOTSET',          @isstruct);    %do we have a previously-used set of model interpolants in memory? 
addParameter(p, 'Sensitivity','NOTSET',          @isstruct);    %parameters for sensitivity-testing mode

%paths
addParameter(p,'DensityPath',        './common/saber_density_filled.mat',@isfile); %path to density data
addParameter(p, 'MasterPath', [LocalDataDir,'/corwin/sampling_project/'], @isdir); %path to density data
addParameter(p,    'OutPath',                                   'NOTSET', @isdir); %output file. will be generated automatically below if set to 'NOTSET' or not supplied
 
%arbitrary values used in the code. Most defaults have been selected via sensivity testing using 3D AIRS data.
addParameter(p,    'MinSignal' ,  0.99,        IsPositive);     %fraction of total signal needed to produce final sample
addParameter(p,    'BlobScale' ,     3,        IsPositive);     %number of standard deviations to compute sensitivity out to (+- from centre)
addParameter(p,   'MinZContrib',  0.02,        IsPositive);     %when rotating, discard vertical levels contributing less fractional weight than this
addParameter(p,      'ZPadding',   0.5,        IsPositive);     %when rotating, vertical padding in decades of pressure

%parse inputs
%%%%%%%%%%%%%%%%%%%%%%%

parse(p,Instrument,ModelName,DayNumber,varargin{:});
Settings = p.Results;

%do some preprocessing based on the above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%is the observational dataset on the valid list?
[~,InstList] = instrument_settings('JUSTLISTING',Settings);
if sum(ismember(InstList,Instrument)) ~= 1;
  disp(['Instrument ''',Instrument,''' not yet supported, stopping'])
  disp(['Valid (case-sensitive) list is:'])
  for iInst=1:1:numel(InstList)
    disp(['--> ',InstList{iInst}])
  end  
  Error = 5; return
end; 

%is the model on the valid list?
ModelList = model_settings([],'JUSTLISTING');
if sum(ismember(ModelList,ModelName)) ~= 1;
  disp(['Model ''',ModelName,''' not yet supported, stopping'])
  disp(['Valid (case-sensitive) list is:'])
  for iModel=1:1:numel(ModelList)
    disp(['--> ',ModelList{iModel}])
  end
  Error = 6; return
end; 

%display what we think we're doing
disp(['This routine is set to sample ',ModelName,' as ',Instrument,' on ',datestr(DayNumber),'.']);

%address sub-daily subsets of the data (e.g. AIRS granules)
if Settings.SubSet == 0; 
  SubSetOutString = '';
else           
  SubSetOutString = ['_subset',sprintf('%06d',Settings.SubSet)];
  disp(['Processing subset ',sprintf('%06d',Settings.SubSet)])
end

%if we're using forecast data, we need to specify how far ahead a forecast
%we want. This will be used by the model selection subroutine to select
%what data to feed back to the sampling parent
if Settings.HoursAhead ~= 0;
  disp(['Using ',num2str(Settings.HoursAhead),' hour forecast'])
  ForecastOutString = ['_fc',sprintf('%06d',Settings.HoursAhead),'hrs'];
else
  ForecastOutString = '';
end

%do we have a previously-used set of model interpolants in memory? if not, create an empty checking variable
if ischar(Settings.OldData); clear OldData; OldData.ModelID = ''; end

%are we including partial blobs that contain NaNs in their volume?
if Settings.IncludeNaNs ~=1;
  disp('Blobs containing any NaNs will be set to NaN');
end


%where shall we put our results? 
if strcmp(Settings.OutPath,'NOTSET')
  Settings.OutPath = [Settings.MasterPath,'/output/',Instrument,'/',ModelName,'/sampled_',num2str(DayNumber),SubSetOutString,ForecastOutString,'.mat'];

  %make sure this directory exists!
  if exist([Settings.MasterPath,'/output/',Instrument              ],'dir') ~= 7; mkdir([Settings.MasterPath,'/output/',Instrument              ]); end
  if exist([Settings.MasterPath,'/output/',Instrument,'/',ModelName],'dir') ~= 7; mkdir([Settings.MasterPath,'/output/',Instrument,'/',ModelName]); end
end


%check if we've already done this day
if exist(Settings.OutPath,'file') ~= 0 && Settings.Clobber == 0
  disp([ModelName,' as ',Instrument,' on ',datestr(DayNumber),' already exists and Clobber is not set, stopping']);
  Error = 3; return;
end

%finally, store the time we started the analysis (especially useful for sensitivity testing!)
RunTime = struct();
RunTime.Start = now;

%tidy
clear IsPositive IsPositiveInteger IsNonNegative p varargin InstList ModelList


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load observation properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get instrument options
%this includes the FineGrid spacings needed below, and the SubetInStrings needed to load multiple-file-per-day datasets
Settings = instrument_settings(Instrument,Settings);

%now identify the file containing the daily grid
ObsGrid = [Settings.MasterPath,'/tracks/',Instrument,'/track_',Instrument,'_',num2str(DayNumber),Settings.SubSetInString,'.mat'];

%and load it
try; ObsGrid = load(ObsGrid);
catch;  %failed.
  disp(['Cannot locate file containing observational track: ',ObsGrid]);
  Error = 1; return
end
%convert obs grid to log space, to make spacings more regular
ObsGrid.Track.Prs = log10(ObsGrid.Track.Prs);
  
disp([Instrument,' track loaded for ',datestr(DayNumber)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep model data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%identify and load the model
%file conventions vary massively, so load the data inside the switch
%child functions will format the data into a single standard
%child output format:
%  struct called Model
%  containing fields:
%  Lon   - 1d. Runs from -180 to +180
%  Lat   - 1d
%  Time  - 1d
%  Prs   - 1d
%  T     - 4d, time x lon x lat x pressure
%this routine then produces a 3d interpolant field for each timestep based on this

%firstly, do we need to bother, or did we already load it?
if strcmp(ModelName,OldData.ModelID) == 1 ...
       && DayNumber >= min(OldData.Time)  ...
       && DayNumber <  max(OldData.Time) %this needs to be < rather than <= as we may have the first timestep of the day in memory from the last day
   
  %same day as last loop - just call the data back up from memory
  Interpolants = OldData.I;
  disp(['Using previously loaded data and interpolants for ',ModelName,' on ',datestr(DayNumber)]);
                 
else
     
  %clear any previous model data we had stored - saves memory
  clear OldData  
  
  %load the model
  [Model,OldData,Settings] = model_settings(DayNumber,ModelName,Settings,ObsGrid);

  %check if we loaded it successfully
  if Model.Error ~=0
    OldData.ModelID = ''; %we need to set this for the return
    disp('Problem loading model data. Stopping')
    Error =2;
    return
  end

  %convert pressure to log-prs, to make maths easier
  Model.Prs = log10(Model.Prs);
  
  %drop parts of the atmosphere that we don't care about
  InPrsRange = find(Model.Prs <= log10(Settings.MaxPrs));
  Model.Prs = Model.Prs(InPrsRange);
  Model.T   = Model.T(:,:,:,InPrsRange);
  clear InPrsRange
  
  %make sure pressure increases in z. needed for griddedinterpolant below.
  if mean(diff(Model.Prs)) < 0
    Model.Prs = flip(Model.Prs,1);
    Model.T   = flip(Model.T,  4);
  end
  
  %%finally, produce gridded interpolant fields for each timestep
  [x,y,z] = ndgrid(Model.Lon,Model.Lat,Model.Prs);
  Interpolants = struct();
  for iTime=1:1:numel(Model.Time)
    F = griddedInterpolant(single(x),single(y),single(z),single(squeeze(Model.T(iTime,:,:,:))));
    Interpolants.(['t',num2str(iTime)]) = F;
  end
  Interpolants.Time = Model.Time;
  clear x y z iTime F

  disp([ModelName,' interpolants prepared for sampling on ',datestr(DayNumber)]);
  
  %finally, store the model ID, time period, pressure and interpolants for future use
  OldData.I = Interpolants;
  OldData.ModelPressureScale = Model.Prs;
  OldData.ModelTimeScale = Model.Time;  
  OldData.ModelID = ModelName;  
  OldData.Time = [min(Model.Time),max(Model.Time)];

  %then dump the model data to save memory
  %we already have the axes and the interpolants, and that's all we want  
  clear Model 
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% That was complicated. Here's an ASCII satellite to relax us again.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                         ooo
%                       / : \
%                       / o0o \
%                 _____"~~~~~~~"_____
%                 \+###|U * * U|###+/
%                 \...!(.>..<)!.../
%                   ^^^^o|   |o^^^^
%               #+=====}:^^^^^:{=====+#
%                .____  .|!!!|.  ____.
%                |#####:/" " "\:#####|
%                |#####=|  O  |=#####|
%                |#####>\_____/<#####|
%                 ^^^^^   | |   ^^^^^
%                         o o
%               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OK, back to work. 
%We want to do sensitivity testing of our grid shapes
%but we also want to use the standard core so bugs don't propagate through
%duplicated and then modified code
%so, this section overrides primary grid-defining variables if we're in
%sensitivity-testing mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.SensTestMode ~= 0

  Sensitivity = Settings.Sensitivity;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %sensitivity-testing mode is active. change grid parameters and output path
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %print out what settings are in force
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp('**********************************************') 
  disp('***********Sensitivity testing mode***********')  
  disp(Sensitivity )
  if isfield(Sensitivity,'FineGrid'); Sensitivity.FineGrid
  end
  disp('**********************************************') 
  disp('**********************************************')  
  
  %first, overwrite the gridding parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if isfield(Sensitivity,'BlobScale'); Settings.BlobScale = Sensitivity.BlobScale; end
  if isfield(Sensitivity,'MinSignal'); Settings.MinSignal = Sensitivity.MinSignal; end
  if isfield(Sensitivity,'MaxPrs');    Settings.MaxPrs    = Sensitivity.MaxPrs;    end
  
  if isfield(Sensitivity,'FineGrid')
    if isfield(Sensitivity.FineGrid,'X');   Settings.FineGrid(1) = Sensitivity.FineGrid.X;   end
    if isfield(Sensitivity.FineGrid,'Y');   Settings.FineGrid(2) = Sensitivity.FineGrid.Y;   end
    if isfield(Sensitivity.FineGrid,'Prs'); Settings.FineGrid(3) = Sensitivity.FineGrid.Prs; end    
  end
  
  %then, change the output path
  %this is required as otherwise it'll overwrite real results
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isfield(Sensitivity,'NewPath')
    Error = 3;
    disp('Sensitivity test output file not set')
    return;
  end
 
  OutPath = [Settings.MasterPath,Sensitivity.NewPath]; 
  Settings.DensityPath = Sensitivity.DensityPath; %we might be running it from a different directory

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now go down to the datapoint level and start sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%final steps before core loop:

%1. if we're rotating in the vertical plane, we need density data.
if sum(ObsGrid.Track.ViewAngleZ,'omitnan') ~= 0
  %load density data (global-mean daily-mean saber-derived density)
  Density = struct();
  Density = load(Settings.DensityPath);
  Density = Density.Results;
  Density.Prs = log10(h2p(Density.HeightScale./1000));
else
  Density = []; %not used, don't waste time making it
end

%2. create a results array and an array for the 'naive' approximation where we just sample the model at measurement centre
FinalT = ObsGrid.Track.Lon.*NaN; 
SimpleT = FinalT;

%3. if we didn't specify a minimum pressure to use, compute this
if Settings.MinPrs == 0; Settings.MinPrs = 10.^min(OldData.ModelPressureScale); end

%done! let's go!
disp(['Commencing sampling, ',num2str(numel(FinalT)),' samples to extract']) 

if Settings.Parallel == 0;  
  %simple for loop. Easier on system memory.
  disp('Using single-threaded mode')
  for    iSample = 1:numel(FinalT); 
    [~,FinalT(iSample),SimpleT(iSample)] = innercore(iSample,ObsGrid,Interpolants,Settings,Density); 
  end
elseif Settings.Parallel == 1; 
  %parfor loop. Often much faster, but needs more memory so won't work for large datasets.
  disp('Using parallelised mode')
  parfor iSample = 1:numel(FinalT);
    [~,FinalT(iSample),SimpleT(iSample)] = innercore(iSample,ObsGrid,Interpolants,Settings,Density);
  end
end

disp(['Sampling complete for ',ModelName,' as ',Instrument,' on ',datestr(DayNumber)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%success! save and return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%call the reconstructor to convert from a list of points to a formatted set of fields
[~,Sampled_Data] = unified_reconstructor(ObsGrid,FinalT,SimpleT);

%save end time so we know how long it all took (especially useful for sensitivity testing!)
RunTime.End = datenum(now);

%write the output file
save(Settings.OutPath,'Sampled_Data','Settings','RunTime');
disp('Saved!');

%and we're done
Error = 0; %we didn't fail!
return


end% this is the end of the main function




















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% inner core
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [InnerError,TSample,TSimple] = innercore(iSample,ObsGrid,Interpolants,Settings,Density)

  %{
  This is the inner loop that does the actual sampling.

  This doesn't need to be a function - it's part of the main programmatic flow - but
  doing it this way allows programmatic switching between parallel and linear flow
  depending on data volume.

  As it's not a "real" function, there's no testing of inputs - I assume the outer
  programme got it right. The errors it spits out are as follows:
  
  1: unknown failure (shouldn't happen)
  2: at least one input is NaN
  3: point is over the top of the model data
  4: point is below model bottom

  All except the first one are normal and routine, and just result in a NaN in the final
  results array. None are read out of the routine into the parent normally (again, we're
  assuming a high level of trust between parent and child), but exist for debugging.
  %}

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  %set some placeholder values
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  %output handling
  InnerError = 1; %assume an error unless the function runs successfully
  TSample = NaN;
  TSimple = NaN;

  %internal use
  Channel = [];
  Wx      = [];   
  Wy      = [];  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %upate progress to screen
  %(only on request, and only in single-threaded mode)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  if Settings.Parallel == 0 && Settings.TextUpdate == 1
    if iSample == 1; tic; end
    if mod(iSample,Settings.ReportEvery) == 0; 
      disp(['Processed ',num2str(round(iSample./numel(ObsGrid.Track.Time).*100.*100)./100),'%; time since last update ',num2str(toc),'s']);tic; end
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %let's go!
  %first, where are we?
  %get geolocation and measurement properties
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Sample = struct();
  Sample.Lon     = ObsGrid.Track.Lon(       iSample);
  Sample.Lat     = ObsGrid.Track.Lat(       iSample);
  Sample.Prs     = ObsGrid.Track.Prs(       iSample); %remember, we already log-ged this
  Sample.Time    = ObsGrid.Track.Time(      iSample);
  
  Sample.AngleH  = ObsGrid.Track.ViewAngleH(iSample);
  Sample.AngleZ  = ObsGrid.Track.ViewAngleZ(iSample);

  Sample.WeightX = ObsGrid.Weight.X(        iSample);
  Sample.WeightY = ObsGrid.Weight.Y(        iSample);
  Sample.WeightZ = ObsGrid.Weight.Z(        iSample);
  
  
  %safety-check that none of the above are NaN
  if isnan(Sample.Lon     + Sample.Lat     + Sample.Prs + Sample.Time ...
         + Sample.AngleH  + Sample.AngleZ                             ...
         + Sample.WeightX + Sample.WeightY + Sample.WeightZ           ...      
          ); InnerError = 2;return; 
  end
   
  %check we're not over the top of the model
  if Sample.Prs < log10(Settings.MinPrs); InnerError = 3;return; end

  %or below the bottom 
  if Sample.Prs > log10(Settings.MaxPrs); InnerError = 4;return; end  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %produce the fine sampling grid
  %this is the grid the actual sampling will be done
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %declare a struct to keep the grids in
  Fine = struct();
  
  %creating the grid in the horizontal is easy
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Fine.X = -Settings.BlobScale*Sample.WeightX : Settings.FineGrid(1) : Settings.BlobScale*Sample.WeightX;
  Fine.Y = -Settings.BlobScale*Sample.WeightY : Settings.FineGrid(2) : Settings.BlobScale*Sample.WeightY;
  
  %for the vertical, we're working in pressure co-ordinates, so this is fiddly
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %tmost instruments are a simple Gaussian in height
  %these have *POSITIVE* Sample.WeightZ values
  %so create a height window accordingly
  if Sample.WeightZ > 0
    %first, find approx height of the measurement-centre pressure level
    z = p2h(10.^Sample.Prs);    
    %next, find height above and below this that we want, and convert to log-P
    z = (Settings.BlobScale*[-1,1].*Sample.WeightZ)+z;
    z = log10(h2p(z));  
  else
    %some instruments are more complex than just a simple Gaussian
    %these are specified with a *NEGATIVE* Sample.WeightZ
    %this |number| tells us which of a pre-defined set of functions to use
    %from those stored in the source file
    
    %identify which channel we want.
    Channel = struct();
    Channel.ID = abs(Sample.WeightZ);
    
    Channel.Prs = log10(ObsGrid.Weight.ZFuncs.PrsScale);
    Channel.W   = ObsGrid.Weight.ZFuncs.Weights(Channel.ID,:); %keep for later  

    %discard parts that contribute very little
    Channel.W(Channel.W < Settings.MinZContrib.*sum(Channel.W(:)),'omitnan') = 0; %less than 2% is a good balance of useful volume and runtime (tested for AIRS 42km)

    %set NaNs to zero
    Channel.W(isnan(Channel.W)) = 0;
    
    %remove low altitudes
    GoodPrs = find(Channel.Prs < log10(Settings.MaxPrs));
    Channel.W   = Channel.W(   GoodPrs);
    Channel.Prs = Channel.Prs(GoodPrs);
    
    %find highest and lowest prs level remaining
    z = Channel.Prs(find(Channel.W > 0));
    z = [max(z),min(z)]+[1,-1].*Settings.ZPadding;
  end
  Fine.Prs = single(z(2):Settings.FineGrid(3):z(1));


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %generate the weighting blob
  %we're going to use the following coord frame for this:
  % X is along-track (i.e. major volume axis in horizontal plane)
  % Y is across-track (i.e minor axis)
  % Z is vertical
  %we'll later rotate the coord frame appropriately for the viewing angle
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %1d along-track weight
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if Sample.WeightX > 0  %simple Gaussian approximation
    Wx = normpdf(Fine.X,0,Sample.WeightX);
  else %not a simple Gaussian - produce from supplied information
    stop %not written yet, as no such cases implemented
  end  
  
  %1d across-track weight
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if Sample.WeightY > 0  %simple Gaussian approximation
    Wy = normpdf(Fine.Y,0,Sample.WeightY);
  else %not a simple Gaussian - produce from supplied information
    stop %not written yet, as no such cases implemented 
  end  

  %1d vertical weight
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if Sample.WeightZ > 0  %simple Gaussian approximation
    z = p2h(10.^Fine.Prs);  
    Wz = normpdf(z,p2h(10.^Sample.Prs),Sample.WeightZ);
  else %not a simple Gaussian - interpolate from supplied information
    Wz = interp1(p2h(10.^Channel.Prs),Channel.W,p2h(10.^Fine.Prs),'linear','extrap');
  end  
  
  %multiply them together to create a 3D blob
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %order of indexing will seem odd if you're unfamiliar with Matlab
  %don't blame me...
  Wx = repmat(Wx',1,numel(Fine.Y),numel(Fine.Prs));
  Wy = repmat(Wy,numel(Fine.X),1,numel(Fine.Prs));
  Wz = repmat(permute(Wz',[2,3,1]),numel(Fine.X),numel(Fine.Y),1); 
  W = Wx.*Wy.*Wz;
  [Fine.Y,Fine.X,Fine.Prs] = meshgrid(Fine.Y,Fine.X,Fine.Prs);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %cut data down to reduce interpolation requirements
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  %right, now we want to choose the minimum number of points that will give us
  %our required weight.
  %this is a time bottleneck in our routine here, but is saving time overall
  %because it reduces the points done in all subsequent steps
  
  %first, sort the weight points
  [WSort,Widx]  = sort(W(:),'descend');
  %find the cumulative sum
  WSort = cumsum(WSort);
  %find where it crosses our required level
  [~,Cidx] = min(abs(WSort-max(WSort).*Settings.MinSignal));
  %and keep only these points
  Points = Widx(1:Cidx);
 
  %pull out these points out of the weight and location arrays
  W        = W(Points);
  Fine.X   = Fine.X(Points);
  Fine.Y   = Fine.Y(Points);  
  Fine.Prs = Fine.Prs(Points); 

  %normalise remaining weights to sum to 1
  W = W./sum(W(:),'omitnan');  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %rotate to appropriate viewing angle
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
  %horizontal plane
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %we want to go from:
  %X along-track, Y across track
  %to
  %X zonal, Y meridional
  %and we have the angle of the measurement c/w from N
  
  if Sample.AngleH ~= 0
    
    A = [Fine.X,Fine.Y]';
    
    R = [cosd(90-Sample.AngleH), -sind(90-Sample.AngleH); ...
         sind(90-Sample.AngleH),  cosd(90-Sample.AngleH)]';
    B = R*A;
    
    Fine.X = B(1,:)';
    Fine.Y = B(2,:)';
    
    %done.
  end
  
  %vertical plane
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  
  if Sample.AngleZ ~= 0
    
    
    %this is more complex. 
    %we have to:
    %1. shift to centre on z=0, to avoid awkward shift terms in the rotation
    %2. carry out the rotation
    %3. undo the shift
    %4. correct for density
    %angles are defined a/c/w from nadir
    
    %extract correct day of density
    [~,dayidx] = min(abs(Density.DayScale-Sample.Time));
    D = Density.Density(dayidx,:);    
     
    %0. retain old pressure coords, for density scaling later
    OldPrs = Fine.Prs;
    
    %1. shift centre (retain location to put it back)
    z = p2h(10.^Fine.Prs);
    [~,centroid] = max(W(:));
    Centre = z(centroid);
    z = z-Centre;
   
    %2. carry out the rotation
    A = [Fine.X,Fine.Y,z]';
    R = [1,                   0,                    0;...
         0,cosd(-Sample.AngleZ),-sind(-Sample.AngleZ);...
         0,sind(-Sample.AngleZ), cosd(-Sample.AngleZ)];
    B = R*A;
    
    Fine.X = B(1,:)';
    Fine.Y = B(2,:)';
    z      = B(3,:)';
    
    %3. undo the shift
    z = z+Centre;
    Fine.Prs = log10(h2p(z));

    %now, correct for density
    %because the bulk of the weight comes from a given density region
    %and flipping round like this changes that

    %we need to:
    %a. unscale the original points by density (because there's an implicit
    %scaling already in the weights we inputted)
    %b. rescale for the new 
    Unscale = interp1(Density.Prs,D,OldPrs);
    Rescale = interp1(Density.Prs,D,Fine.Prs);
    W = W ./ Unscale .* Rescale;

    %renormalise
    W = W./sum(W(:),'omitnan');
    
    %done!
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %now we need to convert these to latitude and longitude locations
  %so we can interpolate to them from the global interpolant
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  %distance from sample-centre in deg of arc and angle in degrees a/c/w from E
  Distance = distdim(quadadd(Fine.X,Fine.Y),'km','deg');
  Direction = atan2d(Fine.Y,Fine.X);
  
  %convert direction to c/w from N
  Direction = wrapTo180(-Direction + 90);
  
  %and hence find the lat and lon of each point
  [Fine.Lat,Fine.Lon] = reckon(Sample.Lat,Sample.Lon,Distance,Direction);

  %I'm getting a very weird bug in some COSMIC profiles near the poles where
  %Fine.Lat is sometimes coming out complex, with a |real| component of 90 and 
  %a near-zero imaginary component (< 0.03)
  %%if this occurs, force it to be real
  if ~isreal(Fine.Lat); Fine.Lat = real(Fine.Lat); end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %interpolate temperature onto the selected points
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %ok, we're there! choose the right interpolant object
  [~,idx] = min(abs(Sample.Time-Interpolants.Time));
  I = Interpolants.(['t',num2str(idx)]);
  
  %and apply it
  Fine.T = I(single(Fine.Lon(:)),single(Fine.Lat(:)),single(Fine.Prs(:)));
  
  %scale by weights, and store
  if Settings.IncludeNaNs == 1;
    TSample = sum(flatten(Fine.T).*flatten(W),'omitnan');
  else
    TSample = sum(flatten(Fine.T).*flatten(W));
  end
  
  %also compute simple T at measurement centre
  TSimple = I(single(Sample.Lon),single(Sample.Lat),single(Sample.Prs));

  %success!
  InnerError = 0;

  end %this is the end of the inner (sampling core) function








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reconstructor - puts the data into the right shape
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Error,Output] = unified_reconstructor(ObsGrid,Final,Simple)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%converts sampled results from list of points to a structure in 
%the same format as the original data
%
%Corwin Wright, c.wright@bath.ac.uk, 03/Feb/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Error = 1; %assume failure unless proved otherwise

%what size should the output be?
DimList = fieldnames(ObsGrid.Recon); NDims = numel(DimList); OutputSize = NaN(NDims,1);
for iDim=1:1:NDims; OutputSize(iDim) = max(ObsGrid.Recon.(DimList{iDim}));end; clear iDim

%order list of dimensions from largest number of elements to smallest (arbitrary, but I want something stable)
[OutputSize,idx] = sort(OutputSize,'desc');DimList = DimList(idx); clear idx


%work out where each data point goes
List = ObsGrid.Recon.(DimList{1});
for iDim=2:1:NDims;  List = [List,ObsGrid.Recon.(DimList{iDim})]; end; clear iDim
[~,Order]  = sortrows(List,NDims:-1:1);

%reshape the metadata
Output = struct();
Fields = {'Lat','Lon','Time','Prs'};
for iField=1:1:numel(Fields)
  Var = ObsGrid.Track.(Fields{iField});
  Var = reshape(Var(Order),OutputSize');
  Output.(Fields{iField}) = Var;
end; clear iField Var

%reshape the output data
Output.T       = reshape( Final(Order),OutputSize');
Output.Tsimple = reshape(Simple(Order),OutputSize');

%done
Error = 0;
return

end

