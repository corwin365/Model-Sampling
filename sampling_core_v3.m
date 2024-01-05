function [Error,OldData] = sampling_core_v3(Instrument,ModelName,DayNumber,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sample models and reanalyses with satellite scan tracks
%
%v1: originally written 2019, 
% 
%v2: extensively updated in February 2023 to make more user-friendly and include 
%    several formerly-manual steps automatically. Not backwards-compatible to 
%    before this rewrite.
%
%v3: modified heavily again to allow for more complex weighting functions. Should
%    be backwards compatible to v2, but central logic has changed a fair bit.
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2023/11/11
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

%the master file path is not optional as it is used upstream. Set it here.
%this is where all the track-input and sampled-output files live
MasterPath = [LocalDataDir,'/corwin/sampling_project'];

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
addParameter(p,      'Clobber',  false,        @islogical);     %overwrite previous output filess (assume no)
addParameter(p,     'Parallel',  false,        @islogical);     %parallel or linear mode? (assumes single-threaded)
addParameter(p,   'TextUpdate',   true,        @islogical);     %display progress of inner sampling loop to screen in linear mode? (yes)
addParameter(p,  'IncludeNaNs',   true,        @islogical);     %if a measurement blob includes NaNs or goes off the edge of the field, do we include it in the result? Useful for e.g. limited-area model runs; if set, output field can include partially-sampled blobs (yes)
addParameter(p,  'GetSettings',  false,        @islogical);     %just get the settings generated here and return it in the 'OldData' field
addParameter(p, 'ScatteredInt',  false,        @islogical);     %use for model data not on a regular lon/lat/prs/t grid. This takes MUCH longer, so don't do it unnecessarily, but can save memory for irregularly-gridded models and may be slightly more accurate for such models.
addParameter(p,  'SaveSingles',   true,        @islogical);     %convert outputs to single() rather than double() to save filespace. Will skip time variables.

%numeric values
addParameter(p,  'ReportEvery',   1000,  IsPositiveInteger);     %if TextUpdate is set, update to screen how often (in samples)?
addParameter(p,      'MinPrs' ,      0,         IsPositive);     %minimum pressure level, i.e. top height. If set to zero, uses model top
addParameter(p,      'MaxPrs' ,   1100,         IsPositive);     %maximum pressure level, i.e. bottom height.
addParameter(p,       'SubSet',      0,  IsPositiveInteger);     %which subfile within the specified day to work on - used e.g. for AIRS granule numbers
addParameter(p,   'HoursAhead',      0,      IsNonNegative);     %if using forecast data, how many hours in advance?

%passed-through data structures
addParameter(p,     'OldData','NOTSET',          @isstruct);    %do we have a previously-used set of model interpolants in memory? 
addParameter(p, 'Sensitivity','NOTSET',          @isstruct);    %parameters for sensitivity-testing mode

%paths
addParameter(p,'DensityPath','./common/saber_density_filled.mat',@isfile);   %path to density data
addParameter(p,    'OutPath',                          'NOTSET', @isfolder); %output file. will be generated automatically below if set to 'NOTSET' or not supplied
 
%arbitrary numbers used in the code. Most defaults have been selected via sensivity testing using 3D AIRS data.
addParameter(p,     'MinSignal',  0.99,        IsPositive);     %fraction of total signal needed to produce final sample
addParameter(p, 'SpecWeightMin',  1e-3,        IsPositive);     %when using specified weighting function, discard any values contributing less than this times the maximum
addParameter(p,     'BlobScale',     3,        IsPositive);     %number of standard deviations to compute sensitivity out to (+- from centre)
addParameter(p,   'MinZContrib',  0.02,        IsPositive);     %when rotating, discard vertical levels contributing less fractional weight than this
addParameter(p,      'ZPadding',   0.5,        IsPositive);     %when rotating, vertical padding in decades of pressure

%parse inputs
%%%%%%%%%%%%%%%%%%%%%%%

parse(p,Instrument,ModelName,DayNumber,varargin{:});
Settings = p.Results; Settings.MasterPath = MasterPath; clear MasterPath

%do some preprocessing based on the above
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ignore the whole routine and just spit out the settings?
if Settings.GetSettings == 1;
  Error = 0;
  OldData = Settings;
  return
end

%is the observational dataset on the valid list?
[~,InstList] = instrument_settings('JUSTLISTING',Settings);
if sum(ismember(InstList,Instrument)) ~= 1;
  disp(['Instrument ''',Instrument,''' not yet supported, stopping'])
  disp( 'Valid (case-sensitive) list is:' )
  for iInst=1:1:numel(InstList)
    disp(['--> ',InstList{iInst}])
  end  
  Error = 5; return
end; 

%is the model on the valid list?
ModelList = model_settings([],'JUSTLISTING');
if sum(ismember(ModelList,ModelName)) ~= 1;
  disp(['Model ''',ModelName,''' not yet supported, stopping'])
  disp( 'Valid (case-sensitive) list is:')
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
if Settings.IncludeNaNs ~=1;  disp('Blobs containing any NaNs will be set to NaN');end

%where shall we put our results? 
if strcmp(Settings.OutPath,'NOTSET')
  Settings.OutPath = [Settings.MasterPath,'/output/',Instrument,'/',ModelName,'/sampled_',num2str(DayNumber),SubSetOutString,ForecastOutString,'.mat'];
  %make sure this directory exists!
  if Settings.SensTestMode ~= 1 && exist([Settings.MasterPath,'/output/',Instrument,'/',ModelName],'dir') ~= 7;
    mkdir([Settings.MasterPath,'/output/',Instrument,'/',ModelName]); end
end

%check if we've already done this day
if exist(Settings.OutPath,'file') ~= 0 && Settings.Clobber == 0
  disp([ModelName,' as ',Instrument,' on ',datestr(DayNumber),' already exists and Clobber is not set, stopping']);
  Error = 3; return;
end

%finally, store the time we started the analysis (especially useful for sensitivity testing)
RunTime = struct();
RunTime.Start = now;

%tidy
clear IsPositive IsPositiveInteger IsNonNegative p varargin InstList ModelList


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load observation properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get instrument options
%this includes the FineGrid spacings needed below, and the SubsetInStrings needed to load multiple-file-per-day datasets
Settings = instrument_settings(Instrument,Settings);

%now identify the file containing the daily grid
ObsGrid = [Settings.MasterPath,'/tracks/',Instrument,'/track_',Instrument,'_',num2str(DayNumber),Settings.SubSetInString,'.mat'];

%and load it
try; ObsGrid = load(ObsGrid);
catch;  %failed.
  disp(['Cannot locate file containing observational track: ',ObsGrid]);
  Error = 1; return
end

%convert obs pressure grid to log space, to make spacings more regular
ObsGrid.Track.Prs = log10(ObsGrid.Track.Prs);

%for backwards compatability: if we don't have a weight format, use 'gaussian'.

if ~isfield(ObsGrid.Weight,'Format'); ObsGrid.Weight.Format = repmat({'gaussian'},size(ObsGrid.Track.Lat)); end

%if any of the instruments we're sampling have a specified weight field, load it now to avoid duplication in-loop
if find(contains(ObsGrid.Weight.Format,'specified_2d'))
  Ws = unique(ObsGrid.Weight.Format);
  for iW=1:1:numel(Ws)
    W = Ws{iW};
    if length(W) > 12 && strcmp(W(1:12),'specified_2d');
      ObsGrid.WeightMatrix.XZ.(W(14:end))     = rCDF(ObsGrid.Params.(W(14:end)).WeightDetails.File);
      ObsGrid.WeightMatrix.XZ.(W(14:end)).avk = smoothn(ObsGrid.WeightMatrix.XZ.(W(14:end)).avk,[1,1,5]);
      ObsGrid.WeightMatrix.Y.(W(14:end))      =      ObsGrid.Params.(W(14:end)).WeightDetails.Y;
    end
  end
  clear Ws iW W
end

  
if Settings.SubSet ~= 0; ExtraInfo = [' subset ',num2str(Settings.SubSet)]; else; ExtraInfo = ''; end
disp([Instrument,' track loaded for ',datestr(DayNumber),ExtraInfo]);
clear ExtraInfo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and prep model data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%identify and load the model
%file conventions vary massively, so load the data inside the child function
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

  %ok, what we do next depends on whether the model is on a lon/lat/prs/t grid (better)
  %or a set of scattered points (much slower)
  if Settings.ScatteredInt ~= 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %gridded data. Phew.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Using gridded interpolation to generate model fields')

    %drop parts of the atmosphere that we don't care about
    InPrsRange = find(Model.Prs <= log10(Settings.MaxPrs) & Model.Prs >= log10(Settings.MinPrs));
    Model.Prs = Model.Prs(InPrsRange);
    Model.T   = Model.T(:,:,:,InPrsRange);
    clear InPrsRange
    
    %limit the analysis to only use height levels present in the model data
    if 10.^max(Model.Prs) < Settings.MaxPrs
      Settings.MaxPrsSet = Settings.MaxPrs;
      Settings.MaxPrs    = round(10.^max(Model.Prs).*10000)./10000;
      disp(['Changing Settings.MaxPrs to model data maximum value: now ',num2str(Settings.MaxPrs),' hPa. Original value retained as Settings.MaxPrsSet'])
    end
    if 10.^min(Model.Prs) > Settings.MinPrs
      Settings.MinPrsSet = Settings.MinPrs;
      Settings.MinPrs    = round(10.^min(Model.Prs).*10000)./10000;
      disp(['Changing Settings.MinPrs to model data minimum value: now ',num2str(Settings.MinPrs),' hPa. Original value retained as Settings.MinPrsSet'])
    end

    %make sure all variables increase montonically
    Fields = {'Time','Lon','Lat','Prs'}; %order is important - this is the same as the Model.T array
    for iField=1:1:numel(Fields)

      %get the order, and re-order the index variable
      [Model.(Fields{iField}),idx] = sort(Model.(Fields{iField}));

      %re-order the main variable along this dimension
      dims  = 1:1:4;
      dims2 = unique([iField,dims],'stable');
      Model.T = permute(Model.T,dims2);
      Model.T = Model.T(idx,:,:,:);
      Model.T = permute(Model.T,[dims(dims < iField)+1,1,dims(dims > iField)]);


    end
    clear Fields iField idx dims dims2


    %make sure pressure increases in z. needed for griddedinterpolant below.
    if mean(diff(Model.Prs)) < 0
      Model.Prs = flip(Model.Prs,1);
      Model.T   = flip(Model.T,  4);
    end

    %produce gridded interpolant fields for each timestep
    [x,y,z] = ndgrid(Model.Lon,Model.Lat,Model.Prs);
    Interpolants = struct();
    for iTime=1:1:numel(Model.Time)
      F = griddedInterpolant(single(x),single(y),single(z),single(squeeze(Model.T(iTime,:,:,:))));
      Interpolants.(['t',num2str(iTime)]) = F;
    end
    Interpolants.Time = Model.Time;
    clear x y z iTime F

  else
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %scattered data. Same code complexity, 
    %but it gonna be slow.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    disp('Using scattered interpolation to generate model fields')

    %produce interpolant for each timestep
    %scatteredinterpolant only takes doubles, which is frustrating...
    for iTime=1:1:numel(Model.Time)
      F = scatteredInterpolant(double(Model.LonI(:)),double(Model.LatI(:)),double(Model.PrsI(:)), ...
                               double(squeeze(Model.TI(iTime,:)))');
      Interpolants.(['t',num2str(iTime)]) = F;
    end
    Interpolants.Time = Model.Time;
  end

  %store the model ID, time period, pressure and interpolants for future use
  OldData.I = Interpolants;
  OldData.ModelPressureScale = Model.Prs;
  OldData.ModelTimeScale = Model.Time;
  OldData.ModelID = ModelName;
  OldData.Time = [min(Model.Time),max(Model.Time)];

  %tell the user we're done
  disp([ModelName,' interpolants prepared for sampling on ',datestr(DayNumber)]);

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
 
  Settings.OutPath = [Settings.MasterPath,Sensitivity.NewPath]; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now go down to the datapoint level and start sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%final steps before core loop:

%1. if we're rotating in the vertical plane, we need density data.
if sum(ObsGrid.Track.ViewAngleZ,'omitnan') ~= 0
  %load density data (global-mean daily-mean saber-derived density)
  Density = load(Settings.DensityPath);
  Density = Density.Results;
  Density.Prs = log10(h2p(Density.HeightScale./1000));
else
  Density = []; %not used, don't waste time making it
end

%2. create a results array and an array for the naive approximation where we just sample the model at measurement centre
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
  %parfor loop. Usually much faster, but needs more memory so won't work for large datasets.
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
[~,Sampled_Data] = unified_reconstructor(Settings,ObsGrid,FinalT,SimpleT);

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
  %set some placeholder values for if the function need to exit
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

  %output handling
  InnerError = 1; %assume an error unless the function runs successfully
  TSample = NaN;
  TSimple = NaN;

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
  %get geolocation and angular parameters
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Sample = struct();
  Sample.Lon     = ObsGrid.Track.Lon(       iSample);
  Sample.Lat     = ObsGrid.Track.Lat(       iSample);
  Sample.Prs     = ObsGrid.Track.Prs(       iSample); %remember, we already log10-ed this
  Sample.Time    = ObsGrid.Track.Time(      iSample);
  
  Sample.AngleH  = ObsGrid.Track.ViewAngleH(iSample);
  Sample.AngleZ  = ObsGrid.Track.ViewAngleZ(iSample);

  %safety-check that none of the above are NaN
  if isnan(Sample.Lon     + Sample.Lat     + Sample.Prs + Sample.Time ...
         + Sample.AngleH  + Sample.AngleZ                             ...   
          ); InnerError = 2;return; 
  end
   
  %check we're not over the top of the model
  if Sample.Prs < log10(Settings.MinPrs); InnerError = 3;return; end

  %or below the bottom 
  if Sample.Prs > log10(Settings.MaxPrs); InnerError = 4;return; end  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %compute the weighting parameters. This varies depending on the format used.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  WeightType = ObsGrid.Weight.Format{iSample};

  if strcmpi(WeightType,'gaussian');

    %'classic' set of 1D gaussians for each dimension. Easy.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %get the weight values
    Sample.WeightX = ObsGrid.Weight.X(iSample);
    Sample.WeightY = ObsGrid.Weight.Y(iSample);
    Sample.WeightZ = ObsGrid.Weight.Z(iSample);

    %safety-check that none of the above are NaN
    if isnan(Sample.WeightX + Sample.WeightY + Sample.WeightZ); InnerError = 2;return;  end

    %compute the fine grid and weights thereon
    [Fine,W] = gaussian_blob(Sample,Settings);

  elseif strcmp(WeightType(1:12),'specified_2d')

    %a specified 2D field in the aLOS and z direction,
    %with a cross-track Gaussian rolloff
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %identify kernel to use
    Instrument = ObsGrid.Weight.Format{iSample}; Instrument = Instrument(14:end);

    %compute the fine grid and weights thereon
    [Fine,W] = specified_2d(Sample,ObsGrid,Instrument,Settings);

  else
    error('Invalid format for weighting functions, terminating')
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %cut data down to reduce interpolation requirements
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  %we want to choose the minimum number of points that will give us our 
  %required weight.
  %this is a local time bottleneck in our routine, but is saving time overall
  %because it reduces the number of points computed in all subsequent steps
  
  %first, sort the weight points
  [WSort,Widx]  = sort(abs(W(:)),'descend');
  %find the cumulative sum
  WSort = cumsum(abs(WSort));
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
    
    
    %Z rotations are more complex than horizontal-plane ones.
    %we have to:
    %1. shift to centre on z=0, to avoid awkward translational terms in the rotation
    %2. carry out the rotation
    %3. undo the shift
    %4. correct for density
    %angles are defined a/c/w from nadir
    
    %extract correct date's density. After the end of the SABER density file it will just use the LAST date present there.
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

    %this has messed up the absolute, but not relative weights
    %so we need to renormalise the distribution
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
  
  %ok, we're there! choose the right interpolant object for the chosen time
  [~,idx] = min(abs(Sample.Time-Interpolants.Time));
  I = Interpolants.(['t',num2str(idx)]);
  
  %and apply it
  if Settings.ScatteredInt ~= 1
    %singles are easier on our memory bottlenecks...
    Fine.T = I(single(Fine.Lon(:)),single(Fine.Lat(:)),single(Fine.Prs(:)));
  else
    %...but the scatteredInterpolant function can't cope with them
    Fine.T = I(double(Fine.Lon(:)),double(Fine.Lat(:)),double(Fine.Prs(:)));
  end
  
  %scale by weights, and store
  if Settings.IncludeNaNs == 1;
    TSample = sum(flatten(Fine.T).*flatten(W),'omitnan');
  else
    TSample = sum(flatten(Fine.T).*flatten(W));
  end
  
  %also compute simple T at measurement centre, for comparison purposes
  if Settings.ScatteredInt ~= 1; TSimple = I(single(Sample.Lon),single(Sample.Lat),single(Sample.Prs));
  else                           TSimple = I(double(Sample.Lon),double(Sample.Lat),double(Sample.Prs));
  end

  

  %success!
  InnerError = 0;

  end %this is the end of the inner (sampling core) function
