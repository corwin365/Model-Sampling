% % function [Error,OldData] = sampling_core_singlethread(Instrument,ModelType,DayNumber,OldData,NoClobber,Sensitivity,SubSet,ForecastHours)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sample models and reanalyses with satellite scan tracks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%possible error states:
%0. successful
%1. problem loading daily obs file
%2. problem loading model file
%3. already done, and noclobber set
%4. sensitivity testing, and output path not set


  %testing parameters
  clearvars -except OldData
  Instrument      = 'airs3d';
  DayNumber       = datenum(2010,10,17);
  ModelType       = 'cesm_ck';
  SubSet = 183;
% % 
% %   %sensitivity-testing testing parameters
% %   Sensitivity.Mode = 1;
% %   Sensitivity.NewPath = 'test.mat';
% %   Sensitivity.FineGrid.X = 20;
% %   Sensitivity.FineGrid.Y = 20;
% %   Sensitivity.FineGrid.Prs = 1./20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% general preliminaries
%Many options here may be overwritten later if we're in sensitivity-testing mode.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display what we think we're doing, for the log file
disp(['This routine is set to sample ',upper(ModelType),' as ',upper(Instrument),' on ',datestr(DayNumber),'.']);


%get main functions
addpath(genpath('../common/'));
addpath(genpath('/home/f/cw785/Matlab/20190122ModelSamplingCode/')); %duplicates above line - specific problem on Balena

%clarify optional inputs
if exist('NoClobber')     ~=1; NoClobber =0 ;       end %handles whether we repeat files done before
if exist('Sensitivity')   ~=1; Sensitivity.Mode =0 ;end %handles whether we're in sensitivity-testing mode
if numel(Sensitivity)     ==0; Sensitivity.Mode =0 ;end %allows us to set var as [] to skip

%very large datsets may need splitting into subdirunal subsets
%the original case this was written for was AIRS-3D, which is in 240
%granules per day. Setting a value here allows the sampling routine to address
%these subsets individually.
if ~exist('SubSet')
  SubSet = 0;
end
if SubSet == 0
  SubSetOutString = '';
else
  if SubSet == 0; SubSetOutString = '';
  else SubSetOutString = ['_subset_',sprintf('%06d',SubSet)];
  end
  disp(['Processing subset ',sprintf('%06d',SubSet)])
end

%get core variables. 
CoreVars = sampling_core_variables;

%declare settings struct - needed for parallelisation
Settings = struct();


%if we're using forecast data, we need to specify how far ahead a forecast
%we want. This will be used by the model selection subroutine to select
%what data to feed back to the sampling parent
if exist('ForecastHours') ~= 0
  if ForecastHours ~= 0;
    Settings.HoursAhead = ForecastHours;
    clear ForecastHours
     disp(['Using ',num2str(Settings.HoursAhead),' hour forecast'])
    ForecastOutString = ['_',sprintf('%03d',Settings.HoursAhead),'hrs'];
  end
end
if ~isfield(Settings,'HoursAhead'); 
  ForecastOutString = '';
  Settings.HoursAhead = 0;
end

%do we have a previously-used set of model interpolants in memory? if not, create an empty checking variable
if ~exist('OldData'); OldData.ModelID = ''; end
if ~isfield(OldData,'ModelID'); OldData.ModelID = ''; end

%we're only interested in the stratosphere, so let's drop the lower troposphere
Settings.MaxPrs = 200;% hPa

%we only want to spend time interpolating onto points that meaningfully
%contribute. so let's define a percentage of the signal we must recover
Settings.MinSignal = 0.99; %fraction

%our blob sizes are specified in st devs. how many times this do we want to include?
%this is +- from centre
Settings.BlobScale = 3;

%where is the density data stored? used for measurements which we rotate in the vertical
DensityPath = '../common/saber_density_filled.mat'; %this version has NaNs filled with the all-time mean for that height, and is linear extrapolated into the t'sphere

%where shall we put our results? 
OutPath = [CoreVars.MasterPath,'/samples/',upper(Instrument),'/',upper(ModelType),'/sampled_',num2str(DayNumber),SubSetOutString,ForecastOutString,'.mat'];

%check if we've already done this day
if exist(OutPath) ~= 0 && NoClobber ~= 0
  Error = 3; 
  disp([upper(ModelType),' as ',upper(Instrument),' on ',datestr(DayNumber),' already exists and NoClobber is set, so skipping']);
  return;
end
clear NoClobber


%finally, store the time we started the analysis (especially useful for sensitivity testing!)
RunTime = struct();
RunTime.Start = now;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load observation properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use try and catch to cover all error type bases in one
try

  %get instrument options
  [Settings,SubSetOutString] = instrument_settings(Instrument,Settings,SubSet);
  
  %now identify the file containing the daily grid
  ObsGrid = [Settings.ObsProperties.Path,'/track_', ...
             Settings.ObsProperties.FileString,'_',num2str(DayNumber),SubSetOutString,'.mat'];
  %and load it
  ObsGrid = load(ObsGrid);

  %convert obs grid to log space, to make spacings more regular
  ObsGrid.Track.Prs = log10(ObsGrid.Track.Prs);
  
catch
  %one of the above failed. escape.
  disp('Problem obtaining observational track.')
  Error = 1;
  return
end


disp([upper(Instrument),' track loaded for ',datestr(DayNumber)]);

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
%this routine then produces a 3dinterpolant field for each timestep based on this


%firstly, do we need to bother, or did we already load it?
if strcmp(ModelType,OldData.ModelID) == 1 ...
       && DayNumber >= min(OldData.Time)  ...
       && DayNumber <  max(OldData.Time) %this needs to be < rather than <= as we may have the first timestep of the day in memory from the last day
   
  %same day as last loop - just call the data back up from memory
  Interpolants = OldData.I;
  disp(['Using previously loaded data and interpolants for ',upper(ModelType),' on ',datestr(DayNumber)]);
                 
else
     

  %clear any previous model data we had stored - saves memory
  clear OldData  
  
  %load the model
  [Model,OldData,Settings] = model_settings(DayNumber,ModelType,Settings,ObsGrid);

  %check if we loaded it successfully
  if Model.Error ~=0
    OldData.ModelID = ''; %we need to set this for the return
    disp('Problem loading model data. Stopping')
    Error =2;
    return
  end

  %convert pressure to log-prs, to make life easier
  Model.Prs = log10(Model.Prs);
  
  
  %finally, drop the lower troposphere - we don't care here, and it takes loads of
  %memory to deal with
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

  disp([upper(ModelType),' interpolants prepared for sampling on ',datestr(DayNumber)]);
  
  %finally, store the model ID, time period, pressure and interpolants for future use
  OldData.I = Interpolants;
  OldData.ModelPressureScale = Model.Prs;
  OldData.ModelTimeScale = Model.Time;  
  OldData.ModelID = ModelType;  
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


if Sensitivity.Mode ~= 0
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %sensitivity-testing mode is active. change grid parameters and output path
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %print out what settings are in force
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  disp('**********************************************') 
  disp('***********Sensitivity testing mode***********')  
  Sensitivity 
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
    if isfield(Sensitivity.FineGrid,'X');   Settings.FineGrid.X   = Sensitivity.FineGrid.X;   end
    if isfield(Sensitivity.FineGrid,'Y');   Settings.FineGrid.Y   = Sensitivity.FineGrid.Y;   end
    if isfield(Sensitivity.FineGrid,'Prs'); Settings.FineGrid.Prs = Sensitivity.FineGrid.Prs; end    
  end
  
  %then, change the output path
  %this is required as otherwise it'll overwrite real results
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if ~isfield(Sensitivity,'NewPath')
    Error = 3;
    disp('Sensitivity test output file not set')
    return;
  end
 
  OutPath = [CoreVars.MasterPath,Sensitivity.NewPath];

  
  DensityPath = Sensitivity.DensityPath; %we might be running it from a different directory

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now we need to go down to the profile level and start sampling
%% everything in this part needs to be parallel-safe, for HPC bulk processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%final steps before core loop:

%1. if we're rotating in the vertical plane, we need density data.
if nansum(ObsGrid.Track.ViewAngleZ) ~= 0
  %load density data (global-mean daily-mean saber-derived density)
  Density = struct();
  Density = load(DensityPath);
  Density = Density.Results;
  Density.Prs = log10(h2p(Density.HeightScale./1000));
else
  Density = []; %var needed, otherwise parfor has trouble classifying it. not used.
end

%2. we need to create a results array
FinalT = ObsGrid.Track.Lon.*NaN; 

%3. also create an array for the 'naive' approximation where we just sample the model at measurement centre
SimpleT = FinalT;

%done! let's go!
disp(['Commencing single-threaded sampling, ',num2str(numel(FinalT)),' samples to extract'])
for iSample = 1:1:numel(FinalT) %replace for with parfor to split over multiple cores
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %this section prints progress to the screen - comment out for non-testing use
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  if iSample == 1; tic; end
% %  if mod(iSample,10000) == 0; disp([round(iSample./numel(FinalT).*100.*100)./100,toc]);tic; end 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %these vars will generate warnings if not declared
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  Channel = [];
  Wx      = [];   
  Wy      = [];
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
          ); continue; end
   
  %check we're not over the top of the model
  if Sample.Prs < min(OldData.ModelPressureScale); continue; end
  %or below our imposed bottom 
  if Sample.Prs > log10(Settings.MaxPrs); continue; end  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %produce the fine sampling grid
  %this is the grid the actual sampling will be done
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %declare a struct to keep the grids in
  Fine = struct();
  
  %creating the grid in the horizontal is easy
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Fine.X = -Settings.BlobScale.*Sample.WeightX : Settings.FineGrid.X : Settings.BlobScale.*Sample.WeightX;
  Fine.Y = -Settings.BlobScale.*Sample.WeightY : Settings.FineGrid.Y : Settings.BlobScale.*Sample.WeightY;
  
  %for the vertical, we're working in pressure co-ordinates, so this is fiddly
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %most instruments are a simple Gaussian in height
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
    
    %load it
    Channel.Prs = log10(ObsGrid.Weight.ZFuncs.PrsScale);
    Channel.W   = ObsGrid.Weight.ZFuncs.Weights(Channel.ID,:); %keep for later  

    %discard parts that contribute very little
    Channel.W(Channel.W < 0.02.*nansum(Channel.W(:))) = 0; %less than 2% is a good balance of useful volume and runtime (tested for AIRS 42km)

    
    %set NaNs to zero
    Channel.W(isnan(Channel.W)) = 0;
    
    %remove low altitudes
    GoodPrs = find(Channel.Prs < log10(Settings.MaxPrs));
    Channel.W   = Channel.W(   GoodPrs);
    Channel.Prs = Channel.Prs(GoodPrs);
    
    %find highest and lowest prs level remaining
    z = Channel.Prs(find(Channel.W > 0));
    z = [max(z),min(z)]+[1,-1].*0.5; %add half a decade each side padding
  end
  Fine.Prs = single(z(2):Settings.FineGrid.Prs:z(1));


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
    stop %not written yet, as no such cases implemented (would be a weird instrument)
  end  
  
  %1d across-track weight
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if Sample.WeightY > 0  %simple Gaussian approximation
    Wy = normpdf(Fine.Y,0,Sample.WeightY);
  else %not a simple Gaussian - produce from supplied information
    stop %not written yet, as no such cases implemented (would be a weird instrument)
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
  [~,Cidx] = min(abs(WSort-max(WSort).*Settings.MinSignal ));
  %and keep only these points
  Points = Widx(1:Cidx);
 
  %pull out these points out of the weight and location arrays
  W        = W(Points);
  Fine.X   = Fine.X(Points);
  Fine.Y   = Fine.Y(Points);  
  Fine.Prs = Fine.Prs(Points); 

  %normalise remaining weights to sum to 1
  W = W./nansum(W(:));  

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
    W = W./nansum(W(:));
    
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
  FinalT(iSample) = nansum(flatten(Fine.T).*flatten(W));
  
  %also compute simple T at measurement centre
  SimpleT(iSample) = I(single(Sample.Lon),single(Sample.Lat),single(Sample.Prs));

end
disp(['Sampling complete for ',upper(ModelType),' as ',upper(Instrument),' on ',datestr(DayNumber)]);



%success! save and return

%create results struct from original geolocation
Output = ObsGrid.Track;

%put our sampled T into the struct
Output.T = FinalT;
Output.TSimple = SimpleT;

%remove some fields to save filespace
Output = rmfield(Output,{'ViewAngleH','ViewAngleZ'});

%save reconstruction data
Recon = ObsGrid.Recon; 

%and runtime (especially useful for sensitivity testing!)
RunTime.End = datenum(now);

%write the file
save(OutPath,'Output','Settings','Recon','RunTime');
disp('Saved!');

%and we're done!
Error = 0; %we didn't fail!
return
